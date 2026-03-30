import json
import re
import subprocess
import tempfile
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import streamlit as st

from gff_to_genbank_patched import gff_fasta_to_genbank

APP_DIR = Path(__file__).resolve().parent
BETWEEN_SCRIPT = APP_DIR / "make_between_regions.py"

st.set_page_config(page_title="FASTA+BED+Guides -> GenBank", layout="centered")
st.title("Annotate haplotypes: BED + guides_even/odd -> combined GFF3 -> GenBank")

st.markdown(
    """
Choose one input mode:

- **Use uploaded custom reference**
- **Fetch hg38 locus by coordinates**

The app then runs:
1) Reference preparation
2) Region annotations -> GFF3
3) Bowtie exact mapping for even + odd guides
4) Between-regions annotation (downstream only, <=20 kb, excludes regions containing N)
5) Combine GFF3s
6) GFF3 + FASTA -> GenBank
"""
)

# -----------------------------
# Helpers
# -----------------------------

COORD_RE = re.compile(r"^(?:chr)?([A-Za-z0-9_]+):(\d+)-(\d+)$")

def run_cmd(cmd, cwd: Path):
    """Run command, return (rc, stdout, stderr)."""
    p = subprocess.run(cmd, cwd=str(cwd), text=True, capture_output=True)
    return p.returncode, p.stdout, p.stderr

def fail(step, stdout="", stderr=""):
    st.error(f"Failed at: {step}")
    if stdout:
        st.subheader("stdout")
        st.code(stdout)
    if stderr:
        st.subheader("stderr")
        st.code(stderr)
    st.stop()

def stop_with_error(msg: str):
    st.error(msg)
    st.stop()

def http_get_text(url: str, headers: dict[str, str] | None = None) -> str:
    req = Request(url, headers=headers or {})
    try:
        with urlopen(req, timeout=60) as resp:
            return resp.read().decode("utf-8")
    except HTTPError as e:
        body = ""
        try:
            body = e.read().decode("utf-8", errors="replace")
        except Exception:
            pass
        stop_with_error(f"HTTP error while fetching remote data: {e.code}\n{body}")
    except URLError as e:
        stop_with_error(f"Network error while fetching remote data: {e}")

def http_get_json(url: str, headers: dict[str, str] | None = None):
    txt = http_get_text(url, headers=headers)
    try:
        return json.loads(txt)
    except json.JSONDecodeError as e:
        stop_with_error(f"Failed to parse JSON response from remote service: {e}\n\nResponse:\n{txt[:4000]}")

def parse_locus_or_stop(coord_text: str) -> tuple[str, int, int]:
    """
    Accept:
      chr7:55019017-55211628
      7:55019017-55211628

    Returns:
      chrom_no_chr, start1, end1
    """
    s = (coord_text or "").strip()
    m = COORD_RE.match(s)
    if not m:
        stop_with_error(
            "Invalid locus coordinates. Please use format like chr7:55019017-55211628 or 7:55019017-55211628."
        )

    chrom, start_s, end_s = m.groups()
    start1 = int(start_s)
    end1 = int(end_s)

    if start1 < 1 or end1 < 1 or end1 < start1:
        stop_with_error(
            "Invalid locus coordinates. Start and end must be positive integers, and end must be >= start."
        )

    # strip 'chr' by regex design above
    return chrom, start1, end1

def make_local_seqid(chrom_no_chr: str, start1: int, end1: int) -> str:
    return f"hg38_{chrom_no_chr}_{start1}_{end1}"

def wrap_fasta(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def fetch_hg38_sequence_ucsc(chrom_no_chr: str, start1: int, end1: int) -> str:
    """
    User input is 1-based closed.
    UCSC sequence API expects:
      chrom = chrN
      start = 0-based
      end   = 1-based
    """
    chrom_ucsc = f"chr{chrom_no_chr}"
    start0 = start1 - 1
    url = (
        "https://api.genome.ucsc.edu/getData/sequence"
        f"?genome=hg38;chrom={chrom_ucsc};start={start0};end={end1}"
    )
    data = http_get_json(url, headers={"Content-Type": "application/json"})

    # UCSC response includes the sequence in the 'dna' field for this endpoint.
    seq = data.get("dna", "")
    if not seq:
        stop_with_error(f"UCSC returned no sequence for hg38 {chrom_ucsc}:{start1}-{end1}.")
    return seq.upper()

def fetch_ensembl_gene_overlaps(chrom_no_chr: str, start1: int, end1: int):
    """
    Use Ensembl overlap endpoint with feature=gene only to keep output uncluttered.
    """
    region = f"{chrom_no_chr}:{start1}-{end1}"
    url = f"https://rest.ensembl.org/overlap/region/human/{region}?feature=gene"
    data = http_get_json(url, headers={"Content-Type": "application/json"})
    if not isinstance(data, list):
        stop_with_error(f"Unexpected Ensembl response for region {region}.")
    return data

def sanitize_attr_value(s: str) -> str:
    return str(s).replace(";", ",").replace("=", "_").replace("\t", " ").strip()

def genes_to_local_gff3(genes, chrom_no_chr: str, start1: int, end1: int, local_seqid: str) -> str:
    """
    Rebase genomic hg38 coordinates onto the extracted local sequence.
    Output only gene-level misc_feature annotations with clean Name/label.
    """
    lines = ["##gff-version 3"]

    for g in genes:
        gstart = g.get("start")
        gend = g.get("end")
        if gstart is None or gend is None:
            continue

        # Clip to requested interval just in case
        clipped_start = max(int(gstart), start1)
        clipped_end = min(int(gend), end1)
        if clipped_start > clipped_end:
            continue

        local_start = clipped_start - start1 + 1
        local_end = clipped_end - start1 + 1

        gene_id = g.get("id") or f"gene_{chrom_no_chr}_{gstart}_{gend}"
        label = (
            g.get("external_name")
            or g.get("gene_name")
            or g.get("description")
            or gene_id
        )
        label = sanitize_attr_value(label)

        strand = g.get("strand")
        if strand == 1:
            strand_char = "+"
        elif strand == -1:
            strand_char = "-"
        else:
            strand_char = "."

        attrs = (
            f"ID={sanitize_attr_value(gene_id)};"
            f"Name={label};"
            f"label={label};"
            f"source=Ensembl;"
            f"genome=hg38;"
            f"orig_coord={chrom_no_chr}:{clipped_start}-{clipped_end}"
        )

        lines.append(
            "\t".join(
                [
                    local_seqid,
                    "Ensembl",
                    "misc_feature",
                    str(local_start),
                    str(local_end),
                    ".",
                    strand_char,
                    ".",
                    attrs,
                ]
            )
        )

    return "\n".join(lines) + "\n"

def bed_to_gff3_text(bed_path: Path) -> str:
    awk_script = r"""BEGIN{OFS="\t"; print "##gff-version 3"}
{
  chrom=$1; start=$2+1; end=$3; name=$4;
  if (name=="" || name==".") name="feature_"NR;
  id=name"_"NR;
  print chrom,"bed","misc_feature",start,end,".",".",".","ID="id";Name="name";label="name
}"""
    rc, out, err = run_cmd(["awk", awk_script, str(bed_path)], bed_path.parent)
    if rc != 0:
        fail("BED -> GFF3 (awk)", out, err)
    return out

def strip_header(text: str) -> str:
    return "\n".join(
        [ln for ln in text.splitlines() if not ln.startswith("##gff-version") and ln.strip() != ""]
    )

# -----------------------------
# UI
# -----------------------------

input_mode = st.radio(
    "Reference source",
    ["Use uploaded custom reference", "Fetch hg38 locus by coordinates"],
)

if input_mode == "Use uploaded custom reference":
    ref_fa = st.file_uploader("Reference FASTA (custom_reference.fa)", type=["fa", "fasta", "fna"])
    regions_bed = st.file_uploader("Regions BED4 (custom_regions.bed)", type=["bed", "txt"])
    locus_coords = None
else:
    locus_coords = st.text_input(
        "Locus coordinates",
        placeholder="chr7:55019017-55211628 or 7:55019017-55211628",
    )
    ref_fa = None
    regions_bed = None

guides_even = st.file_uploader("Guides EVEN FASTA (guides_even.fa)", type=["fa", "fasta", "fna"])
guides_odd = st.file_uploader("Guides ODD FASTA (guides_odd.fa)", type=["fa", "fasta", "fna"])

with st.expander("Advanced", expanded=False):
    max_between_len = st.number_input(
        "Max between-region length (bp)", min_value=1, value=20000, step=1000
    )
    keep_intermediates = st.checkbox("Show & allow download of intermediate files", value=True)

if input_mode == "Use uploaded custom reference":
    ready = bool(ref_fa and regions_bed and guides_even and guides_odd)
else:
    ready = bool(locus_coords and guides_even and guides_odd)

run_btn = st.button("Run", type="primary", disabled=not ready)

# -----------------------------
# Main pipeline
# -----------------------------

if run_btn:
    with tempfile.TemporaryDirectory() as td:
        wd = Path(td)

        # Shared output/intermediate paths
        ref_path = wd / "custom_reference.fa"
        even_path = wd / "guides_even.fa"
        odd_path = wd / "guides_odd.fa"

        regions_gff = wd / "custom_regions.gff3"
        even_mapped = wd / "guides_even.mapped"
        odd_mapped = wd / "guides_odd.mapped"
        even_gff = wd / "guides_even.gff3"
        odd_gff = wd / "guides_odd.gff3"
        between_even = wd / "between_regions_even.gff3"
        between_odd = wd / "between_regions_odd.gff3"
        combined_gff = wd / "combined.gff3"
        gbk_path = wd / "custom.gbk"

        # Only used in custom-upload mode
        bed_path = wd / "custom_regions.bed"

        log_area = st.empty()
        logs = []

        def log(msg):
            logs.append(msg)
            log_area.code("\n".join(logs))

        st.info("Running pipeline...")
        log("Starting...")

        # Write guides
        even_path.write_bytes(guides_even.getvalue())
        odd_path.write_bytes(guides_odd.getvalue())

        fetched_ref_for_download = False
        fetched_regions_for_download = False

        if input_mode == "Use uploaded custom reference":
            # Existing behavior
            log("1) write uploaded reference and BED")
            ref_path.write_bytes(ref_fa.getvalue())
            bed_path.write_bytes(regions_bed.getvalue())

            log("2) BED -> GFF3")
            regions_gff.write_text(bed_to_gff3_text(bed_path), encoding="utf-8")

        else:
            # New hg38 mode
            log("1) parse locus coordinates")
            chrom_no_chr, start1, end1 = parse_locus_or_stop(locus_coords)
            local_seqid = make_local_seqid(chrom_no_chr, start1, end1)

            log(f"2) fetch hg38 sequence for {chrom_no_chr}:{start1}-{end1}")
            seq = fetch_hg38_sequence_ucsc(chrom_no_chr, start1, end1)
            ref_text = f">{local_seqid}\n{wrap_fasta(seq)}\n"
            ref_path.write_text(ref_text, encoding="utf-8")
            fetched_ref_for_download = True

            log("3) fetch Ensembl gene annotations")
            genes = fetch_ensembl_gene_overlaps(chrom_no_chr, start1, end1)

            log("4) convert Ensembl genes -> local GFF3")
            regions_gff.write_text(
                genes_to_local_gff3(genes, chrom_no_chr, start1, end1, local_seqid),
                encoding="utf-8",
            )
            fetched_regions_for_download = True

        # Bowtie index
        log("5) bowtie-build")
        rc, out, err = run_cmd(["bowtie-build", str(ref_path), "refidx"], wd)
        if rc != 0:
            fail("bowtie-build", out, err)

        def map_guides_to_gff(guides_fa_path: Path, mapped_out: Path, gff_out: Path, label: str):
            log(f"6) bowtie map {label}")
            rc, out, err = run_cmd(
                ["bowtie", "-f", "-v", "0", "-a", "--best", "--strata", "refidx", str(guides_fa_path)],
                wd,
            )
            if rc != 0:
                fail(f"bowtie map ({label})", out, err)
            mapped_out.write_text(out, encoding="utf-8")

            log(f"7) mapped -> GFF3 {label}")
            # Bowtie default tabular: qname strand rname off0 seq qual ...
            gawk_script = r"""BEGIN{OFS="\t"; print "##gff-version 3"}
{
  qname=$1; strand=$2; rname=$3; off0=$4; seq=$5;
  if (rname=="*" || off0=="" ) next;

  len=length(seq);
  start1=off0+1;
  end1=off0+len;

  id=qname "." rname "." start1 "." strand;
  attrs="ID=" id ";Name=" qname ";label=" qname ";source=bowtie";
  print rname,"bowtie","misc_feature",start1,end1,".",strand,".",attrs
}"""
            rc, out, err = run_cmd(["gawk", gawk_script, str(mapped_out)], wd)
            if rc != 0:
                fail(f"mapped -> GFF3 ({label})", out, err)
            gff_out.write_text(out, encoding="utf-8")

        map_guides_to_gff(even_path, even_mapped, even_gff, "even")
        map_guides_to_gff(odd_path, odd_mapped, odd_gff, "odd")

        # Between-regions
        log("8) between-regions even")
        rc, out, err = run_cmd(
            [
                "python",
                str(BETWEEN_SCRIPT),
                "--gff",
                str(even_gff),
                "--fasta",
                str(ref_path),
                "--out",
                str(between_even),
                "--maxlen",
                str(int(max_between_len)),
            ],
            wd,
        )
        if rc != 0:
            fail("make_between_regions.py (even)", out, err)
        log(out.strip() or "between even done")

        log("9) between-regions odd")
        rc, out, err = run_cmd(
            [
                "python",
                str(BETWEEN_SCRIPT),
                "--gff",
                str(odd_gff),
                "--fasta",
                str(ref_path),
                "--out",
                str(between_odd),
                "--maxlen",
                str(int(max_between_len)),
            ],
            wd,
        )
        if rc != 0:
            fail("make_between_regions.py (odd)", out, err)
        log(out.strip() or "between odd done")

        # Combine GFF3s
        log("10) combine GFF3s")
        combined = ["##gff-version 3"]
        combined.append(strip_header(regions_gff.read_text(encoding="utf-8")))
        combined.append(strip_header(even_gff.read_text(encoding="utf-8")))
        combined.append(strip_header(odd_gff.read_text(encoding="utf-8")))
        combined.append(strip_header(between_even.read_text(encoding="utf-8")))
        combined.append(strip_header(between_odd.read_text(encoding="utf-8")))
        combined_gff.write_text("\n".join([c for c in combined if c != ""]) + "\n", encoding="utf-8")

        # GenBank conversion
        log("11) GFF3+FASTA -> GenBank (patched)")
        try:
            gff_fasta_to_genbank(str(combined_gff), str(ref_path), str(gbk_path))
        except Exception as e:
            st.error("Failed at: GenBank conversion")
            st.exception(e)
            st.stop()

        st.success("Done!")

        st.download_button(
            "Download GenBank (custom.gbk)",
            data=gbk_path.read_bytes(),
            file_name="custom.gbk",
            mime="application/octet-stream",
        )
        st.download_button(
            "Download combined GFF3 (combined.gff3)",
            data=combined_gff.read_bytes(),
            file_name="combined.gff3",
            mime="text/plain",
        )

        if keep_intermediates:
            st.subheader("Intermediates")

            if input_mode == "Use uploaded custom reference":
                st.download_button(
                    "custom_regions.gff3",
                    regions_gff.read_bytes(),
                    file_name="custom_regions.gff3",
                )
            else:
                st.download_button(
                    "custom_reference.fa",
                    ref_path.read_bytes(),
                    file_name="custom_reference.fa",
                )
                st.download_button(
                    "custom_regions.gff3",
                    regions_gff.read_bytes(),
                    file_name="custom_regions.gff3",
                )

            st.download_button("guides_even.gff3", even_gff.read_bytes(), file_name="guides_even.gff3")
            st.download_button("guides_odd.gff3", odd_gff.read_bytes(), file_name="guides_odd.gff3")
            st.download_button(
                "between_regions_even.gff3",
                between_even.read_bytes(),
                file_name="between_regions_even.gff3",
            )
            st.download_button(
                "between_regions_odd.gff3",
                between_odd.read_bytes(),
                file_name="between_regions_odd.gff3",
            )
            st.download_button("guides_even.mapped", even_mapped.read_bytes(), file_name="guides_even.mapped")
            st.download_button("guides_odd.mapped", odd_mapped.read_bytes(), file_name="guides_odd.mapped")
