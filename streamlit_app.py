import subprocess
import tempfile
from pathlib import Path

import streamlit as st

st.set_page_config(page_title="FASTA+BED+Guides → GenBank", layout="centered")
st.title("Annotate haplotypes: BED + guides_even/odd → combined GFF3 → GenBank")

st.markdown(
    """
Upload:
- **custom_reference.fa** (FASTA)
- **custom_regions.bed** (BED4)
- **guides_even.fa** (FASTA)
- **guides_odd.fa** (FASTA)

This app runs:
1) BED → GFF3  
2) Bowtie exact mapping for even + odd guides  
3) Between-regions annotation (downstream only, ≤20 kb, excludes regions containing N)  
4) Combine GFF3s  
5) GFF3 + FASTA → GenBank
"""
)

ref_fa = st.file_uploader("Reference FASTA (custom_reference.fa)", type=["fa", "fasta", "fna"])
regions_bed = st.file_uploader("Regions BED4 (custom_regions.bed)", type=["bed", "txt"])
guides_even = st.file_uploader("Guides EVEN FASTA (guides_even.fa)", type=["fa", "fasta", "fna"])
guides_odd = st.file_uploader("Guides ODD FASTA (guides_odd.fa)", type=["fa", "fasta", "fna"])

with st.expander("Advanced", expanded=False):
    max_between_len = st.number_input("Max between-region length (bp)", min_value=1, value=20000, step=1000)
    keep_intermediates = st.checkbox("Show & allow download of intermediate files", value=True)

run_btn = st.button(
    "Run",
    type="primary",
    disabled=not (ref_fa and regions_bed and guides_even and guides_odd),
)

def run_cmd(cmd, cwd: Path):
    """Run command, return (rc, stdout, stderr)."""
    p = subprocess.run(cmd, cwd=str(cwd), text=True, capture_output=True)
    return p.returncode, p.stdout, p.stderr

def fail(step, stdout, stderr):
    st.error(f"Failed at: {step}")
    if stdout:
        st.subheader("stdout")
        st.code(stdout)
    if stderr:
        st.subheader("stderr")
        st.code(stderr)
    st.stop()

if run_btn:
    with tempfile.TemporaryDirectory() as td:
        wd = Path(td)

        # Write uploads to temp directory
        ref_path = wd / "custom_reference.fa"
        bed_path = wd / "custom_regions.bed"
        even_path = wd / "guides_even.fa"
        odd_path = wd / "guides_odd.fa"

        ref_path.write_bytes(ref_fa.getvalue())
        bed_path.write_bytes(regions_bed.getvalue())
        even_path.write_bytes(guides_even.getvalue())
        odd_path.write_bytes(guides_odd.getvalue())

        # Outputs / intermediates
        regions_gff = wd / "custom_regions.gff3"
        even_mapped = wd / "guides_even.mapped"
        odd_mapped = wd / "guides_odd.mapped"
        even_gff = wd / "guides_even.gff3"
        odd_gff = wd / "guides_odd.gff3"
        between_even = wd / "between_regions_even.gff3"
        between_odd = wd / "between_regions_odd.gff3"
        combined_gff = wd / "combined.gff3"
        gbk_path = wd / "custom.gbk"

        log_area = st.empty()
        logs = []

        def log(msg):
            logs.append(msg)
            log_area.code("\n".join(logs))

        st.info("Running pipeline…")
        log("Starting...")

        # 1) BED -> GFF3
        log("1) BED -> GFF3")
        awk_script = r"""BEGIN{OFS="\t"; print "##gff-version 3"}
{
  chrom=$1; start=$2+1; end=$3; name=$4;
  if (name=="" || name==".") name="feature_"NR;
  id=name"_"NR;
  print chrom,"bed","misc_feature",start,end,".",".",".","ID="id";Name="name";label="name
}"""
        rc, out, err = run_cmd(["awk", awk_script, str(bed_path)], wd)
        if rc != 0:
            fail("BED -> GFF3 (awk)", out, err)
        regions_gff.write_text(out, encoding="utf-8")

        # 2) bowtie-build
        log("2) bowtie-build")
        rc, out, err = run_cmd(["bowtie-build", str(ref_path), "refidx"], wd)
        if rc != 0:
            fail("bowtie-build", out, err)

        # Helper to map one guide file and produce a guides.gff3
        def map_guides_to_gff(guides_fa_path: Path, mapped_out: Path, gff_out: Path, label: str):
            log(f"3) bowtie map {label}")
            rc, out, err = run_cmd(
                ["bowtie", "-f", "-v", "0", "-a", "--best", "--strata", "refidx", str(guides_fa_path)],
                wd,
            )
            if rc != 0:
                fail(f"bowtie map ({label})", out, err)
            mapped_out.write_text(out, encoding="utf-8")

            log(f"4) mapped -> GFF3 {label}")
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

        # 5) between-regions
        log("5) between-regions even")
        rc, out, err = run_cmd(
            ["python", "make_between_regions.py",
             "--gff", str(even_gff),
             "--fasta", str(ref_path),
             "--out", str(between_even),
             "--maxlen", str(int(max_between_len))],
            wd,
        )
        if rc != 0:
            fail("make_between_regions.py (even)", out, err)
        log(out.strip() or "between even done")

        log("6) between-regions odd")
        rc, out, err = run_cmd(
            ["python", "make_between_regions.py",
             "--gff", str(odd_gff),
             "--fasta", str(ref_path),
             "--out", str(between_odd),
             "--maxlen", str(int(max_between_len))],
            wd,
        )
        if rc != 0:
            fail("make_between_regions.py (odd)", out, err)
        log(out.strip() or "between odd done")

        # 6) combine all gffs (avoid multiple headers)
        log("7) combine GFF3s")
        def strip_header(text: str) -> str:
            return "\n".join([ln for ln in text.splitlines() if not ln.startswith("##gff-version") and ln.strip() != ""])

        combined = ["##gff-version 3"]
        combined.append(strip_header(regions_gff.read_text(encoding="utf-8")))
        combined.append(strip_header(even_gff.read_text(encoding="utf-8")))
        combined.append(strip_header(odd_gff.read_text(encoding="utf-8")))
        combined.append(strip_header(between_even.read_text(encoding="utf-8")))
        combined.append(strip_header(between_odd.read_text(encoding="utf-8")))
        combined_gff.write_text("\n".join([c for c in combined if c != ""]) + "\n", encoding="utf-8")

        # 7) gff-to-genbank
        log("8) GFF3+FASTA -> GenBank")
        rc, out, err = run_cmd(["gff-to-genbank", str(combined_gff), str(ref_path)], wd)
        if rc != 0:
            fail("gff-to-genbank", out, err)
        gbk_path.write_text(out, encoding="utf-8")

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
            st.download_button("custom_regions.gff3", regions_gff.read_bytes(), file_name="custom_regions.gff3")
            st.download_button("guides_even.gff3", even_gff.read_bytes(), file_name="guides_even.gff3")
            st.download_button("guides_odd.gff3", odd_gff.read_bytes(), file_name="guides_odd.gff3")
            st.download_button("between_regions_even.gff3", between_even.read_bytes(), file_name="between_regions_even.gff3")
            st.download_button("between_regions_odd.gff3", between_odd.read_bytes(), file_name="between_regions_odd.gff3")
            st.download_button("guides_even.mapped", even_mapped.read_bytes(), file_name="guides_even.mapped")
            st.download_button("guides_odd.mapped", odd_mapped.read_bytes(), file_name="guides_odd.mapped")
