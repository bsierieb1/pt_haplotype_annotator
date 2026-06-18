import io
import json
import re
import subprocess
import tempfile
import zipfile
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import streamlit as st

from Bio import SeqIO
from dna_features_viewer import BiopythonTranslator, GraphicFeature, GraphicRecord

from gff_to_genbank_patched import gff_fasta_to_genbank

APP_DIR = Path(__file__).resolve().parent
BETWEEN_SCRIPT = APP_DIR / "make_between_regions.py"

st.set_page_config(page_title="PT annotator", layout="centered")
st.title("PureTarget tile annotator")

st.markdown(
    """
This tool takes a set of gRNAs designed for a PureTarget assay, maps them to a reference, and predicts tiles. Disclaimer: This is not an official PacBio tool.

Choose one input mode:

- **Provide the list of targeted genes or coordinates of targeted loci (human hg38 only)**
- **Upload custom target region FASTA**

The app then runs:
1) Bowtie exact mapping for required guide set and optional odd guide set
2) Tile annotation (downstream only, <=20 kb by default, excludes regions containing N)
3) Combine all annotations
4) Render a locus preview and package downloadable outputs

In human hg38 mode, the app will also check for common SNPs overlapping the guide or its PAM. Guides with non-NGG PAMs will also be flagged.
"""
)

# -----------------------------
# Helpers
# -----------------------------

COORD_RE = re.compile(r"^(?:chr)?([A-Za-z0-9_]+):(\d+)-(\d+)$")
CANONICAL_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "M"}
GENE_FLANK_BP = 20000
UCSC_COMMON_SNP_TRACK = "dbSnp155Common"


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


def canonicalize_chrom(chrom: str) -> str | None:
    chrom = (chrom or "").strip()
    if chrom.lower().startswith("chr"):
        chrom = chrom[3:]
    chrom = chrom.upper()
    if chrom == "MT":
        chrom = "M"
    return chrom if chrom in CANONICAL_CHROMS else None


def parse_locus(coord_text: str) -> tuple[str, int, int]:
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
        raise ValueError(
            f"Invalid locus coordinates '{coord_text}'. Use format like chr7:55019017-55211628 or 7:55019017-55211628."
        )

    chrom, start_s, end_s = m.groups()
    chrom = canonicalize_chrom(chrom)
    if chrom is None:
        raise ValueError(
            f"Invalid locus coordinates '{coord_text}'. Only canonical human chromosomes 1-22, X, Y, and M are supported."
        )

    start1 = int(start_s)
    end1 = int(end_s)

    if start1 < 1 or end1 < 1 or end1 < start1:
        raise ValueError(
            f"Invalid locus coordinates '{coord_text}'. Start and end must be positive integers, and end must be >= start."
        )

    return chrom, start1, end1


def fetch_hg38_gene_by_symbol(symbol: str) -> dict:
    gene_symbol = (symbol or "").strip()
    url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene_symbol}"
    data = http_get_json(url, headers={"Content-Type": "application/json"})
    if not isinstance(data, dict) or not data:
        raise ValueError(f"'{symbol}' could not be resolved as a human gene symbol.")

    seq_region_name = data.get("seq_region_name")
    start1 = data.get("start")
    end1 = data.get("end")
    display_name = data.get("display_name") or data.get("id") or gene_symbol

    chrom = canonicalize_chrom(str(seq_region_name or ""))
    if chrom is None:
        raise ValueError(
            f"Gene symbol '{symbol}' resolved to a non-canonical locus ('{seq_region_name}'), which is not supported."
        )
    if start1 is None or end1 is None:
        raise ValueError(f"Gene symbol '{symbol}' resolved, but coordinates were missing.")

    return {
        "display_name": str(display_name),
        "chrom": chrom,
        "start1": int(start1),
        "end1": int(end1),
    }


def resolve_locus_or_gene(line: str) -> tuple[str, int, int]:
    try:
        return parse_locus(line)
    except ValueError:
        gene = fetch_hg38_gene_by_symbol(line)
        return (
            gene["chrom"],
            max(1, gene["start1"] - GENE_FLANK_BP),
            gene["end1"] + GENE_FLANK_BP,
        )


def merge_overlapping_loci(loci: list[tuple[str, int, int]]) -> list[tuple[str, int, int]]:
    if not loci:
        return []

    merged = []
    for chrom, start1, end1 in sorted(loci, key=lambda x: (x[0], x[1], x[2])):
        if not merged:
            merged.append([chrom, start1, end1])
            continue

        last_chrom, last_start1, last_end1 = merged[-1]
        if chrom == last_chrom and start1 <= last_end1:
            merged[-1][2] = max(last_end1, end1)
        else:
            merged.append([chrom, start1, end1])

    return [(chrom, start1, end1) for chrom, start1, end1 in merged]



def parse_multi_loci_or_stop(multiline_text: str) -> list[tuple[str, int, int]]:
    """
    One interval or gene symbol per line. Deduplicate exact repeated resolved loci, then merge overlapping loci.
    """
    raw_lines = [ln.strip() for ln in (multiline_text or "").splitlines()]
    lines = [ln for ln in raw_lines if ln]

    if not lines:
        stop_with_error(
            "Please enter at least one locus coordinate or gene symbol, one per line."
        )

    parsed = []
    seen = set()

    for i, line in enumerate(lines, start=1):
        try:
            locus = resolve_locus_or_gene(line)
        except ValueError as e:
            stop_with_error(f"Error on line {i}: {e}")
        if locus not in seen:
            seen.add(locus)
            parsed.append(locus)

    return merge_overlapping_loci(parsed)


def make_local_seqid(chrom_no_chr: str, start1: int, end1: int) -> str:
    return f"hg38_{chrom_no_chr}_{start1}_{end1}"


def make_locus_slug(chrom_no_chr: str, start1: int, end1: int) -> str:
    return f"hg38_{chrom_no_chr}_{start1}_{end1}"


def wrap_fasta(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


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


def fetch_ensembl_gene_model(gene_id: str) -> dict:
    """
    Fetch transcripts/exons for a gene. MANE annotations are included when
    Ensembl has them, so previews can prefer a representative isoform.
    """
    url = f"https://rest.ensembl.org/lookup/id/{gene_id}?expand=1;mane=1"
    data = http_get_json(url, headers={"Content-Type": "application/json"})
    if not isinstance(data, dict):
        stop_with_error(f"Unexpected Ensembl transcript response for gene {gene_id}.")
    return data


def strip_stable_id_version(stable_id: str | None) -> str:
    return str(stable_id or "").split(".", 1)[0]


def mane_rank(transcript: dict) -> int:
    mane_entries = transcript.get("MANE") or []
    mane_types = {
        str(entry.get("type") or "")
        for entry in mane_entries
        if isinstance(entry, dict)
    }
    if "MANE_Select" in mane_types:
        return 0
    if mane_entries:
        return 1
    return 2


def select_main_transcript(gene_model: dict) -> dict | None:
    transcripts = [
        tx
        for tx in gene_model.get("Transcript", [])
        if isinstance(tx, dict) and tx.get("Exon")
    ]
    if not transcripts:
        return None

    canonical_id = strip_stable_id_version(gene_model.get("canonical_transcript"))

    def transcript_sort_key(tx: dict):
        tx_id = strip_stable_id_version(tx.get("id"))
        is_canonical = tx_id == canonical_id or bool(tx.get("is_canonical"))
        is_protein_coding = str(tx.get("biotype") or "") == "protein_coding"
        exon_count = len(tx.get("Exon") or [])
        try:
            length = int(tx.get("length") or 0)
        except (TypeError, ValueError):
            length = 0
        return (
            mane_rank(tx),
            0 if is_canonical else 1,
            0 if is_protein_coding else 1,
            -exon_count,
            -length,
            tx_id,
        )

    return sorted(transcripts, key=transcript_sort_key)[0]


def _extract_track_records(obj, track_name: str):
    if isinstance(obj, list):
        return obj
    if isinstance(obj, dict):
        if track_name in obj and isinstance(obj[track_name], list):
            return obj[track_name]
        for value in obj.values():
            if isinstance(value, list) and value and isinstance(value[0], dict):
                return value
    return []


def fetch_ucsc_common_snps(chrom_no_chr: str, start1: int, end1: int):
    chrom_ucsc = f"chr{chrom_no_chr}"
    start0 = start1 - 1
    url = (
        "https://api.genome.ucsc.edu/getData/track"
        f"?genome=hg38;track={UCSC_COMMON_SNP_TRACK};chrom={chrom_ucsc};start={start0};end={end1}"
    )
    data = http_get_json(url, headers={"Content-Type": "application/json"})
    records = _extract_track_records(data, UCSC_COMMON_SNP_TRACK)
    if not isinstance(records, list):
        return []
    return records


def ucsc_common_snps_to_local_gff3(snps, chrom_no_chr: str, start1: int, end1: int, local_seqid: str):
    lines = ["##gff-version 3"]
    local_records = []
    locus_start0 = start1 - 1

    for snp in snps:
        chrom_start0 = snp.get("chromStart")
        chrom_end1 = snp.get("chromEnd")
        if chrom_start0 is None or chrom_end1 is None:
            continue

        chrom_start0 = int(chrom_start0)
        chrom_end1 = int(chrom_end1)
        local_start1 = chrom_start0 - locus_start0 + 1
        local_end1 = chrom_end1 - locus_start0
        if local_start1 < 1 or local_end1 < local_start1:
            continue

        snp_id = str(snp.get("name") or f"snp_{chrom_no_chr}_{chrom_start0+1}")
        ref = sanitize_attr_value(snp.get("ref") or "")
        alts = sanitize_attr_value(snp.get("alts") or "")
        snp_class = sanitize_attr_value(snp.get("class") or "")
        maf = sanitize_attr_value(snp.get("minorAlleleFreq") or "")
        minor = sanitize_attr_value(snp.get("minorAllele") or "")
        major = sanitize_attr_value(snp.get("majorAllele") or "")

        attrs = [
            f"ID={sanitize_attr_value(snp_id)}",
            f"Name={sanitize_attr_value(snp_id)}",
            f"label={sanitize_attr_value(snp_id)}",
            "source=UCSC_common_dbSNP",
            "genome=hg38",
            f"orig_coord={chrom_no_chr}:{chrom_start0+1}-{chrom_end1}",
        ]
        if ref:
            attrs.append(f"ref={ref}")
        if alts:
            attrs.append(f"alts={alts}")
        if snp_class:
            attrs.append(f"class={snp_class}")
        if maf:
            attrs.append(f"minorAlleleFreq={maf}")
        if major:
            attrs.append(f"majorAllele={major}")
        if minor:
            attrs.append(f"minorAllele={minor}")

        lines.append(
            "\t".join(
                [
                    local_seqid,
                    "UCSC",
                    "variation",
                    str(local_start1),
                    str(local_end1),
                    ".",
                    ".",
                    ".",
                    ";".join(attrs),
                ]
            )
        )

        local_records.append(
            {
                "name": snp_id,
                "local_start1": local_start1,
                "local_end1": local_end1,
                "chrom": chrom_no_chr,
                "genomic_start1": chrom_start0 + 1,
                "genomic_end1": chrom_end1,
                "ref": str(snp.get("ref") or ""),
                "alts": str(snp.get("alts") or ""),
                "minorAlleleFreq": str(snp.get("minorAlleleFreq") or ""),
                "freqSourceCount": str(snp.get("freqSourceCount") or ""),
            }
        )

    return "\n".join(lines) + "\n", local_records


def summarize_minor_allele_freq(maf_text: str, freq_source_count_text: str) -> str:
    values = []
    for token in str(maf_text or "").split(","):
        token = token.strip()
        if not token:
            continue
        try:
            value = float(token)
        except ValueError:
            continue
        if value == float("inf") or value == float("-inf"):
            continue
        if value != value:
            continue
        values.append(value)

    reported = len(values)
    try:
        total = int(str(freq_source_count_text or "").strip())
    except ValueError:
        total = reported

    if reported == 0:
        return "MAF unavailable"

    max_value = max(values)
    return f"MAF reported by {reported}/{total} sources; max={max_value:.3f}"


def guide_pam_interval(feature: dict) -> tuple[int, int] | None:
    start1 = int(feature.get("start1", 0))
    end1 = int(feature.get("end1", 0))
    strand = str(feature.get("strand") or ".")
    if start1 < 1 or end1 < start1:
        return None
    if strand == "+":
        return end1 + 1, end1 + 3
    if strand == "-":
        pam_start1 = max(1, start1 - 3)
        pam_end1 = start1 - 1
        if pam_end1 < pam_start1:
            return None
        return pam_start1, pam_end1
    return None


def extract_interval_sequence(seq: str, start1: int, end1: int) -> str:
    if start1 < 1 or end1 < start1 or end1 > len(seq):
        return ""
    return seq[start1 - 1:end1].upper()


def format_warning_lines(items: list[dict], message_builder, indent: str = "  - ") -> list[str]:
    grouped = {}
    order = []
    for item in items:
        key = (item.get("pool", "?"), item.get("guide", "guide"))
        if key not in grouped:
            grouped[key] = []
            order.append(key)
        grouped[key].append(item)

    lines = []
    for pool, guide in order:
        lines.append(f"{guide} ({pool})")
        for item in grouped[(pool, guide)]:
            lines.append(indent + message_builder(item))
    return lines


def find_guides_with_non_ngg_pam(guide_gff_paths, ref_seq: str):
    warnings = []
    for pool_name, gff_path in guide_gff_paths:
        for feat in parse_gff_features(gff_path):
            guide_name = feat["attrs"].get("Name") or feat["attrs"].get("label") or feat["attrs"].get("ID") or "guide"
            pam_interval = guide_pam_interval(feat)
            if pam_interval is None:
                continue
            pam_start1, pam_end1 = pam_interval
            pam_seq = extract_interval_sequence(ref_seq, pam_start1, pam_end1)
            if len(pam_seq) != 3:
                continue

            strand = feat.get("strand", ".")
            is_bad_pam = False

            if strand == "+":
                is_bad_pam = pam_seq[1:] != "GG"
            elif strand == "-":
                is_bad_pam = pam_seq[:2] != "CC"
            else:
                continue

            if is_bad_pam:
                warnings.append(
                    {
                        "pool": pool_name,
                        "guide": guide_name,
                        "guide_start1": feat["start1"],
                        "guide_end1": feat["end1"],
                        "guide_strand": strand,
                        "pam_start1": pam_start1,
                        "pam_end1": pam_end1,
                        "pam_seq": pam_seq,
                    }
                )
    return warnings


class PureTargetGenBankTranslator(BiopythonTranslator):
    def compute_feature_label(self, feature):
        qualifiers = getattr(feature, "qualifiers", {}) or {}
        source = str((qualifiers.get("source") or [""])[0]).lower()
        feature_type = str(getattr(feature, "type", "")).lower()

        if source == "ensembl" and feature_type == "exon":
            for key in ("label", "Name"):
                value = qualifiers.get(key)
                if value:
                    return str(value[0])[:40]
            return None

        for key in ("label", "Name", "name", "gene"):
            value = qualifiers.get(key)
            if value:
                return str(value[0])[:40]
        feature_id = qualifiers.get("ID")
        if feature_id:
            return str(feature_id[0])[:40]
        return None

    def compute_feature_color(self, feature):
        qualifiers = getattr(feature, "qualifiers", {}) or {}
        source = str((qualifiers.get("source") or [""])[0]).lower()
        feature_type = str(getattr(feature, "type", "")).lower()

        if source == "ensembl":
            if feature_type == "exon":
                return "#7d8da6"
            return "#c7cedb"
        if source == "bowtie":
            name = str((qualifiers.get("Name") or qualifiers.get("label") or [""])[0]).lower()
            if "even" in name:
                return "#4f81bd"
            if "odd" in name:
                return "#c0504d"
            return "#6d9eeb"
        if feature_type == "misc_feature":
            name = str((qualifiers.get("Name") or qualifiers.get("label") or [""])[0]).upper()
            if name == "TILE":
                return "#9bbb59"
        return "#b7b7b7"

    def compute_feature_box_color(self, feature):
        return self.compute_feature_color(feature)

    def compute_feature_linewidth(self, feature):
        qualifiers = getattr(feature, "qualifiers", {}) or {}
        source = str((qualifiers.get("source") or [""])[0]).lower()
        if source == "ensembl":
            return 0.8
        if source == "bowtie":
            return 1.0
        return 0.6

    def compute_filtered_features(self, features):
        kept = []
        for feature in features:
            qualifiers = getattr(feature, "qualifiers", {}) or {}
            source = str((qualifiers.get("source") or [""])[0]).lower()
            feature_type = str(getattr(feature, "type", "")).lower()

            if source in {"ensembl", "bowtie"}:
                kept.append(feature)
                continue

            if feature_type == "misc_feature":
                name = str((qualifiers.get("Name") or qualifiers.get("label") or [""])[0]).upper()
                if name == "TILE":
                    kept.append(feature)

        return kept


def render_genbank_preview_png(gbk_path: Path, out_path: Path):
    with open(gbk_path, "r", encoding="utf-8") as handle:
        record = SeqIO.read(handle, "genbank")

    translator = PureTargetGenBankTranslator()
    graphic_record = translator.translate_record(record)

    seq_len = len(record.seq)
    figure_width = max(14, min(40, seq_len / 2500))
    figure_height = 6

    fig, ax = plt.subplots(1, 1, figsize=(figure_width, figure_height))
    graphic_record.plot(ax=ax, strand_in_label_threshold=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def gff_preview_label(feature: dict) -> str | None:
    attrs = feature.get("attrs", {}) or {}
    source = str(feature.get("source") or "").lower()
    feature_type = str(feature.get("type") or "").lower()

    if source == "ensembl" and feature_type == "exon":
        for key in ("label", "Name"):
            value = attrs.get(key)
            if value:
                return str(value)[:40]
        return None

    for key in ("label", "Name", "name", "gene"):
        value = attrs.get(key)
        if value:
            return str(value)[:40]

    feature_id = attrs.get("ID")
    if feature_id:
        return str(feature_id)[:40]
    return None


def gff_preview_color(feature: dict) -> str:
    attrs = feature.get("attrs", {}) or {}
    source = str(feature.get("source") or "").lower()
    feature_type = str(feature.get("type") or "").lower()

    if source == "ensembl":
        if feature_type == "exon":
            return "#7d8da6"
        return "#c7cedb"
    if source == "bowtie":
        name = str(attrs.get("Name") or attrs.get("label") or "").lower()
        if "even" in name:
            return "#4f81bd"
        if "odd" in name:
            return "#c0504d"
        return "#6d9eeb"
    if feature_type == "misc_feature":
        name = str(attrs.get("Name") or attrs.get("label") or "").upper()
        if name == "TILE":
            return "#9bbb59"
    return "#b7b7b7"


def keep_gff_preview_feature(feature: dict) -> bool:
    attrs = feature.get("attrs", {}) or {}
    source = str(feature.get("source") or "").lower()
    feature_type = str(feature.get("type") or "").lower()

    if source in {"ensembl", "bowtie"}:
        return True

    if feature_type == "misc_feature":
        name = str(attrs.get("Name") or attrs.get("label") or "").upper()
        return name == "TILE"

    return False


def gff_strand_to_int(strand: str) -> int:
    if strand == "+":
        return 1
    if strand == "-":
        return -1
    return 0


def render_gff_preview_png(gff_path: Path, fasta_path: Path, out_path: Path):
    with open(fasta_path, "r", encoding="utf-8") as handle:
        record = SeqIO.read(handle, "fasta")

    seq_len = len(record.seq)
    graphic_features = []
    for feature in parse_gff_features(gff_path):
        if not keep_gff_preview_feature(feature):
            continue

        start0 = max(0, int(feature["start1"]) - 1)
        end0 = min(seq_len, int(feature["end1"]))
        if end0 <= start0:
            continue

        graphic_features.append(
            GraphicFeature(
                start=start0,
                end=end0,
                strand=gff_strand_to_int(feature.get("strand", ".")),
                color=gff_preview_color(feature),
                label=gff_preview_label(feature),
            )
        )

    graphic_record = GraphicRecord(sequence_length=seq_len, features=graphic_features)
    figure_width = max(14, min(40, seq_len / 2500))
    figure_height = 6

    fig, ax = plt.subplots(1, 1, figsize=(figure_width, figure_height))
    graphic_record.plot(ax=ax, strand_in_label_threshold=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def maybe_render_genbank_preview(locus_dir: Path, gbk_path: Path):
    preview_path = locus_dir / "custom.png"
    try:
        render_genbank_preview_png(gbk_path, preview_path)
    except Exception as e:
        return None, f"{type(e).__name__}: {e}"
    return preview_path, None


def maybe_render_gff_preview(locus_dir: Path, gff_path: Path, fasta_path: Path):
    preview_path = locus_dir / "custom.png"
    try:
        render_gff_preview_png(gff_path, fasta_path, preview_path)
    except Exception as e:
        return None, f"{type(e).__name__}: {e}"
    return preview_path, None


def normalize_non_stranded_tiles_in_genbank(gbk_path: Path):
    with open(gbk_path, "r", encoding="utf-8") as handle:
        record = SeqIO.read(handle, "genbank")

    changed = False
    for feature in record.features:
        qualifiers = getattr(feature, "qualifiers", {}) or {}
        feature_type = str(getattr(feature, "type", "")).lower()
        name = str((qualifiers.get("Name") or qualifiers.get("label") or [""])[0]).upper()
        if feature_type == "misc_feature" and name == "TILE":
            if getattr(feature.location, "strand", None) is not None:
                feature.location.strand = None
                changed = True

    if changed:
        with open(gbk_path, "w", encoding="utf-8") as handle:
            SeqIO.write(record, handle, "genbank")


def show_preview_section(preview_path: Path | None, *, title: str | None = None, warning: str | None = None):
    if title:
        st.subheader(title)
    if warning:
        st.warning(f"Annotation outputs were generated, but preview rendering failed: {warning}")
        return
    if preview_path is None or not preview_path.exists():
        return
    st.image(str(preview_path), caption=preview_path.name, use_container_width=True)


def render_multi_locus_preview_gallery(preview_rows: list[dict]):
    available = [row for row in preview_rows if row.get("preview_path") and Path(row["preview_path"]).exists()]
    if not available:
        failed = [row for row in preview_rows if row.get("preview_error")]
        if failed:
            st.warning("Annotated files were created, but preview PNGs could not be rendered for any locus.")
        return

    st.subheader("Locus preview")
    for row in available:
        st.markdown(f"**{row['locus']}**")
        st.image(str(row["preview_path"]), caption=f"{row['locus']} preview", use_container_width=True)
        st.download_button(
            f"Download preview PNG ({row['locus']})",
            data=Path(row["preview_path"]).read_bytes(),
            file_name=f"{row['locus_slug']}_preview.png",
            mime="image/png",
            key=f"download_preview_{row['locus_slug']}",
        )

    failed = [row for row in preview_rows if row.get("preview_error")]
    for row in failed[:5]:
        st.warning(f"Preview rendering failed for {row['locus']}: {row['preview_error']}")
    if len(failed) > 5:
        st.warning(f"... and {len(failed) - 5} more preview rendering failures")


def parse_gff_features(path: Path):
    features = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            attrs = {}
            for item in cols[8].split(";"):
                if not item:
                    continue
                if "=" in item:
                    k, v = item.split("=", 1)
                    attrs[k] = v
            features.append(
                {
                    "seqid": cols[0],
                    "source": cols[1],
                    "type": cols[2],
                    "start1": int(cols[3]),
                    "end1": int(cols[4]),
                    "score": cols[5],
                    "strand": cols[6],
                    "phase": cols[7],
                    "attrs": attrs,
                }
            )
    return features


def find_guide_snp_overlaps(guide_gff_paths, snp_records):
    guide_overlaps = []
    pam_overlaps = []
    for pool_name, gff_path in guide_gff_paths:
        for feat in parse_gff_features(gff_path):
            guide_name = feat["attrs"].get("Name") or feat["attrs"].get("label") or feat["attrs"].get("ID") or "guide"
            pam_interval = guide_pam_interval(feat)
            for snp in snp_records:
                maf_summary = summarize_minor_allele_freq(
                    snp.get("minorAlleleFreq", ""),
                    snp.get("freqSourceCount", ""),
                )
                common_info = {
                    "pool": pool_name,
                    "guide": guide_name,
                    "guide_start1": feat["start1"],
                    "guide_end1": feat["end1"],
                    "guide_strand": feat["strand"],
                    "snp": snp["name"],
                    "snp_local_start1": snp["local_start1"],
                    "snp_local_end1": snp["local_end1"],
                    "snp_genomic_coord": f"{snp['chrom']}:{snp['genomic_start1']}-{snp['genomic_end1']}",
                    "maf_summary": maf_summary,
                }
                if snp["local_start1"] <= feat["end1"] and snp["local_end1"] >= feat["start1"]:
                    guide_overlaps.append(common_info.copy())
                if pam_interval is not None:
                    pam_start1, pam_end1 = pam_interval
                    if snp["local_start1"] <= pam_end1 and snp["local_end1"] >= pam_start1:
                        pam_info = common_info.copy()
                        pam_info["pam_start1"] = pam_start1
                        pam_info["pam_end1"] = pam_end1
                        pam_overlaps.append(pam_info)
    return guide_overlaps, pam_overlaps


def sanitize_attr_value(s: str) -> str:
    return str(s).replace(";", ",").replace("=", "_").replace("\t", " ").strip()


def strand_to_gff(strand) -> str:
    if strand == 1 or str(strand) == "1":
        return "+"
    if strand == -1 or str(strand) == "-1":
        return "-"
    return "."


def clipped_local_interval(
    feature_start1,
    feature_end1,
    locus_start1: int,
    locus_end1: int,
) -> tuple[int, int, int, int] | None:
    if feature_start1 is None or feature_end1 is None:
        return None

    clipped_start = max(int(feature_start1), locus_start1)
    clipped_end = min(int(feature_end1), locus_end1)
    if clipped_start > clipped_end:
        return None

    local_start = clipped_start - locus_start1 + 1
    local_end = clipped_end - locus_start1 + 1
    return clipped_start, clipped_end, local_start, local_end


def gene_display_label(gene: dict, fallback_id: str) -> str:
    label = (
        gene.get("external_name")
        or gene.get("gene_name")
        or gene.get("display_name")
        or gene.get("description")
        or fallback_id
    )
    return sanitize_attr_value(label)


def gene_span_to_local_gff3_line(
    gene: dict,
    chrom_no_chr: str,
    start1: int,
    end1: int,
    local_seqid: str,
) -> str | None:
    interval = clipped_local_interval(gene.get("start"), gene.get("end"), start1, end1)
    if interval is None:
        return None

    clipped_start, clipped_end, local_start, local_end = interval
    gene_id = str(gene.get("id") or f"gene_{chrom_no_chr}_{gene.get('start')}_{gene.get('end')}")
    label = gene_display_label(gene, gene_id)

    attrs = (
        f"ID={sanitize_attr_value(gene_id)};"
        f"Name={label};"
        f"label={label};"
        f"source=Ensembl;"
        f"genome=hg38;"
        f"orig_coord={chrom_no_chr}:{clipped_start}-{clipped_end}"
    )

    return "\t".join(
        [
            local_seqid,
            "Ensembl",
            "misc_feature",
            str(local_start),
            str(local_end),
            ".",
            strand_to_gff(gene.get("strand")),
            ".",
            attrs,
        ]
    )


def transcript_exons_to_local_gff3_lines(
    gene: dict,
    transcript: dict,
    chrom_no_chr: str,
    start1: int,
    end1: int,
    local_seqid: str,
) -> list[str]:
    gene_id = str(gene.get("id") or transcript.get("Parent") or "gene")
    transcript_id = str(transcript.get("id") or f"{gene_id}.main_transcript")
    label = gene_display_label(gene, gene_id)
    tx_label = sanitize_attr_value(transcript.get("display_name") or transcript_id)
    strand_char = strand_to_gff(transcript.get("strand", gene.get("strand")))
    exons = transcript.get("Exon") or []
    exons = sorted(
        exons,
        key=lambda exon: int(exon.get("start") or 0),
        reverse=(strand_char == "-"),
    )

    visible_exons = []
    for exon_number, exon in enumerate(exons, start=1):
        interval = clipped_local_interval(exon.get("start"), exon.get("end"), start1, end1)
        if interval is None:
            continue
        clipped_start, clipped_end, local_start, local_end = interval
        visible_exons.append(
            {
                "exon": exon,
                "exon_number": exon_number,
                "clipped_start": clipped_start,
                "clipped_end": clipped_end,
                "local_start": local_start,
                "local_end": local_end,
            }
        )

    lines = []
    mane_entries = transcript.get("MANE") or []
    mane_type = ""
    for entry in mane_entries:
        if isinstance(entry, dict) and entry.get("type"):
            mane_type = sanitize_attr_value(entry.get("type"))
            break

    for visible_index, item in enumerate(visible_exons, start=1):
        exon = item["exon"]
        exon_id = str(
            exon.get("id")
            or f"{transcript_id}.exon{item['exon_number']}"
        )
        attrs = [
            f"ID={sanitize_attr_value(transcript_id)}.exon{item['exon_number']}",
            f"gene={label}",
            f"transcript_id={sanitize_attr_value(transcript_id)}",
            f"transcript_name={tx_label}",
            f"source=Ensembl",
            f"genome=hg38",
            f"orig_coord={chrom_no_chr}:{item['clipped_start']}-{item['clipped_end']}",
            f"exon_number={item['exon_number']}",
            f"ensembl_exon_id={sanitize_attr_value(exon_id)}",
        ]
        if visible_index == 1:
            attrs.append(f"Name={label}")
            attrs.append(f"label={label}")
        if mane_type:
            attrs.append(f"mane={mane_type}")

        lines.append(
            "\t".join(
                [
                    local_seqid,
                    "Ensembl",
                    "exon",
                    str(item["local_start"]),
                    str(item["local_end"]),
                    ".",
                    strand_char,
                    ".",
                    ";".join(attrs),
                ]
            )
        )

    return lines


def genes_to_local_gff3(genes, chrom_no_chr: str, start1: int, end1: int, local_seqid: str) -> str:
    """
    Rebase hg38 gene annotations onto the extracted sequence.
    Prefer the MANE Select transcript, then canonical transcript, and write
    visible exons so the preview shows exon/intron structure.
    """
    lines = ["##gff-version 3"]

    for gene in genes:
        gene_id = gene.get("id")
        exon_lines = []
        if gene_id:
            gene_model = fetch_ensembl_gene_model(str(gene_id))
            transcript = select_main_transcript(gene_model)
            if transcript is not None:
                exon_lines = transcript_exons_to_local_gff3_lines(
                    gene,
                    transcript,
                    chrom_no_chr,
                    start1,
                    end1,
                    local_seqid,
                )

        if exon_lines:
            lines.extend(exon_lines)
            continue

        span_line = gene_span_to_local_gff3_line(
            gene,
            chrom_no_chr,
            start1,
            end1,
            local_seqid,
        )
        if span_line is not None:
            lines.append(span_line)

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


def count_gff_features(path: Path) -> int:
    count = 0
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                count += 1
    return count


def write_zip_from_directory(root_dir: Path) -> bytes:
    bio = io.BytesIO()
    with zipfile.ZipFile(bio, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(root_dir.rglob("*")):
            if path.is_file():
                arcname = path.relative_to(root_dir)
                zf.write(path, arcname=str(arcname))
    bio.seek(0)
    return bio.getvalue()


def hg38_chrom_name(chrom_no_chr: str) -> str:
    return f"chr{chrom_no_chr}"


def format_gff3_attrs(attrs: dict) -> str:
    items = []
    for key, value in attrs.items():
        if key:
            items.append(f"{sanitize_attr_value(key)}={sanitize_attr_value(value)}")
    return ";".join(items) if items else "."


def write_genome_rebased_gff3(
    input_gff: Path,
    output_gff: Path,
    chrom_no_chr: str,
    locus_start1: int,
    *,
    source_override: str | None = None,
    feature_name_override: str | None = None,
    id_prefix_override: str | None = None,
):
    chrom = hg38_chrom_name(chrom_no_chr)
    lines = ["##gff-version 3"]

    for feature in parse_gff_features(input_gff):
        local_start1 = int(feature["start1"])
        local_end1 = int(feature["end1"])
        genome_start1 = locus_start1 + local_start1 - 1
        genome_end1 = locus_start1 + local_end1 - 1
        attrs = dict(feature.get("attrs", {}) or {})

        if source_override:
            source = source_override
        else:
            source = feature["source"]

        if feature_name_override:
            attrs["Name"] = feature_name_override
            attrs["label"] = feature_name_override

        if id_prefix_override and attrs.get("ID"):
            attrs["ID"] = re.sub(r"^[^.]+", id_prefix_override, attrs["ID"], count=1)
        elif source == "bowtie":
            guide_name = attrs.get("Name") or attrs.get("label") or attrs.get("ID") or "guide"
            attrs["ID"] = f"{guide_name}.{chrom}.{genome_start1}.{feature['strand']}"

        attrs["genome"] = "hg38"
        attrs["local_seqid"] = feature["seqid"]
        attrs["local_coord"] = f"{local_start1}-{local_end1}"

        lines.append(
            "\t".join(
                [
                    chrom,
                    source,
                    feature["type"],
                    str(genome_start1),
                    str(genome_end1),
                    feature["score"],
                    feature["strand"],
                    feature["phase"],
                    format_gff3_attrs(attrs),
                ]
            )
        )

    output_gff.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_hg38_igv_export_files(
    export_locus_dir: Path,
    outputs: dict,
    chrom_no_chr: str,
    locus_start1: int,
) -> dict:
    export_locus_dir.mkdir(parents=True, exist_ok=True)

    exported = {
        "guides_even": export_locus_dir / "guides_even.gff3",
        "tiles_even": export_locus_dir / "tiles_even.gff3",
        "guides_odd": None,
        "tiles_odd": None,
    }

    write_genome_rebased_gff3(
        outputs["even_gff"],
        exported["guides_even"],
        chrom_no_chr,
        locus_start1,
    )
    write_genome_rebased_gff3(
        outputs["even_tiles"],
        exported["tiles_even"],
        chrom_no_chr,
        locus_start1,
        source_override="tiles",
        feature_name_override="tiles",
        id_prefix_override="TILES",
    )

    if outputs.get("odd_gff") is not None and Path(outputs["odd_gff"]).exists():
        exported["guides_odd"] = export_locus_dir / "guides_odd.gff3"
        write_genome_rebased_gff3(
            outputs["odd_gff"],
            exported["guides_odd"],
            chrom_no_chr,
            locus_start1,
        )

    if outputs.get("odd_tiles") is not None and Path(outputs["odd_tiles"]).exists():
        exported["tiles_odd"] = export_locus_dir / "tiles_odd.gff3"
        write_genome_rebased_gff3(
            outputs["odd_tiles"],
            exported["tiles_odd"],
            chrom_no_chr,
            locus_start1,
            source_override="tiles",
            feature_name_override="tiles",
            id_prefix_override="TILES",
        )

    return exported


def map_guides_to_gff(guides_fa_path: Path, mapped_out: Path, gff_out: Path, label: str, wd: Path, log):
    log(f"bowtie map {label}")
    rc, out, err = run_cmd(
        ["bowtie", "-f", "-v", "0", "-a", "--best", "--strata", "refidx", str(guides_fa_path)],
        wd,
    )
    if rc != 0:
        fail(f"bowtie map ({label})", out, err)
    mapped_out.write_text(out, encoding="utf-8")

    log(f"mapped -> GFF3 {label}")
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


def run_single_locus_pipeline(
    locus_dir: Path,
    ref_path: Path,
    regions_gff: Path,
    even_path: Path,
    odd_path: Path | None,
    max_between_len: int,
    log,
    make_genbank_output: bool = True,
):
    even_mapped = locus_dir / "guides_even.mapped"
    odd_mapped = locus_dir / "guides_odd.mapped"
    even_gff = locus_dir / "guides_even.gff3"
    odd_gff = locus_dir / "guides_odd.gff3"
    even_tiles = locus_dir / "between_regions_even.gff3"
    odd_tiles = locus_dir / "between_regions_odd.gff3"
    combined_gff = locus_dir / "combined.gff3"
    gbk_path = locus_dir / "custom.gbk"
    preview_png_path = locus_dir / "custom.png"

    log("bowtie-build")
    rc, out, err = run_cmd(["bowtie-build", str(ref_path), "refidx"], locus_dir)
    if rc != 0:
        fail("bowtie-build", out, err)

    map_guides_to_gff(even_path, even_mapped, even_gff, "even", locus_dir, log)
    if odd_path is not None:
        map_guides_to_gff(odd_path, odd_mapped, odd_gff, "odd", locus_dir, log)

    log("predict even tiles")
    rc, out, err = run_cmd(
        [
            "python",
            str(BETWEEN_SCRIPT),
            "--gff",
            str(even_gff),
            "--fasta",
            str(ref_path),
            "--out",
            str(even_tiles),
            "--maxlen",
            str(int(max_between_len)),
        ],
        locus_dir,
    )
    if rc != 0:
        fail("make_between_regions.py (even)", out, err)
    log(out.strip() or "between even done")

    if odd_path is not None:
        log("predict odd tiles")
        rc, out, err = run_cmd(
            [
                "python",
                str(BETWEEN_SCRIPT),
                "--gff",
                str(odd_gff),
                "--fasta",
                str(ref_path),
                "--out",
                str(odd_tiles),
                "--maxlen",
                str(int(max_between_len)),
            ],
            locus_dir,
        )
        if rc != 0:
            fail("make_between_regions.py (odd)", out, err)
        log(out.strip() or "between odd done")

    log("combine GFF3s")
    combined = ["##gff-version 3"]
    combined.append(strip_header(regions_gff.read_text(encoding="utf-8")))
    combined.append(strip_header(even_gff.read_text(encoding="utf-8")))
    if odd_path is not None:
        combined.append(strip_header(odd_gff.read_text(encoding="utf-8")))
    combined.append(strip_header(even_tiles.read_text(encoding="utf-8")))
    if odd_path is not None:
        combined.append(strip_header(odd_tiles.read_text(encoding="utf-8")))
    combined_gff.write_text("\n".join([c for c in combined if c != ""]) + "\n", encoding="utf-8")

    if make_genbank_output:
        log("GFF3+FASTA -> GenBank")
        try:
            gff_fasta_to_genbank(str(combined_gff), str(ref_path), str(gbk_path))
        except Exception as e:
            st.error("Failed at: GenBank conversion")
            st.exception(e)
            st.stop()

        normalize_non_stranded_tiles_in_genbank(gbk_path)
        preview_path, preview_error = maybe_render_genbank_preview(locus_dir, gbk_path)
    else:
        log("render preview PNG")
        preview_path, preview_error = maybe_render_gff_preview(locus_dir, combined_gff, ref_path)

    return {
        "even_mapped": even_mapped,
        "odd_mapped": odd_mapped if odd_path is not None else None,
        "even_gff": even_gff,
        "odd_gff": odd_gff if odd_path is not None else None,
        "even_tiles": even_tiles,
        "odd_tiles": odd_tiles if odd_path is not None else None,
        "combined_gff": combined_gff,
        "gbk_path": gbk_path if make_genbank_output else None,
        "preview_png_path": preview_path or preview_png_path,
        "preview_error": preview_error,
    }


# -----------------------------
# UI
# -----------------------------

input_mode = st.radio(
    "Reference source",
    ["Fetch hg38 locus by coordinates", "Use uploaded custom reference"],
)

if input_mode == "Use uploaded custom reference":
    ref_fa = st.file_uploader("Reference FASTA (custom_reference.fa)", type=["fa", "fasta", "fna"])
    regions_bed = st.file_uploader("Regions BED4 (custom_regions.bed)", type=["bed", "txt"])
    locus_coords = None
else:
    locus_coords = st.text_area(
        "Locus coordinates or gene names (one per line)",
        placeholder="EGFR\nchr7:55019017-55211628\nMYC",
        height=150,
    )
    ref_fa = None
    regions_bed = None

guides_even = st.file_uploader("Guides FASTA (for single-tile design) or even guides FASTA (for multi-tile design)", type=["fa", "fas", "fasta", "fna"])
guides_odd = st.file_uploader("Odd guides FASTA (for multi-tile design; leave empty for single-tile design)", type=["fa", "fas", "fasta", "fna"])

with st.expander("Advanced", expanded=False):
    max_between_len = st.number_input(
        "Max fragment length (bp)", min_value=1, value=20000, step=1000
    )
    keep_intermediates = st.checkbox("Show & allow download of intermediate files", value=False)

if input_mode == "Use uploaded custom reference":
    ready = bool(ref_fa and regions_bed and guides_even)
else:
    ready = bool(locus_coords and guides_even)

run_btn = st.button("Run", type="primary", disabled=not ready)

# -----------------------------
# Main pipeline
# -----------------------------

if run_btn:
    with tempfile.TemporaryDirectory() as td:
        wd = Path(td)

        # Shared guide files
        even_path = wd / "guides_even.fa"
        odd_path = wd / "guides_odd.fa" if guides_odd is not None else None
        even_path.write_bytes(guides_even.getvalue())
        if odd_path is not None:
            odd_path.write_bytes(guides_odd.getvalue())

        log_area = st.empty()
        logs = []

        def log(msg):
            logs.append(msg)
            log_area.code("\n".join(logs))

        st.info("Running pipeline...")
        log("Starting...")

        if input_mode == "Use uploaded custom reference":
            # Existing single-locus behavior
            ref_path = wd / "custom_reference.fa"
            bed_path = wd / "custom_regions.bed"
            regions_gff = wd / "custom_regions.gff3"

            ref_path.write_bytes(ref_fa.getvalue())
            bed_path.write_bytes(regions_bed.getvalue())

            log("BED -> GFF3")
            regions_gff.write_text(bed_to_gff3_text(bed_path), encoding="utf-8")

            outputs = run_single_locus_pipeline(
                locus_dir=wd,
                ref_path=ref_path,
                regions_gff=regions_gff,
                even_path=even_path,
                odd_path=odd_path,
                max_between_len=int(max_between_len),
                log=log,
            )

            st.success("Done!")

            st.download_button(
                "Download GenBank (custom.gbk)",
                data=outputs["gbk_path"].read_bytes(),
                file_name="custom.gbk",
                mime="application/octet-stream",
            )
            st.download_button(
                "Download combined GFF3 (combined.gff3)",
                data=outputs["combined_gff"].read_bytes(),
                file_name="combined.gff3",
                mime="text/plain",
            )

            show_preview_section(
                outputs.get("preview_png_path"),
                title="Locus preview",
                warning=outputs.get("preview_error"),
            )
            if outputs.get("preview_png_path") and Path(outputs["preview_png_path"]).exists():
                st.download_button(
                    "Download preview PNG (custom.png)",
                    data=Path(outputs["preview_png_path"]).read_bytes(),
                    file_name="custom.png",
                    mime="image/png",
                )

            if keep_intermediates:
                st.subheader("Intermediates")
                st.download_button(
                    "custom_regions.gff3",
                    regions_gff.read_bytes(),
                    file_name="custom_regions.gff3",
                )
                st.download_button(
                    "guides_even.gff3",
                    outputs["even_gff"].read_bytes(),
                    file_name="guides_even.gff3",
                )
                if outputs.get("odd_gff") is not None and Path(outputs["odd_gff"]).exists():
                    st.download_button(
                        "guides_odd.gff3",
                        outputs["odd_gff"].read_bytes(),
                        file_name="guides_odd.gff3",
                    )
                st.download_button(
                    "between_regions_even.gff3",
                    outputs["even_tiles"].read_bytes(),
                    file_name="between_regions_even.gff3",
                )
                if outputs.get("odd_tiles") is not None and Path(outputs["odd_tiles"]).exists():
                    st.download_button(
                        "between_regions_odd.gff3",
                        outputs["odd_tiles"].read_bytes(),
                        file_name="between_regions_odd.gff3",
                    )
                st.download_button(
                    "guides_even.mapped",
                    outputs["even_mapped"].read_bytes(),
                    file_name="guides_even.mapped",
                )
                if outputs.get("odd_mapped") is not None and Path(outputs["odd_mapped"]).exists():
                    st.download_button(
                        "guides_odd.mapped",
                        outputs["odd_mapped"].read_bytes(),
                        file_name="guides_odd.mapped",
                    )

        else:
            # Multi-locus hg38 mode
            loci = parse_multi_loci_or_stop(locus_coords)
            results_root = wd / "results"
            results_root.mkdir(parents=True, exist_ok=True)
            export_root = wd / "igv_gff3_outputs"
            export_root.mkdir(parents=True, exist_ok=True)

            summary_rows = []
            preview_rows = []

            for idx, (chrom_no_chr, start1, end1) in enumerate(loci, start=1):
                locus_slug = make_locus_slug(chrom_no_chr, start1, end1)
                locus_dir = results_root / locus_slug
                locus_dir.mkdir(parents=True, exist_ok=True)

                ref_path = locus_dir / "custom_reference.fa"
                regions_gff = locus_dir / "custom_regions.gff3"

                log(f"[{idx}/{len(loci)}] Processing {chrom_no_chr}:{start1}-{end1}")

                local_seqid = make_local_seqid(chrom_no_chr, start1, end1)

                log(f"[{locus_slug}] fetch hg38 sequence")
                seq = fetch_hg38_sequence_ucsc(chrom_no_chr, start1, end1)
                ref_text = f">{local_seqid}\n{wrap_fasta(seq)}\n"
                ref_path.write_text(ref_text, encoding="utf-8")

                log(f"[{locus_slug}] fetch Ensembl gene annotations")
                genes = fetch_ensembl_gene_overlaps(chrom_no_chr, start1, end1)

                log(f"[{locus_slug}] fetch UCSC common SNPs")
                raw_snps = fetch_ucsc_common_snps(chrom_no_chr, start1, end1)

                log(f"[{locus_slug}] convert Ensembl main transcripts -> local GFF3")
                regions_gff.write_text(
                    genes_to_local_gff3(genes, chrom_no_chr, start1, end1, local_seqid),
                    encoding="utf-8",
                )

                _, snp_records = ucsc_common_snps_to_local_gff3(
                    raw_snps, chrom_no_chr, start1, end1, local_seqid
                )

                outputs = run_single_locus_pipeline(
                    locus_dir=locus_dir,
                    ref_path=ref_path,
                    regions_gff=regions_gff,
                    even_path=even_path,
                    odd_path=odd_path,
                    max_between_len=int(max_between_len),
                    log=lambda msg, prefix=locus_slug: log(f"[{prefix}] {msg}"),
                    make_genbank_output=False,
                )

                export_paths = write_hg38_igv_export_files(
                    export_locus_dir=export_root / locus_slug,
                    outputs=outputs,
                    chrom_no_chr=chrom_no_chr,
                    locus_start1=start1,
                )

                guide_gff_paths = [("even", outputs["even_gff"])]
                if outputs.get("odd_gff") is not None and Path(outputs["odd_gff"]).exists():
                    guide_gff_paths.append(("odd", outputs["odd_gff"]))

                guide_snp_overlaps, pam_snp_overlaps = find_guide_snp_overlaps(
                    guide_gff_paths,
                    snp_records,
                )
                non_ngg_pam_warnings = find_guides_with_non_ngg_pam(
                    guide_gff_paths,
                    seq,
                )
                if guide_snp_overlaps or pam_snp_overlaps:
                    st.warning(f"Common SNP warning for {chrom_no_chr}:{start1}-{end1}")

                    shown_count = 0

                    for item in guide_snp_overlaps[:10]:
                        line = f"{item['guide']} ({item['pool']}) overlaps {item['snp']} at {item['snp_genomic_coord']}"
                        if item.get("maf_summary"):
                            line += f" [{item['maf_summary']}]"
                        st.warning(line)
                        shown_count += 1

                    for item in pam_snp_overlaps[:10]:
                        line = f"{item['guide']} ({item['pool']}) has a common SNP in PAM: {item['snp']} at {item['snp_genomic_coord']}"
                        if item.get("maf_summary"):
                            line += f" [{item['maf_summary']}]"
                        st.warning(line)
                        shown_count += 1

                    n_total = len(guide_snp_overlaps) + len(pam_snp_overlaps)
                    if n_total > shown_count:
                        st.warning(f"... and {n_total - shown_count} more")

                if non_ngg_pam_warnings:
                    st.warning(f"Non-NGG PAM warning for {chrom_no_chr}:{start1}-{end1}")

                    shown_count = 0
                    for item in non_ngg_pam_warnings[:10]:
                        st.warning(
                            f"{item['guide']} ({item['pool']}) has non-NGG PAM"
                        )
                        shown_count += 1

                    if len(non_ngg_pam_warnings) > shown_count:
                        st.warning(f"... and {len(non_ngg_pam_warnings) - shown_count} more")

                summary_rows.append(
                    {
                        "locus": f"{chrom_no_chr}:{start1}-{end1}",
                        "seq_length_bp": len(seq),
                        "genes": len(genes),
                        "common_snps": len(snp_records),
                        "even_guide_hits": count_gff_features(outputs["even_gff"]),
                        "odd_guide_hits": count_gff_features(outputs["odd_gff"]) if outputs.get("odd_gff") is not None and Path(outputs["odd_gff"]).exists() else 0,
                        "even_tiles": count_gff_features(outputs["even_tiles"]),
                        "odd_tiles": count_gff_features(outputs["odd_tiles"]) if outputs.get("odd_tiles") is not None and Path(outputs["odd_tiles"]).exists() else 0,
                        "guide_outputs": ", ".join(
                            name
                            for name, path in [
                                ("guides_even.gff3", export_paths.get("guides_even")),
                                ("guides_odd.gff3", export_paths.get("guides_odd")),
                            ]
                            if path is not None
                        ),
                        "tile_outputs": ", ".join(
                            name
                            for name, path in [
                                ("tiles_even.gff3", export_paths.get("tiles_even")),
                                ("tiles_odd.gff3", export_paths.get("tiles_odd")),
                            ]
                            if path is not None
                        ),
                    }
                )
                preview_rows.append(
                    {
                        "locus": f"{chrom_no_chr}:{start1}-{end1}",
                        "locus_slug": locus_slug,
                        "preview_path": outputs.get("preview_png_path"),
                        "preview_error": outputs.get("preview_error"),
                    }
                )

            zip_bytes = write_zip_from_directory(export_root)

            st.success(f"Done! Processed {len(summary_rows)} locus/loci.")
            st.subheader("Per-locus summary")
            st.dataframe(summary_rows, use_container_width=True)
            render_multi_locus_preview_gallery(preview_rows)

            st.download_button(
                "Download IGV GFF3 outputs (.zip)",
                data=zip_bytes,
                file_name="hg38_igv_gff3_outputs.zip",
                mime="application/zip",
            )

            if keep_intermediates:
                st.caption(
                    "The ZIP contains one folder per locus with genome-coordinate guides_even.gff3, optional guides_odd.gff3, tiles_even.gff3, and optional tiles_odd.gff3 for IGV."
                )
