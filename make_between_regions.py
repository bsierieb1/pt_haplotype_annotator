#!/usr/bin/env python3
import argparse
from pathlib import Path

def read_fasta(path: str) -> dict[str, str]:
    seqs = {}
    name = None
    parts = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(parts).upper()
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if name is not None:
            seqs[name] = "".join(parts).upper()
    return seqs

def parse_attrs(attr_str: str) -> dict[str, str]:
    d = {}
    for item in attr_str.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d

def read_gff3(path: str):
    feats = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = cols
            feats.append({
                "seqid": seqid,
                "source": source,
                "type": ftype,
                "start": int(start),
                "end": int(end),
                "strand": strand,
                "attrs": parse_attrs(attrs),
            })
    return feats

def has_N(refseq: str, start1: int, end1: int) -> bool:
    return "N" in refseq[start1-1:end1]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True, help="Input guides.gff3 (from bowtie mappings)")
    ap.add_argument("--fasta", required=True, help="Reference FASTA for N-check")
    ap.add_argument("--out", required=True, help="Output GFF3")
    ap.add_argument("--maxlen", type=int, default=20000, help="Max allowed gap length (default 20000)")
    ap.add_argument("--feature-type", default="misc_feature", help="GFF3 type for output features")
    ap.add_argument("--name", default="TILE", help="Name/label for output features")
    args = ap.parse_args()

    ref = read_fasta(args.fasta)
    feats = read_gff3(args.gff)

    by_seqid = {}
    for f in feats:
        by_seqid.setdefault(f["seqid"], []).append(f)

    out_lines = ["##gff-version 3"]
    kept = skipped_len = skipped_n = skipped_none = 0

    for seqid, flist in by_seqid.items():
        if seqid not in ref:
            raise SystemExit(f"ERROR: seqid '{seqid}' is in GFF3 but not in FASTA headers.")

        plus = [f for f in flist if f["strand"] == "+"]
        minus = [f for f in flist if f["strand"] == "-"]
        minus.sort(key=lambda x: x["start"])

        for p in plus:
            p_id = p["attrs"].get("ID", f"plus_{p['start']}_{p['end']}")
            found_any = False

            for m in minus:
                if m["start"] <= p["end"]:
                    continue  # not downstream of this + feature
                found_any = True

                gs = p["end"] + 1
                ge = m["start"] - 1
                if gs > ge:
                    continue

                glen = ge - gs + 1
                if glen > args.maxlen:
                    skipped_len += 1
                    continue

                if has_N(ref[seqid], gs, ge):
                    skipped_n += 1
                    continue

                m_id = m["attrs"].get("ID", f"minus_{m['start']}_{m['end']}")
                feature_id = f"BETWEEN.{p_id}.{m_id}.{gs}-{ge}"
                attrs = (
                    f"ID={feature_id};Name={args.name};label={args.name};"
                    f"plus={p_id};minus={m_id};len={glen}"
                )
                out_lines.append("\t".join([
                    seqid, "between", args.feature_type, str(gs), str(ge),
                    ".", ".", ".", attrs
                ]))
                kept += 1

            if not found_any:
                skipped_none += 1

    Path(args.out).write_text("\n".join(out_lines) + "\n", encoding="utf-8")

    # Small summary for logs
    print(f"Wrote {kept} between-features to {args.out}")
    print(f"Skipped (too long): {skipped_len}")
    print(f"Skipped (contains N): {skipped_n}")
    print(f"Plus guides with no downstream minus: {skipped_none}")

if __name__ == "__main__":
    main()
