# gff_to_genbank_patched.py
from __future__ import annotations

import re
from typing import Iterable

from Bio import SeqIO
from BCBio import GFF


def _fix_ncbi_id(records: Iterable):
    # Keep behavior similar to gff-to-genbank: ensure record.id is GenBank-safe
    for rec in records:
        rec.id = re.sub(r"[^A-Za-z0-9_.:-]", "_", rec.id)
        yield rec


def _extract_regions(gff_iter: Iterable):
    """
    gff-to-genbank modifies features in-place and ensures strand is present.
    The upstream bug is accessing feature.strand; strand lives on location.
    """
    for rec in gff_iter:
        # Ensure every feature has a strand via location.strand
        for feat in rec.features:
            try:
                _ = feat.location.strand
            except Exception:
                # If location missing or weird, force unknown
                feat.location.strand = None
        yield rec


def _check_gff(gff_iter: Iterable):
    # pass-through; upstream tool used this to validate
    for rec in gff_iter:
        yield rec


def gff_fasta_to_genbank(gff_path: str, fasta_path: str, out_path: str):
    from Bio import SeqIO
    from BCBio import GFF
    import re

    # Load FASTA records
    with open(fasta_path, "r", encoding="utf-8") as fasta_handle:
        fasta_records = list(SeqIO.parse(fasta_handle, "fasta"))
    base_dict = {r.id: r for r in fasta_records}

    # IMPORTANT: parse into a LIST while gff_handle is open (avoid closed-file generator)
    with open(gff_path, "r", encoding="utf-8") as gff_handle:
        annotated_records = list(GFF.parse(gff_handle, base_dict=base_dict))

    # Fix IDs to be GenBank-safe
    for rec in annotated_records:
        rec.id = re.sub(r"[^A-Za-z0-9_.:-]", "_", rec.id)

    # Write GenBank
    with open(out_path, "w", encoding="utf-8") as out_handle:
        SeqIO.write(annotated_records, out_handle, "genbank")
