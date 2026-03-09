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
    """
    Read FASTA + GFF3 and write GenBank.
    """
    # Parse GFF3 annotations onto FASTA sequences:
    with open(fasta_path, "r") as fasta_handle:
        fasta_records = list(SeqIO.parse(fasta_handle, "fasta"))

    with open(gff_path, "r") as gff_handle:
        annotated_records = GFF.parse(gff_handle, base_dict={r.id: r for r in fasta_records})

    records = _check_gff(_fix_ncbi_id(_extract_regions(annotated_records)))

    with open(out_path, "w") as out_handle:
        SeqIO.write(records, out_handle, "genbank")
