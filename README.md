# PureTarget tile annotator

**Testing app for annotating PureTarget loci from guide sets.**  
This is **not an official PacBio tool** and should be used with caution.

## What the app does

This Streamlit app takes a set of gRNAs designed for a PureTarget assay, maps them to a reference, predicts tiles, and produces annotated outputs including GFF3 and GenBank files.

## Input modes

The app currently supports two modes.

### 1) Provide the list of targeted genes or coordinates of targeted loci (human hg38 only)

Enter one item per line. Each line can be:

- a genomic interval, for example:
  - `chr7:55019017-55211628`
  - `7:140424943-140624564`
- a human gene symbol, for example:
  - `EGFR`
  - `CYP2D6`

Gene symbols are resolved against **human hg38 / GRCh38** using Ensembl.  
For gene-name inputs, the app uses the **full gene body plus 20 kb on each side**.

If multiple resolved loci overlap, they are **merged into a single locus** before downstream processing. This is useful for neighboring targets such as gene families or duplicated loci.

For each resolved locus, the app:

1. fetches the hg38 reference sequence from UCSC
2. fetches overlapping **gene-level** annotations from Ensembl
3. maps even and odd guides by exact match with Bowtie
4. predicts tiles between opposing guides
5. combines all annotations
6. generates a GenBank file

### 2) Upload custom target region FASTA

Upload:

- custom target region FASTA
- BED4 file with region annotations
- guides_even FASTA
- guides_odd FASTA

In this mode, the app:

1. converts BED4 annotations to GFF3
2. maps even and odd guides by exact match with Bowtie
3. predicts tiles between opposing guides
4. combines all annotations
5. generates a GenBank file

## Pipeline summary

Depending on mode, the app runs the following steps:

1. reference preparation
2. annotation preparation
3. Bowtie exact mapping for even + odd guides
4. tile annotation
5. combine all annotations
6. GFF3 + FASTA to GenBank conversion

## Tile prediction logic

Tiles are generated from guide mappings using the bundled `make_between_regions.py` helper.

Current behavior:

- only guides facing each other are paired
- maximum tile length is **20 kb by default** (there's a control to adjust it)
- regions containing `N` are excluded

## Outputs

### Custom reference mode

Primary outputs:

- `combined.gff3`
- `custom.gbk`

Optional intermediate outputs:

- `custom_regions.gff3`
- `guides_even.gff3`
- `guides_odd.gff3`
- `between_regions_even.gff3`
- `between_regions_odd.gff3`
- `guides_even.mapped`
- `guides_odd.mapped`

### hg38 multi-locus mode

The app creates one folder per processed locus and packages them into a ZIP archive.

Each locus folder contains:

- `custom_reference.fa`
- `custom_regions.gff3`
- `guides_even.gff3`
- `guides_odd.gff3`
- `between_regions_even.gff3`
- `between_regions_odd.gff3`
- `combined.gff3`
- `custom.gbk`
- Bowtie mapping intermediates

The app also displays a per-locus summary table in Streamlit.

## Input expectations

### Guide FASTA files

Two guide FASTA files are required:

- `guides_even.fa`
- `guides_odd.fa`

Guides are mapped by exact match using Bowtie.

### BED file

Custom-reference mode expects a BED4-style file:

- column 1: chromosome / sequence ID
- column 2: start
- column 3: end
- column 4: feature name

## Dependencies

### System packages

From `packages.txt`:

- `bowtie`
- `gawk`

### Python packages

From `requirements.txt`:

- `streamlit>=1.31`
- `biopython`
- `bcbio-gff`

## Project files

Main files in this repo:

- `streamlit_app.py` — main Streamlit app
- `make_between_regions.py` — tile prediction helper
- `gff_to_genbank_patched.py` — patched GFF3/FASTA to GenBank conversion helper
- `requirements.txt` — Python dependencies
- `packages.txt` — system dependencies for deployment
- `runtime.txt` — Python runtime version

## Notes and limitations

- This is a **testing / internal-style utility**, not a production-supported PacBio application.
- hg38 mode is intended for **human canonical loci on GRCh38/hg38**.
- Gene-name mode depends on remote Ensembl lookups.
- Sequence retrieval in hg38 mode depends on the UCSC API.
- Only **gene-level** Ensembl annotations are included in the locus GFF3 to keep outputs uncluttered.
- If multiple requested loci overlap after resolution, they are merged before processing.
