# PureTarget tile annotator

Testing app for annotating PureTarget loci from guide sets.  
Not an official PacBio tool.

## What it does

This Streamlit app maps PureTarget guides, predicts tiles, generates annotated GFF3 and GenBank outputs, and renders a preview PNG from the non-SNP GenBank file.

## Input modes

### 1) hg38 loci
Enter one target per line as either:
- a coordinate interval, for example `chr7:55019017-55211628`
- a human gene symbol, for example `EGFR`

Gene symbols are resolved with Ensembl and expanded to **full gene body ±20 kb**.  
Overlapping resolved loci are merged before processing.

For each locus, the app:
- fetches hg38 sequence from UCSC
- fetches gene annotations from Ensembl
- maps even and odd guides with Bowtie
- predicts tiles
- adds common SNP annotations from UCSC
- writes GFF3 and GenBank outputs
- renders a static PNG preview from the non-SNP GenBank output

It also reports:
- guides overlapping common SNPs (including common SNPs in the PAM)
- guides with a non-NGG PAM

### 2) Custom reference
Upload:
- reference FASTA
- BED4 annotations
- `guides_even.fa`
- `guides_odd.fa`

The app converts BED to GFF3, maps guides, predicts tiles, combines annotations, writes GFF3 + GenBank output, and renders a PNG preview from the non-SNP GenBank output.

## Tile logic

Tiles are predicted from mapped guides using `make_between_regions.py`.

Rules:
- only inward-facing guide pairs are used
- default max fragment length is 20 kb
- regions containing `N` are skipped

## Outputs

### Custom reference mode
Main outputs:
- `combined.gff3`
- `custom.gbk`
- `custom.png`

Optional intermediates:
- `custom_regions.gff3`
- `guides_even.gff3`
- `guides_odd.gff3`
- `between_regions_even.gff3`
- `between_regions_odd.gff3`
- `guides_even.mapped`
- `guides_odd.mapped`

### hg38 mode
Outputs are packaged as one folder per locus in a ZIP.

Typical files per locus:
- `custom_reference.fa`
- `custom_regions.gff3`
- `guides_even.gff3`
- `guides_odd.gff3`
- `between_regions_even.gff3`
- `between_regions_odd.gff3`
- `combined.gff3`
- `custom.gbk`
- `custom.png`
- `common_snps.gff3`
- `combined_with_common_snps.gff3`
- `custom_with_common_snps.gbk`

## Requirements

System packages:
- `bowtie`
- `gawk`

Python packages:
- `streamlit`
- `biopython`
- `bcbio-gff`
- `dna_features_viewer`

## Main files

- `streamlit_app.py` — main app
- `make_between_regions.py` — tile prediction helper
- `gff_to_genbank_patched.py` — GFF3/FASTA to GenBank helper

## Notes

- hg38 mode supports canonical human chromosomes only
- gene lookups use Ensembl
- sequence and common SNP retrieval use UCSC
- only gene-level annotations are included from Ensembl
