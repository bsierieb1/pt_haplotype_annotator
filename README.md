# PureTarget haplotype annotator
### not an official PacBio app; very early version; use at your own risk!

Uploads:
- haplotype FASTA
- region annotations BED4
- "even" set of guides FASTA
- "odd" set of guides FASTA

Outputs:
- GFF3 file with all annotations combined
- GenBank (.gbk) file with the sequence and all annotations combined

System deps (packages.txt): bowtie, gawk  
Python deps (requirements.txt): streamlit, gff-to-genbank, biopython==1.81
