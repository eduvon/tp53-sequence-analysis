# TP53 Sequence Analysis

A bioinformatics project analyzing human TP53 sequence and clinical
variants using Python, Biopython, ClinVar, and command-line genomics tools.

## Project Goals

This project explores the relationship between genetic variation in the
TP53 gene and its biological and clinical consequences.

The project is being developed as a hands-on study of bioinformatics,
with an emphasis on working with real genomic data and building
reproducible analysis workflows.

## Current Work

- Basic DNA sequence analysis
- Reverse-complement analysis
- Basic FASTA and sequence analysis
- TP53 transcript and coding-sequence analysis
- DNA-to-protein translation using Biopython
- Retrieval of TP53 variants from a ClinVar VCF
- Parsing VCF records with Python
- Analysis of ClinVar clinical significance classifications
- Genomic distribution of TP53 variants

## Data

Clinical variant data are obtained from ClinVar.

The large ClinVar VCF files are not included in this repository.
See the project documentation for information about obtaining the data.

## Tools

- Python
- Biopython
- pandas
- Matplotlib
- tabix / htslib
- ClinVar

## Project Status

This project is currently under development.

The next stage will annotate TP53 variants with transcript-level and
protein-level consequences and investigate the biological distribution
of clinically significant variants.