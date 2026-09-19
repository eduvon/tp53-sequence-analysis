# TP53 Variant Analysis

A bioinformatics project investigating the relationship between genetic
variation in human TP53, predicted molecular consequences, and clinical
significance using Python, ClinVar, Ensembl VEP, and command-line
genomics tools.

## Research Question

How do molecular consequences of TP53 genetic variants relate to their
clinical significance in ClinVar?

TP53 is used as a case study for building an end-to-end computational
biology workflow from public genomic data through variant annotation,
statistical analysis, and biological interpretation.

## Project Goals

This project is being developed as a hands-on study of computational
biology and bioinformatics, with an emphasis on:

- Working with real genomic data
- Understanding relationships between genes, transcripts, and proteins
- Building reproducible genomic analysis workflows
- Interpreting variant annotations biologically
- Applying statistical methods to biological data

## Current Work

- Basic DNA sequence analysis
- Reverse-complement analysis
- FASTA and sequence analysis
- TP53 transcript and coding-sequence analysis
- DNA-to-protein translation using Biopython
- Retrieval of TP53-region variants from ClinVar
- Parsing VCF records with Python
- Analysis of ClinVar clinical significance classifications
- Genomic distribution of TP53 variants
- Variant annotation with Ensembl VEP
- Transcript-level consequence analysis
- MANE Select transcript analysis
- Variant-level consequence classification
- Statistical association analysis using Fisher's exact test
- Odds ratios and 95% confidence intervals
- Benjamini-Hochberg false-discovery-rate correction
- Forest-plot visualization of consequence associations

## Data

Clinical variant data are obtained from ClinVar.

VEP annotations are generated using Ensembl Variant Effect Predictor
against the GRCh38 reference assembly.

The analysis considers both:

- All TP53 transcript annotations
- The TP53 MANE Select transcript (`NM_000546.6`)

The large ClinVar VCF and VEP output files are not included in this
repository. See the project documentation for information about
obtaining or generating the data.

## Analysis

The primary analysis focuses on four molecular consequence categories:

- Missense
- Frameshift
- Stop gained
- Splice site

These categories were selected because they represent molecular
changes that can directly affect protein sequence or RNA processing.

For the main statistical analysis, variants classified by ClinVar as
Pathogenic or Uncertain are compared with respect to the presence or
absence of each consequence type on the MANE Select transcript.

### Preliminary Results

ClinVar returned 3,876 records in the TP53 genomic region. One record
represented a `no_sequence_alteration` assertion rather than a
sequence-changing variant, leaving 3,875 variants annotated by VEP.

The preliminary statistical analysis shows different consequence
profiles between the ClinVar Pathogenic and Uncertain groups. Frameshift,
stop-gained, and splice-site consequences were associated with higher
odds in the Pathogenic group, while missense consequences were
associated with higher odds in the Uncertain group.

All four associations remained statistically significant after
Benjamini-Hochberg adjustment for multiple comparisons.

## Statistical Results

Odds ratios compare the odds of each consequence in ClinVar Pathogenic
variants with the odds in Uncertain variants.

![TP53 consequence odds ratios](results/figures/tp53_consequence_odds_ratios.png)

The points represent estimated odds ratios, horizontal lines represent
95% confidence intervals, and the vertical reference line at OR = 1
represents no association.

## Tools

- Python
- Biopython
- pandas
- SciPy
- statsmodels
- Matplotlib
- tabix / htslib
- Docker
- ClinVar
- Ensembl VEP
- MANE

## Project Status

The project is currently under development.

### Completed

- TP53 sequence and transcript characterization
- ClinVar variant extraction and parsing
- VEP variant annotation
- Transcript-aware consequence analysis
- MANE Select analysis
- Variant-level consequence classification
- Statistical association testing
- Odds-ratio visualization

### Next Steps

- Biological interpretation of the observed consequence patterns
- Protein-level analysis of TP53 variants
- Additional validation and quality checks
- Further visualization
- Final documentation and reproducibility improvements