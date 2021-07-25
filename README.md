# BPHunter
Genome-wide detection of human variants that disrupt intronic branchpoints *(LICENSE: CC BY-NC-ND 4.0)*

## Introduction
- The search for pathogenic candidate variants in next-generation sequencing (NGS) data typically focuses on non-synonymous variants within coding sequence or variants in essential splice sites, while mostly ignoring other non-coding intronic variants. 

- RNA splicing, as a necessary step for protein-coding gene expression in eukaryotic cells, operates its spliceosome mostly within introns to define the exon-intron boundaries and hence the coding sequences. Introns probably harbor a substantially larger number of pathogenic mutations than has so far been appreciated. 

- Intronic branchpoint (BP) is recognized by the spliceosome at the beginning of the splicing process, and constitutes a vulnerability of splicing by its potential mutations. BP mutations may potentially result in aberrant splicing consequences (exon skiping, intron retention), which could be deleterious.

- We developed BPHunter as a genome-wide computational approach to systematically detect intronic variants that may disrupt BP recognition in NGS data, in an efficient and informative manner. Its one-line command that can be easily implemented into NGS analysis. We also provided [BPHunter webserver](http://hgidsoft.rockefeller.edu/BPHunter) for users with less computational expertise.

## Usage
### Dependency
The code is written in [python3](https://www.python.org/downloads/), and requires [bedtools](https://bedtools.readthedocs.io/en/latest/) installed.

### Reference datasets
Due to the file size is limited at max. 25MB in GitHub, please download the [BPHunter reference datasets](http://hgidsoft.rockefeller.edu/BPHunter/reference_datasets.html) and put them into your BPHunter folder.

### File Format
**Input:** Variants in VCF format, with 5 mandatory and tab-delimited fields (CHROM, POS, ID, REF, ALT), where ID field will be ignored in running BPHunter.
  - The 44 published pathogenic mutations are provided as the example of input data. *(Data_BPMut.vcf)*

**Output:** Variants detected by BPHunter that may disrupt BP thus splicing, with following annotation, in a file with surfix '.bphunter.out':
  - CHROM, POS, REF, ALT, STRAND, VAR_TYPE
  - GENE, BP_NAME, BP_TYPE, BP_RANK, HIT_POS, DIST_3SS
  - MAF, GERP, PHYLOP, ENERGY, CONSENSUS, EV_SCORE, SOURCE
  - INTRON_TYPE, TRANSCRIPT_#INTRON

### Command
```
python BPHunter.py -i variants.vcf
```
```
python BPHunter.py -i variants.vcf -g hg19/hg38 -c no/yes
```

### Parameters
Parameter | Type | Description | Default
----------|------|-------------|--------------
*-i*|file|input variants in VCF file|N.A.
*-g*|str|human reference genome assembly (hg19/hg38)|hg19
*-c*|str|canonical transcripts? (no/yes)|no

## References
- *Zhang P. et al.* Genome-wide detection of human variants that disrupt intronic branchpoints. (2021)

## Contact
> **Author:** Peng Zhang, Ph.D.

> **Email:** pzhang@rockefeller.edu

> **Laboratory:** St. Giles Laboratory of Human Genetics of Infectious Diseases

> **Institution:** The Rockefeller University, New York, NY, USA
