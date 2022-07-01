# BPHunter
Genome-wide detection of human variants that disrupt intronic branchpoints

## Introduction
- The search for pathogenic candidate variants in next-generation sequencing (NGS) data typically focuses on non-synonymous variants within coding sequence or variants in essential splice sites, while mostly ignoring non-coding intronic variants. 

- RNA splicing, as a necessary step for protein-coding gene expression in eukaryotic cells, operates its spliceosome mostly within introns to define the exon-intron boundaries and hence the coding sequences. Introns probably harbor a substantially larger number of pathogenic mutations than has so far been appreciated. 

- Intronic branchpoint (BP) is recognized by spliceosome in the beginning of the splicing process, and constitutes a vulnerability of splicing by its potential mutations. BP mutations may potentially result in aberrant splicing consequences (exon skipping, intron retention), which could be deleterious to the gene product.

- BPHunter is a genome-wide computational approach to systematically detect intronic variants that may disrupt BP recognition in NGS data, in an efficient, systematic and informative manner. Its single-line command that can be easily implemented into NGS analysis. We also provided a [BPHunter webserver](http://hgidsoft.rockefeller.edu/BPHunter) for users with less computational expertise.

## Usage
### Dependency
The code is written in [python3](https://www.python.org/downloads/), and requires [bedtools](https://bedtools.readthedocs.io/en/latest/) installed.

### Reference datasets
Due to the file size is limited at max. 25MB in GitHub, please download the [BPHunter reference datasets](http://hgidsoft.rockefeller.edu/BPHunter/standalone.html) and put them into your BPHunter folder.

### File Format
**Input:** Variants in VCF format, with 5 mandatory and tab-delimited fields (CHROM, POS, ID, REF, ALT).
  - The 48 published pathogenic mutations are provided as the example of input data. (Data_BPMut.vcf)

**Output:** Variants detected by BPHunter that may disrupt BP thus splicing, in a file with suffix '.bphunter.txt'. The output includes the following annotation.
  - CHROM, POS, ID, REF, ALT (exactly the same as input)
  - STRAND, VAR_TYPE, GENE, BP_NAME, BP_RANK, HIT_POS, BP_3SS_DIST, CONSENSUS, EVI, SOURCE_LIST, MAF, GERP, PHYLOP, IVS_TYPE, IVS_LENGTH, TRANSCRIPT_#IVS, BPHunter_SCORE

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
*-i*|file|variants in VCF format, with 5 fields (CHROM, POS, ID, REF, ALT)|N.A.
*-g*|str|human reference genome assembly (hg19/hg38)|hg19
*-c*|str|canonical transcripts? (no/yes)|no

### BPHunter Scoring Scheme
<img src="https://hgidsoft.rockefeller.edu/BPHunter/img/BPHunter_Scoring.png" width="70%" height="70%">
<!--START_SECTION:update_image-->
<!--END_SECTION:update_image-->

## References
- *Zhang P. et al.* Genome-wide detection of human variants that disrupt intronic branchpoints. (2022)

## Contact
> **Developer:** Peng Zhang, Ph.D.

> **Email:** pzhang@rockefeller.edu

> **Laboratory:** St. Giles Laboratory of Human Genetics of Infectious Diseases

> **Institution:** The Rockefeller University, New York, NY, USA
