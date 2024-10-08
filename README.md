# collapse-bubble

(under development)

SV merging for pangenome VCF generated by `minigraph-cactus` pipeline and decomposed by `vcfwave`. It uses `Truvari`'s API to merge similar SVs **within** the same bubble, e.g., highly polymorphic INS at the same position.


## Overview
Input:
- multiallelic graph VCF generated by `minigraph-cactus` pipeline and processed by `vcfbub`

Preprocess:
- split into biallelic: `bcftools norm -m -any`
- annotate VCF and assign unique variant ID: `scripts/annotate_var_id.py`
- decompose VCF: `vcfwave`
- annotate, normalize and sort VCF: `bcftools +fill-tags wave.vcf.gz -- -t AC,AN,AF | bcftools norm -f ref.fa | bcftools sort`

SV merging (`scripts/collapse_bubble.py`)

1. Merged VCF (unsorted):
- `INFO/ID_LIST`: list of collapsed SVs
- `INFO/TYPE`: type of variant (SNP, MNP, INS, DEL, INV, COMPLEX)
- `INFO/REFLEN`: len(ref)
- `INFO/SVLEN`: len(alt) - len(ref)
- `FORMAT/GT`: merged genotypes

2. SV merging results in tsv format, e.g.:
```
CHROM   POS      Bubble_ID   Variant_ID             Collapse_ID
chr1    10039    >123>456    >123>456.COMPLEX.1_1   >123>456.COMPLEX.1_1
chr1    10040    >123>456    >123>456.COMPLEX.1_2   >123>456.COMPLEX.1_1
```

Postprocess:
- sort VCF: `bcftools sort merged.vcf.gz`
- merge duplicated variants: `scripts/merge_duplicates.py`


## Dependency
```
truvari
pysam
pandas
numpy
```