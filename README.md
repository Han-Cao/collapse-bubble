# collapse-bubble

(under development)

SV merging for pangenome VCF generated by `minigraph-cactus` pipeline and decomposed by `vcfwave`. It merge similar SVs **within** the same bubble, e.g., highly polymorphic INS at the same position.

- [Overview](#overview)
- [Dependency](#dependency)
- [Input](#input)
- [Analysis pipeline](#analysis-pipeline)
    1. [Preprocessing](#1-preprocessing)
    2. [SV mering](#2-sv-merging)
    3. [Postprocessing](#3-postprocessing)


## Overview

This pipeline uses [Truvari](https://github.com/ACEnglish/truvari)'s API to merge SVs and genotypes. As compared to `truvari collapse`, it is optimized for pangenome VCF:

1. Only SVs within the same bubble are compared for merging.
2. When merging genotypes, the phase of genotypes is checked and retained:
    - 0|0 and 1|0: consistent, merge into 1|0
    - 0|1 and 1|0: consistent, merge into 1|1
    - 1|0 and 1|0: inconsistent, don't merge if found in any sample
3. If a SV has large REF and ALT alleles, its size is determined by `REFLEN` (length of REF) instead of `SVLEN` (length difference between REF and ALT). This works better for large inversions or complex SVs. For example, 100bp INV and 1000bp INV have the same `SVLEN` but differnt `REFLEN`

## Dependency

All scripts have been tested with Python 3.10. Please install the following python modules to run the scripts.

```
truvari
pysam
pandas
numpy
```

We also use the below tools to process VCF:
- [bcftools](https://github.com/samtools/bcftools): `+fill-tags` plugin is used, please export `BCFTOOLS_PLUGINS` to configure it
- [vcfwave](https://github.com/vcflib/vcflib/blob/master/doc/vcfwave.md): version 1.0.10 is used

## Input
- Multiallelic graph VCF generated by [minigraph-cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) and processed by [vcfbub](https://github.com/pangenome/vcfbub). The ID field of this VCF should be the unique bubble ID, e.g. `<73488<73517` . This is default output VCF of `minigraph-cactus` pipeline.
- Reference genome fasta file for normalization

## Analysis pipeline

### 1. preprocessing:

The preprocessing step aims to generate a VCF file based on the following criteria:

- decomposed by `vcfwave`
- biallelic, sorted VCF
- all variants have unique ID
- only variants from the same bubble have identical `INFO/BUBBLE_ID`
- AC/AN/AF align with the genotypes

If you already have such a VCF, it can be directly used without preprocessing. To start with the multiallelic graph VCF, we will do:
- split multiallelic records into biallelic, update AC, AN, AF (`bcftools`)
- annotate variants' bubble ID, rename variants with unique names (`annotate_var_id.py`)
- decompose VCF (`vcfwave`)
- left-align and sort (`bcftools`)

```
# suppose the name of input VCF is "mc.vcf.gz"

# split into biallelic
bcftools norm -m -any mc.vcf.gz -Oz -o mc.biallele.vcf.gz

# annotate VCF and assign unique variant ID
python scripts/annotate_var_id.py \
-i mc.biallele.vcf.gz \
-o mc.biallele.uniq_id.vcf.gz

# add AC/AN/AF and decompose VCF
# optional: first drop INFO/AT (useless after decomposition, suggested by cactus)
bcftools annotate -x INFO/AT -Ou mc.biallele.uniq_id.vcf.gz | \
bcftools +fill-tags -- -t AC,AN,AF | \
vcfwave -I 1000 | bgzip -c > mc.biallele.uniq_id.vcfwave.vcf.gz

# normalization and sort, update AC/AN/AF
# note: vcfwave likely report incorrect AC/AN/AF, we re-calculate it here
bcftools norm -f ref.fa mc.biallele.uniq_id.vcfwave.vcf.gz | \
bcftools +fill-tags -- -t AC,AN,AF | \
bcftools sort --max-mem 4G -Oz -o mc.biallele.uniq_id.vcfwave.sort.vcf.gz
```

`annotate_var_id.py`:

After spliting the multiallelic graph VCF into a biallelic VCF, all variants from the same bubble (i.e., same multiallelic record) should have the same value in ID field.This script copies this ID to `INFO/BUBBLE_ID`, and renames variants with a unique ID in format of [BUBBLE_ID].[TYPE].[No.].

```
usage: annotate_var_id.py [-h] -i VCF -o VCF

options:
  -i VCF, --input VCF   Input VCF
  -o VCF, --output VCF  Output VCF
```

### 2. SV merging 

To perform SV merging, run `collapse_bubble.py`:
```
python collapse_bubble.py \
-i mc.biallele.uniq_id.vcfwave.sort.vcf.gz \
-o merged.vcf.gz \
--map merged.mapping.txt
```

This `Truvari` API to merge SVs have the same `INFO/BUBBLE_ID`, it output 2 files:

**1. VCF**:

The output VCF is not sorted and has merged variants and genotypes. Similar SVs are merged into the one with highest `MAF`. The following fields of the merged SV are added or modified:

- `INFO/ID_LIST`: comma-separated list of variants merged into this variants
- `INFO/TYPE`: type of variant (SNP, MNP, INS, DEL, INV, COMPLEX)
- `INFO/REFLEN`: len(ref)
- `INFO/SVLEN`: len(alt) - len(ref)
- `FORMAT/GT`: merged genotypes

For example:
```
#input:
chr1  10039  >123>456.INS.1  A  ATTTTTT  AC=2;AN=6;AF=0.333;BUBBLE_ID=>123>456  0|1  1|0  0|0
chr1  10039  >123>456.INS.2  A  ATTTTTG  AC=3;AN=6;AF=0.500;BUBBLE_ID=>123>456  1|0  0|1  1|0

#output
chr1  10039  >123>456.INS.2  A  ATTTTTG  AC=3;AN=6;AF=0.500;BUBBLE_ID=>123>456;ID_LIST=>123>456.INS.1;TYPE=INS;REFLEN=1;SVLEN=6  1|1  1|1  1|0
```

**2. SV merging table**:

A tsv file to map the original SVs (`Variant_ID`) to merged SVs (`Collapse_ID`), for example:
   
```
CHROM   POS      Bubble_ID   Variant_ID       Collapse_ID
chr1    10039    >123>456    >123>456.INS.1   >123>456.INS.2
chr1    10039    >123>456    >123>456.INS.2   >123>456.INS.2
```

Moreover, VCF INFO fields can also be added as columns by specifying `--info`. If `--info SVLEN` is used, the output SVLEN in the tsv file will be REFLEN for COMPLEX and INV to indicate the value used for comparison. While in the output VCF, INFO/SVLEN is always set as `len(alt) - len(ref)` for all SVs.

**Arguments**:

```
usage: collapse_bubble.py [-h] -i VCF -o VCF -m TSV [--chr CHR] [--info TAG] [-l 50] [-r 100] [-p 0.9] [-P 0.9] [-O 0.9]

Input / Output arguments:
  -i VCF, --invcf VCF   Input VCF
  -o VCF, --outvcf VCF  Output VCF
  -m TSV, --map TSV     Write SV mapping table to this file. Default: None
  --chr CHR             chromosome to work on. Default: all
  --info TAG            Comma-separated INFO/TAG list to include in the output map. Default: None

Collapse arguments:
  -l 50, --min-len 50   Minimum allele length of variants to be included, 
                        defined as max(len(alt), len(ref)). Default: 50
  -r 100, --refdist 100
                        Max reference location distance. Default: 100
  -p 0.9, --pctseq 0.9  Min percent sequence similarity (REF for DEL, ALT for other SVs). Default: 0.9
  -P 0.9, --pctsize 0.9
                        Min percent size similarity (SVLEN for INS, DEL; REFLEN for INV, COMPLEX). Default: 0.9
  -O 0.9, --pctovl 0.9  Min pct reciprocal overlap. Default: 0.9
```

### 3. Postprocessing:

The merged VCF is not sorted, so we first sort it by `bcftools sort merged.vcf.gz -Oz -o merged.sort.vcf.gz`. Moreover, INDELs from differnt bubbles can have the same position and alleles after left-align. We can use `merge_duplicates.py` to further merge duplicated variants:

```
merge_duplicates.py \
-i merged.sort.vcf.gz \
-o merged.sort.dedup.vcf.gz
```

If 2 variants are truly duplicated, we anticipate that they should not have any genotype conflict on the same haplotype (e.g., 1|0 vs 1|0). Therefore, duplicated variants with haplotype conflict will not be merged and both of them will write to output by default. If `--keep first` is specified, only the first variant will write to output (like `bcftools norm --rm-dup exact`). Moreover, if 2 duplicated variants have both missing and non-missing genotypes on the same haplotype (e.g., .|0 vs 1|0), they can be merged with warning (default), without warning, or not merged just like haplotype conflicts. This can be adjusted by setting `--missing warn|merge|conflict`.

**Arguments**:

```
usage: merge_duplicates.py [-h] -i VCF -o VCF [--keep all|first] [--missing warn|merge|conflict]

options:
  -i VCF, --invcf VCF   Input VCF, sorted and phased
  -o VCF, --outvcf VCF  Output VCF
  --keep all|first      For duplicated variants with haplotype conflict (e.g., 1|0 vs 1|0), 
                        output all of them (all) or only the first one (first).
  --missing warn|merge|conflict
                        For duplicated variants with missing conflict (e.g., .|0 vs 1|0), 
                        merge them with warning (warn), without warning (merge), or treat as haplotype 
                        conflict (conflict).
```
