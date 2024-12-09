# collapse-bubble

SV merging for pangenome VCF generated by `minigraph-cactus` pipeline and decomposed by `vcfwave`. It merges similar SVs **within** the same bubble, e.g., highly polymorphic INS at the same position.

- [Overview](#overview)
- [Dependency](#dependency)
- [Input](#input)
- [Analysis pipeline](#analysis-pipeline)
    1. [Preprocessing](#1-preprocessing)
    2. [SV merging](#2-sv-merging)
    3. [Postprocessing](#3-postprocessing)

**Note**: I will consolidate the entire pipeline into a single script once I have time. For now, please follow the step-by-step instructions.

## Overview

This pipeline uses [Truvari](https://github.com/ACEnglish/truvari)'s API to merge SVs and genotypes. As compared to `truvari collapse`, it is optimized for pangenome VCF:

1. Only SVs within the same bubble are compared for merging.
2. When merging genotypes, the phase of genotypes is checked and retained:
    - 0|0 and 1|0: consistent, merge into 1|0
    - 0|1 and 1|0: consistent, merge into 1|1
    - 1|0 and 1|0: inconsistent, don't merge if found in any sample
3. If a SV has large REF and ALT alleles, its size is determined by `REFLEN` (length of REF) instead of `SVLEN` (length difference between REF and ALT). This works better for large inversions or complex SVs. For example, 100bp INV and 1000bp INV have the same `SVLEN` but different `REFLEN`

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
- [vcfwave](https://github.com/vcflib/vcflib/blob/master/doc/vcfwave.md): please use `vcfwave` v1.0.12 or later, the previous versions can output incorrect genotypes for a few variants.

## Input
- Multiallelic graph VCF generated by [minigraph-cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) and processed by [vcfbub](https://github.com/pangenome/vcfbub). The ID field of this VCF should be the unique bubble ID, e.g. `<73488<73517` . This is default output VCF of `minigraph-cactus` pipeline.
- Reference genome fasta file for normalization

## Analysis pipeline

### 1. Preprocessing:

The preprocessing step aims to generate a VCF file meet the following requirements:

- decomposed by `vcfwave`
- biallelic, sorted VCF
- all variants have unique ID
- only variants from the same bubble have identical `INFO/BUBBLE_ID`
- AC/AN/AF have been updated based on the genotypes

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

This script uses `Truvari` API to merge SVs have the same `INFO/BUBBLE_ID`, it output 2 files:

**1. VCF**:

The output VCF is not sorted and has merged variants and genotypes. Similar SVs are merged into the one with the highest `MAF`. The following fields of the merged SV are added or modified:

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

The merged VCF is not sorted, so we need to sort it. Moreover, variants from differnt bubbles can have the same position and alleles after left-align. We can use `merge_duplicates.py` to further merge duplicated variants. This step can significantly improve the genotype concordance between variants called from pangenome graph and linear reference genome.

```
# sort VCF
bcftools sort merged.vcf.gz -Oz -o merged.sort.vcf.gz

# merge duplicated variants with the same POS, REF, ALT
python merge_duplicates.py \
-i merged.sort.vcf.gz \
-o merged.sort.dedup.vcf.gz
```

If there are two identical alleles on the same haplotype, they will be concatenated to create a new variant. For example, two insertions (A → AAA) in sample 2 result in a longer insertion (A → AAAAA).


Input:
```
#CHROM   POS   ID     REF   ALT     CHM13   Sample1   Sample2   Sample3
chr1     5     var1   A     AAA     0       0|0       1|0       0|1
chr1     5     var2   A     AAA     1       1|.       0|0       0|0
chr1     5     var3   A     AAA     .       0|0       1|0       0|0
```

Output:

```
#CHROM   POS   ID     REF   ALT     CHM13   Sample1   Sample2   Sample3
chr1     5     var1   A     AAA     1       1|.       0|0       0|1
chr1     5     var3   A     AAAAA   .       0|0       1|0       0|0
```

Merging missing genotypes with non-missing ones can also yield a missing genotype in one of the variants, depending on whether the merge is with reference or alternate alleles. If the option `--merge-mis-as-ref`` is used, missing genotypes will be treated as reference alleles during the merging process. Merging two missing genotypes will always result in a missing genotype.

**Arguments**:

```
usage: merge_duplicates.py [-h] -i VCF -o VCF [--merge-mis-as-ref]

options:
  -i VCF, --invcf VCF   Input VCF, sorted and phased
  -o VCF, --outvcf VCF  Output VCF
  --merge-mis-as-ref    Convert missing to ref when merging missing genotypes with non-missing genotypes
```
