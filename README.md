# collapse-bubble

[![pytest](https://github.com/Han-Cao/collapse-bubble/actions/workflows/pytest.yml/badge.svg)](https://github.com/Han-Cao/collapse-bubble/actions)
[![Coverage](https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/han-cao/5baaf521e0e161712d3493546f6a8876/raw/collapse-bubble-cobertura-coverage.json)](https://github.com/Han-Cao/collapse-bubble)

Variant merging for VCF deconstructed from a pangenome graph, including vertical merging for overlapping/duplicated variants and horizontal merging for similar structural variantions (SVs).

- [Overview](#overview)
- [Dependency](#dependency)
- [Input](#input)
- [Analysis pipeline](#analysis-pipeline)
    1. [Preprocessing](#1-preprocessing)
    2. [Merge overlapping variants](#2-merge-overlapping-variants)
    3. [SV merging](#3-sv-merging)
- [Todo list](#todo-list)
- [Acknowledgement](#acknowledgement)

## Overview

**Vertical variant merging**

Variant decomposition and normalization are usually used to simplify variants deconstructed from a pangenome graph and make variants comparable between different callsets. However, these analysis can produce overlapping or even duplicated variant records. The `collapse-bubble` pipeline can concatenate overlapping variants and merge the genotypes of duplicated records to generate a deduplicated non-overlapping VCF.

**Horizontal SV merging**

A pangenome VCF includes highly similar SVs with even 1 base difference. Therefore, SV merging is important to remove redundant SVs. The `collapse-bubble` pipeline uses [Truvari](https://github.com/ACEnglish/truvari)'s engine to merge SVs. As compared to `truvari collapse`, it is optimized for pangenome VCF by using bubbles and haplotypes to improve SV merging.

For more information on how it works, please refer to the [documentation](docs/).

**Note**: `collapse-bubble`  is currently tested with VCFs from the [minigraph-cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) pipeline. But it should theoretically support VCFs generated by [vg deconstruct](https://github.com/vgteam/vg/wiki/VCF-export-with-vg-deconstruct) from phased assembly-based pangenome graphs.

## Dependency

All scripts have been tested with Python 3.10. Please run the following command to install the required modules:

```
pip install -r requirements.txt
```

The following tools are not used by `collapse-bubble` scripts but are required to prepare the input VCF:
- [bcftools](https://github.com/samtools/bcftools): The `+fill-tags` plugin is used. Please make sure `BCFTOOLS_PLUGINS` is configured.
- [vcfwave](https://github.com/vcflib/vcflib/blob/master/doc/vcfwave.md): Please use `vcfwave` v1.0.12 or later, as earlier versions may output incorrect genotypes for some mutli-allelic variants.

## Input
- **Multiallelic graph VCF**: Generated by the [vg deconstruct](https://github.com/vgteam/vg/wiki/VCF-export-with-vg-deconstruct) and processed by [vcfbub](https://github.com/pangenome/vcfbub). The ID field of this VCF are bubble IDs, e.g., `<73488<73517`. This is the default output VCF of the [minigraph-cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) pipeline.
- **Reference genome FASTA file**: Used for normalization.

## Analysis pipeline

Example script to run the pipeline from the default output VCF of `minigraph-cactus`:
```
##### 1. Preprocessing #####
# split into biallelic
bcftools norm -m -any mc.vcf.gz -Oz -o mc.biallele.vcf.gz

# annotate VCF and assign unique variant ID
python scripts/annotate_var_id.py \
-i mc.biallele.vcf.gz \
-o mc.biallele.uniq_id.vcf.gz

# drop INFO/AT and decompose by vcfwave
bcftools annotate -x INFO/AT mc.biallele.uniq_id.vcf.gz | \
vcfwave -I 1000 | bgzip -c > mc.biallele.uniq_id.vcfwave.vcf.gz

# normalize variants, update AC/AN/AF, and sort
bcftools norm -f ref.fa mc.biallele.uniq_id.vcfwave.vcf.gz | \
bcftools +fill-tags -- -t AC,AN,AF | \
bcftools sort --max-mem 4G -Oz -o mc.biallele.uniq_id.vcfwave.sort.vcf.gz

##### 2. Merge overlapping variants #####
python scripts/merge_duplicates.py \
-i mc.biallele.uniq_id.vcfwave.sort.vcf.gz \
-o mc.biallele.uniq_id.vcfwave.sort.merge_dup.vcf.gz \
-c repeat \
-t ID

##### 3. SV merging #####
python collapse_bubble.py \
-i mc.biallele.uniq_id.vcfwave.sort.merge_dup.vcf.gz \
-o mc.biallele.uniq_id.vcfwave.sort.merge_dup.merge_sv.vcf.gz \
--map mc.biallele.uniq_id.vcfwave.sort.merge_dup.merge_sv.mapping
```

### 1. Preprocessing:

The preprocessing step generates a VCF file that meets the following requirements:
- Decomposed by `vcfwave`.
- Biallelic and sorted.
- All variants have unique IDs.
- Variants from the same bubble have identical `INFO/BUBBLE_ID` annotation.
- AC, AN, and AF have been updated based on the genotypes.


If you already have such a VCF, it can be used directly without preprocessing. To start with the multiallelic graph VCF, the following steps are required:

1. Split multiallelic records into biallelic and update AC, AN, AF using `bcftools`.
2. Annotate variants' bubble ID and generate unique variant IDs using `annotate_var_id.py`.
3. Decompose the VCF using `vcfwave`.
4. Left-align and sort using `bcftools`.

```
# suppose the name of input VCF is "mc.vcf.gz"

# split into biallelic
bcftools norm -m -any mc.vcf.gz -Oz -o mc.biallele.vcf.gz

# annotate VCF and assign unique variant ID
python scripts/annotate_var_id.py \
-i mc.biallele.vcf.gz \
-o mc.biallele.uniq_id.vcf.gz

# drop INFO/AT (optional, suggested by cactus) and decompose by vcfwave
bcftools annotate -x INFO/AT mc.biallele.uniq_id.vcf.gz | \
vcfwave -I 1000 | bgzip -c > mc.biallele.uniq_id.vcfwave.vcf.gz

# for merge_duplicates.py -c repeat:
# fast normalization, update AC/AN/AF, and sort in one step
bcftools norm -f ref.fa mc.biallele.uniq_id.vcfwave.vcf.gz | \
bcftools +fill-tags -- -t AC,AN,AF | \
bcftools sort --max-mem 4G -Oz -o mc.biallele.uniq_id.vcfwave.sort.vcf.gz

# for merge_duplicates.py -c position:
# normalization, update AC/AN/AF, and stable sort
bcftools norm -f ref.fa mc.biallele.uniq_id.vcfwave.vcf.gz | \
bcftools +fill-tags -Oz -o tmp.vcf.gz -- -t AC,AN,AF
(bcftools view -h tmp.vcf.gz ; bcftools view -H tmp.vcf.gz | sort -s -k1,1d -k2,2n) | bgzip > mc.biallele.uniq_id.vcfwave.sort.vcf.gz
```

`annotate_var_id.py`:

This script assign unique IDs in format of [BUBBLE_ID].[TYPE].[No.] to each variants. The original variant ID (i.e., bubble ID) is stored in `INFO/BUBBLE_ID`. If the VCF has been processed by `vcfwave`, the separator (`_`) between bubble ID and suffix can be customized by `--suffix-sep _`.


```
usage: annotate_var_id.py [-h] -i VCF -o VCF [--suffix-sep SUFFIX_SEP]

Annotate and assign unique variant ID for pangenome VCF

options:
  -i VCF, --input VCF   Input VCF
  -o VCF, --output VCF  Output VCF
  --suffix-sep SUFFIX_SEP
                        Separator between bubble ID and suffix, e.g., "_" for vcfwave processed VCF (default: None)
```


### 2. Merge overlapping variants

After variant decomposition and left align, the VCF contains overlapping variants at the same position (mostly SNPs and small INDELs). For example:

```
chr1   100   var1   C   G     1|0
chr1   100   var2   C   G     0|1
chr1   100   var3   C   CAA   1|0
chr1   100   var4   C   CAA   1|0
```
In this example:
- `var1` and `var2` are duplicates, as they share the same `POS`, `REF`, and `ALT`. This is mainly due to variant decomposition.
- `var1` and `var3`/`var4` overlap on the **first haplotype**, as there are three alternative alleles at the same `POS`. This is caused by left align.

`merge_duplicates.py` can clean up duplicated and overlapping variants:

```
python scripts/merge_duplicates.py \
-i mc.biallele.uniq_id.vcfwave.sort.vcf.gz \
-o mc.biallele.uniq_id.vcfwave.sort.merge_dup.vcf.gz \
-c repeat \
-t ID
```

- It first concatenates overlapping **tandem repeats** (specified by `-c repeat`) using the algorithm described in the [documentation](docs/merge_duplicates.md). For example, `var3` and `var4` are concatenated into `C CAAAA`.
- After concatenating all overlapping variants, it merges duplicates into a single record and also updates the phased genotypes.
- `-t ID` tracks how the overlapping variants are concatenated (`INFO/CONCAT`) and how duplicates are merged (`INFO/DUP`). These information is useful for downstream SV merging.

Output:
```
chr1   100   var1   C   G     1|1
chr1   100   chr1:100_0   C   CAAAA 1|0
```

**Note**: When using `merge_duplicates.py -c position`, it concatenates any overlapping variants at the same position. This method reconstructs the local haplotypes and significantly increasing polymorphism, which may not be suitable for merging tandem repeats. Additionally, it requires the input VCF sorted by `CHROM` and `POS` only (not guaranteed by recent `bcftools`). Please see [documentation](docs/merge_duplicates.md) for more details.


**Arguments**:
```
usage: merge_duplicates.py [-h] -i VCF -o VCF [-c {position,repeat,none}] [-m MAX_REPEAT] [-t {ID,AT}] [--merge-mis-as-ref] [--keep-order] [--debug]

Merge duplicated variants in phased VCF

options:
  -i VCF, --invcf VCF   Input VCF, sorted and phased
  -o VCF, --outvcf VCF  Output VCF
  -c {position,repeat,none}, --concat {position,repeat,none}
                        Concatenate variants when they have identical "position" (default) or "repeat" motif, "none" to skip
  -m MAX_REPEAT, --max-repeat MAX_REPEAT
                        Maximum size a variant to search for repeat motif (default: None)
  -t {ID,AT}, --track {ID,AT}
                        Track how variants are merged by "ID" or "AT" (default: None)
  --merge-mis-as-ref    Convert missing to ref when merging missing genotypes with non-missing genotypes
  --keep-order          keep the order of variants in the input VCF (default: sort by chr, pos, alleles)
  --debug               Debug mode
```

### 3. SV merging 

To perform SV merging, run `collapse_bubble.py`:
```
python collapse_bubble.py \
-i mc.biallele.uniq_id.vcfwave.sort.merge_dup.vcf.gz \
-o mc.biallele.uniq_id.vcfwave.sort.merge_dup.merge_sv.vcf.gz \
--map mc.biallele.uniq_id.vcfwave.sort.merge_dup.merge_sv.mapping
```

This will generate 3 output files:

**1. VCF**:

The output VCF is not sorted and has merged variants and genotypes. Similar SVs are merged into the most common one. The following fields of the merged SV are added:

- `INFO/ID_LIST`: comma-separated list of variants merged into this variants
- `INFO/TYPE`: type of variant (SNP, MNP, INS, DEL, INV, COMPLEX)
- `INFO/REFLEN`: len(ref)
- `INFO/SVLEN`: len(alt) - len(ref)

For example:
```
#input:
chr1  10039  >123>456.INS.1  A  ATTTTTT  AC=2;AN=6;AF=0.333;BUBBLE_ID=>123>456  0|1  1|0  0|0
chr1  10039  >123>456.INS.2  A  ATTTTTG  AC=3;AN=6;AF=0.500;BUBBLE_ID=>123>456  1|0  0|1  1|0

#output
chr1  10039  >123>456.INS.2  A  ATTTTTG  AC=5;AN=6;AF=0.833;BUBBLE_ID=>123>456;ID_LIST=>123>456.INS.1;TYPE=INS;REFLEN=1;SVLEN=6  1|1  1|1  1|0
```

**2. SV merging table**:

A TSV file mapping original SVs (`Variant_ID`) to merged SVs (`Collapse_ID`). For example:
   
| Chrom | Position | Bubble_ID          | Variant_ID                      | Collapse_ID                     | PctSeqSimilarity | PctSizeSimilarity | PctRecOverlap | SizeDiff | StartDistance | EndDistance | TruScore |   |   |
|-------|----------|--------------------|---------------------------------|---------------------------------|------------------|-------------------|---------------|----------|---------------|-------------|----------|---|---|
| chr22 | 16386947 | >38058649>38058909 | >38058649>38058909.DEL.33       | >38058649>38058909.DEL.34       | 0.997            | 1.000             | 0.990         | 0        | 3             | 3           | 99.5     |   |   |
| chr22 | 16386970 | >38058649>38058909 | >38058649>38058909.COMPLEX.41_2 | >38058649>38058909.COMPLEX.40_2 | 0.992            | 1.000             | 1.000         | 0        | 0             | 0           | 99.7     |   |   |
| chr22 | 16386973 | >38058649>38058909 | >38058649>38058909.INS.48       | >38058649>38058909.INS.51       | 0.950            | 0.954             | 0.948         | -7       | 0             | 0           | 95.0     |   |   |
| chr22 | 16387000 | >38058649>38058909 | >38058649>38058909.INS.56       | >38058649>38058909.INS.55       | 0.996            | 1.000             | 1.000         | 0        | 0             | 0           | 99.9     |   |   |
| chr22 | 16387057 | >38058649>38058909 | >38058649>38058909.INS.81       | >38058649>38058909.COMPLEX.78_2 | 1.000            | 1.000             | 1.000         | 0        | 0             | 0           | 100.0    |   |   |

In this table, the first row indicates the SV `>38058649>38058909.DEL.33` at `chr22:16386947` from bubble `>38058649>38058909` is merged into `>38058649>38058909.DEL.34`. The results of SV comparison performed by `Truvari` are included in columns from `PctSeqSimilarity` to `TruScore`.

Additionally, VCF INFO fields can be included as separate columns by specifying `--info`. If `--info SVLEN` is used, the output SVLEN in the tsv file will be REFLEN for COMPLEX and INV to indicate the value used for comparison. While in the output VCF, INFO/SVLEN is always calculated by `len(alt) - len(ref)`.

**3. Similar SV pairs with conflicting genotypes**

A TSV file listing SVs (`Variant_ID`) that are similar to another SV (`Collapse_ID`) but have conflicting genotypes.

**Arguments**:

```
usage: collapse_bubble.py [-h] -i VCF -o VCF -m PREFIX [--chr CHR] [--info TAG] [-l 50] [-r 100] [-p 0.9] [-P 0.9] [-O 0.9] [--debug]

Collapse biallelic SVs within the same bubble in VCF

Input / Output arguments:
  -i VCF, --invcf VCF   Input VCF
  -o VCF, --outvcf VCF  Output VCF
  -m PREFIX, --map PREFIX
                        Write collapsed and conflicting SV tables to PREFIX.collapse.txt and PREFIX.conflict.txt.
  --chr CHR             chromosome to work on. Default: all
  --info TAG            Comma-separated INFO/TAG list to include in the output map. Default: None

Collapse arguments:
  -l 50, --min-len 50   Minimum allele length of variants to be included, defined as max(len(alt), len(ref)). Default: 50
  -r 100, --refdist 100
                        Max reference location distance. Default: 100
  -p 0.9, --pctseq 0.9  Min percent sequence similarity (REF for DEL, ALT for other SVs). Default: 0.9
  -P 0.9, --pctsize 0.9
                        Min percent size similarity (SVLEN for INS, DEL; REFLEN for INV, COMPLEX). Default: 0.9
  -O 0.9, --pctovl 0.9  Min pct reciprocal overlap. Default: 0.9

Other arguments:
  --debug               Debug mode
```

## Todo list

Major:
- [ ] workflow script to run all analysis in the pipeline.
- [x] new SV merging pipeline that concatenates overlapping variants before SV merging, which will improve variant merging for tandem repeats and more compatible with the latest MC pipeline.

Minor:
- [x] update AC, AN, AF in the output VCF
- [x] support ploidy > 2


## Acknowledgement

- We thank Adam English for [Truvari](https://github.com/ACEnglish/truvari) and for sharing the insightful [script](https://github.com/ACEnglish/truvari/issues/228#issuecomment-2308535253) that inspired `collapse_bubble.py`.
- We thank Glenn Hickey for valuable suggestions on merging overlapping variants.
