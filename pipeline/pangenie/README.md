# Use `collapse-bubble` to merge `PanGenie` results

This pipeline integrate `collapse-bubble` into the `PanGenie` genotyping pipeline. Please make sure you have read the [Pangenie pipeline](https://github.com/eblerjana/genotyping-pipelines/tree/main/prepare-vcf-MC) and installed required `PanGenie` scripts.

**Important**: This pipeline is under development and may change in the future. Some of the scripts reply on estimated Hardy-Weinberg Equilibrium (HWE) and allele frequencies (AF) in the population. It may not work well for small datasets. Please use this pipeline with caution.

## Prepare PanGenie reference

Input
- `Minigraph-Cactus` raw VCF: `mc.raw.vcf.gz`
- `Minigraph-Cactus` GFA: `mc.gfa`

Output
- Graph VCF with collapsed ID annotation: `mc.pangenie.collapse_id.vcf.gz`
- Collapsed biallelic VCF: `mc.pangenie.biallelic.uniqid.vcfwave.merge_dup.collapse.vcf.gz`
- ID annotation file: `mc.pangenie.collapse_id.txt.gz`
- PanGenie index: `mc.pangenie.index`

```bash
# vcfbub
vcfbub -l 0 -r 100000 -i mc.raw.vcf.gz > mc.vcfbub.r100k.vcf

# PanGenie's decomposition
python3 /path/to/pangenie/scripts/annotate_vcf.py \
-vcf mc.vcfbub.r100k.vcf \
-gfa mc.gfa \
-o mc.pangenie

# sort graph VCF
bcftools sort -o mc.pangenie.sort.vcf mc.pangenie.vcf
bcftools sort --write-index -Oz -o mc.pangenie.sort.vcf.gz mc.pangenie.sort.vcf
# sort and normalize biallelic VCF
bcftools norm -f ref.fa -Ou mc.pangenie.biallelic.vcf | bcftools sort --write-index -Oz -o mc.pangenie.biallelic.sort.vcf.gz

# annotate ID for biallelic VCF
python /path/to/collapse-bubble/scripts/annotate_var_id.py \
-i mc.pangenie.biallelic.sort.vcf.gz \
-o mc.pangenie.biallelic.uniqid.vcf.gz
tabix mc.pangenie.biallelic.uniqid.vcf.gz

# Pangenie's decomposition method cannot decompose all variants
# We use vcfwave to further decompose it
bcftools annotate -x INFO/AT -Ou mc.pangenie.biallelic.uniqid.vcf.gz | \
bcftools +fill-tags -- -t AC,AN,AF | \
vcfwave -I 1000 | bgzip -c > mc.pangenie.biallelic.uniqid.vcfwave.vcf.gz

bcftools norm -f ref.fa mc.pangenie.biallelic.uniqid.vcfwave.vcf.gz | \
bcftools sort --write-index -Oz -o mc.pangenie.biallelic.uniqid.vcfwave.sort.vcf.gz

# merge duplicates and overlap
python /path/to/collapse-bubble/scripts/merge_duplicates.py \
-i mc.pangenie.biallelic.uniqid.vcfwave.sort.vcf.gz \
-o mc.pangenie.biallelic.uniqid.vcfwave.sort.merge_dup.vcf.gz \
-c repeat \
--track ID \
tabix mc.pangenie.biallelic.uniqid.vcfwave.merge_dup.vcf.gz

# SV merging
python /path/to/collapse-bubble/scripts/collapse_bubble.py \
-i mc.pangenie.biallelic.uniqid.vcfwave.merge_dup.vcf.gz \
-o mc.pangenie.biallelic.uniqid.vcfwave.merge_dup.collapse.vcf.gz \
--map mc.pangenie.biallelic.uniqid.vcfwave.merge_dup.collapse.mapping

# sort
# PanGenie biallelic VCF: mc.pangenie.biallelic.uniqid.vcfwave.merge_dup.collapse.sort.vcf.gz
bcftools sort \
-m 4G -Oz -o mc.pangenie.biallelic.uniqid.vcfwave.merge_dup.collapse.sort.vcf.gz \
mc.pangenie.biallelic.uniqid.vcfwave.merge_dup.collapse.vcf.gz
tabix mc.pangenie.biallelic.uniqid.vcfwave.merge_dup.collapse.sort.vcf.gz

# annotate graph VCF with collapsed ID
# PanGenie graph VCF: mc.pangenie.collapse_id.vcf
python /path/to/collapse-bubble/pipeline/pangenie/annotate_graph_id.py \
--graph-vcf mc.pangenie.sort.vcf.gz \
--wave-vcf mc.pangenie.biallelic.uniqid.vcfwave.sort.vcf.gz \
--mergedup-vcf mc.pangenie.biallelic.uniqid.vcfwave.sort.merge_dup.vcf.gz \
--mapping mc.pangenie.biallelic.uniqid.vcfwave.merge_dup.collapse.mapping.collapse.txt \
-o mc.pangenie.collapse_id.vcf

# Prepare ID annotation
bcftools query -f '%CHROM\t%POS\t%INFO/ID\n' mc.pangenie.collapse_id.vcf | \
bgzip > mc.pangenie.collapse_id.txt.gz
tabix -s1 -b2 -e2 mc.pangenie.collapse_id.txt.gz

# prepare PanGenie index
PanGenie-index \
-v mc.pangenie.sort.vcf \
-r ref.fa \
-o mc.pangenie.index

PanGenie-index \
-v mc.pangenie.collapse_id.vcf \
-r ref.fa \
-o mc.pangenie.collapse_id.index
```

## Run PanGenie

Run PanGenie, and convert to biallelic VCF with or without SV merging:
```bash
# Genotyping
PanGenie \
-i <(zcat fastq.gz) \
-f mc.pangenie.index \
-o sample.pangenie.vcf \
-s sample

# convert to biallelic VCF without SV merging
cat sample.pangenie.vcf | \
python3 convert-to-biallelic.py mc.pangenie.biallelic.sort.vcf.gz | \
bgzip > sample.pangenie.biallelic.raw.vcf.gz

# convert to biallelic VCF with SV merging
# this need to annotate pangenie output with collapsed ID
bcftools annotate \
-a mc.pangenie.collapse_id.txt.gz \
-c CHROM,POS,INFO/ID \
--pair-logic "all" \
sample.pangenie.vcf | \
python3 convert-to-biallelic.py \
mc.pangenie.biallelic.uniqid.vcfwave.merge_dup.collapse.sort.vcf.gz | \
bgzip > sample.pangenie.biallelic.collapse.vcf.gz

# merge across samples
# this is necessary for merge_duplicates_pangenie.py
bcftools merge *raw.vcf.gz -Oz -o population.pangenie.biallelic.raw.vcf.gz
bcftools merge *collapse.vcf.gz -Oz -o population.pangenie.biallelic.collapse.vcf.gz
```

If you don't want the raw VCF without SV merging, you can directly run from the annotated graph VCF:

```bash
# Genotyping
PanGenie \
-i <(zcat fastq.gz) \
-f mc.pangenie.collapse_id.index \
-o sample.pangenie.vcf \
-s sample

# convert to biallelic VCF with SV merging
cat sample.pangenie.vcf | \
python3 convert-to-biallelic.py \
mc.pangenie.biallelic.uniqid.vcfwave.merge_dup.collapse.sort.vcf.gz | \
bgzip > sample.pangenie.biallelic.collapse.vcf.gz

# merge across samples
# this is necessary for merge_duplicates_pangenie.py
bcftools merge *raw.vcf.gz -Oz -o population.pangenie.biallelic.raw.vcf.gz
bcftools merge *collapse.vcf.gz -Oz -o population.pangenie.biallelic.collapse.vcf.gz
```

## Merge duplicates in biallelic PanGenie output

When converting PanGenie's results to biallelic, similar SVs (i.e., merged by `collapse_bubble.py`) within the same bubble will be merged into a single record:

```
Input:
T TAAA,TAAAA  0/0  0/1  1/2

Original PanGenie pipeline output:
T TAAA        0/0  0/1  1/0
T TAAAA       0/0  0/0  0/1

Collapse-bubble + PanGenie output:
T TAAA        0/0  0/1  1/1
```

However, if similar SVs reside in different bubbles, duplicated VCF records will be generated:
```
Input:
T TAAA         0/0  0/1  1/1
T TAAAA        0/0  0/1  0/0

Original PanGenie pipeline output:
T TAAA         0/0  0/1  1/1
T TAAAA        0/0  0/1  0/0

Collapse-bubble + PanGenie output:
T TAAA         0/0  0/1  1/1
T TAAA         0/0  0/1  0/0
```

The duplicated records can be merged by `merge_duplicates_pangenie.py`:

```
# IMPORTANT:
# this step relies on estimated AF and HWE, please only run it on population VCFs.
python /path/to/collapse-bubble/pipeline/pangenie/merge_duplicates_pangenie.py \
-i population.pangenie.biallelic.collapse.vcf.gz \
-o population.pangenie.biallelic.collapse.merge_dup.vcf.gz \
-r mc.pangenie.biallelic.uniqid.vcfwave.merge_dup.collapse.vcf.gz
```

This script can merge heterozygous genotypes into homozygous if possible, for example:

```
Input:
T TAAA         0/0  0/1  1/1
T TAAA         0/0  0/1  0/0

merge_duplicates_pangenie.py output:
T TAAA         0/0  1/1  1/1
```

Merging heterozygous genotypes is not suitable for all cases. If two bubbles share similar unique k-mers, PanGenie may genotype some SVs redundantly. Incorrectly merging two redundant heterozygous genotypes into one homozygous genotype would introduce errors. To determine whether merging is appropriate, we apply the following quality checks:

1. **Conflict check**: If the proportion of conflicting genotypes between two candidate SVs exceeds the threshold `--max-conflict`, the genotypes are not merged.
```
T TAAA    0/1  1/1  1/1
T TAAA    0/1  0/0  1/0
Conflict: No   No   Yes
```

2. **HWE check**: If genotype merging is not appropriate, it increases the discrepancy from HWE. We compute HWE p‑values for both the unmerged and merged scenarios. If `HWE_P(unmerge)` / `HWE_P(merge)` > `--hwe-ratio`, then merging is rejected. Since accurate genotypes should yield high HWE p‑values, a small `--hwe-ratio` is generally not required based on our test.

3. **Frequency check**: After merging, we compare the resulting AF against the reference panel frequencies (`mc.pangenie.biallelic.uniqid.vcfwave.merge_dup.collapse.vcf.gz`). If merging heterozygous genotypes increases the AF discrepancy from the reference panel, the merge is rejected. This ensures merging is performed only when it genuinely improves genotype accuracy.