# Use `collapse-bubble` to merge `PanGenie` results

This pipeline integrate `collapse-bubble` into the `PanGenie` genotyping pipeline. Please make sure you have read the [Pangenie pipeline](https://github.com/eblerjana/genotyping-pipelines/tree/main/prepare-vcf-MC) and installed required `PanGenie` scripts.

**Important**: This pipeline is under development and may change in the future. Use at your own risk.

## Prepare PanGenie reference

Input
- `Minigraph-Cactus` raw VCF: `mc.raw.vcf.gz`
- `Minigraph-Cactus` GFA: `mc.gfa`

Output
- Collapsed biallelic VCF: `mc.pangenie.biallelic.uniqid.vcfwave.merge_dup.collapse.vcf.gz`
- Graph VCF with collapsed ID annotation: `mc.pangenie.collapse_id.vcf.gz`

```bash
# vcfbub
vcfbub -l 0 -r 100000 -i mc.raw.vcf.gz > mc.vcfbub.r100k.vcf

# PanGenie's decomposition
python3 /path/to/pangenie/scripts/annotate_vcf.py \
-vcf mc.vcfbub.r100k.vcf \
-gfa mc.gfa \
-o mc.pangenie

# sort graph VCF
bcftools sort --write-index -Oz -o mc.pangenie.sort.vcf.gz mc.pangenie.vcf
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

bcftools sort \
-m 4G -Oz -o mc.pangenie.biallelic.uniqid.vcfwave.merge_dup.collapse.sort.vcf.gz \
mc.pangenie.biallelic.uniqid.vcfwave.merge_dup.collapse.vcf.gz
tabix mc.pangenie.biallelic.uniqid.vcfwave.merge_dup.collapse.sort.vcf.gz

# annotate graph VCF with collapsed ID
# this can be used as PanGenie reference Graph VCF
python /path/to/collapse-bubble/pipeline/pangenie/annotate_graph_id.py \
--graph-vcf mc.pangenie.sort.vcf.gz \
--wave-vcf mc.pangenie.biallelic.uniqid.vcfwave.sort.vcf.gz \
--mergedup-vcf mc.pangenie.biallelic.uniqid.vcfwave.sort.merge_dup.vcf.gz \
--mapping mc.pangenie.biallelic.uniqid.vcfwave.merge_dup.collapse.mapping.collapse.txt \
-o mc.pangenie.collapse_id.vcf.gz
tabix mc.pangenie.collapse_id.vcf.gz
```

## Merge duplicates in biallelic PanGenie output

```
python /path/to/collapse-bubble/pipeline/pangenie/merge_duplicates_pangenie.py \
-i pangenie.biallelic.vcf.gz \
-o pangenie.biallelic.merge_dup.vcf.gz \
-r mc.pangenie.biallelic.uniqid.vcfwave.merge_dup.collapse.vcf.gz
```