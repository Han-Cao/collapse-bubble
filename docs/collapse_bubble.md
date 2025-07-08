# Collapse bubble

The aim of structural variantion (SV) merging is to remove redundant SVs by merging SVs with similar position, size, and sequences. `collapse_bubble.py` uses [Truvari](https://github.com/ACEnglish/truvari)'s engine to compare and merge SVs. Moreover, it has the following optimizations to better fit pangenome VCFs:

1. Only SVs within the same bubble or bubble clusters are compared for merging.
2. The phased haplotype is checked when merging SVs:
3. For SVs with both REF and ALT alleles > 1 bp, the SV size is determined by `REFLEN` (length of REF) instead of `SVLEN` (length difference between REF and ALT).

## Table of contents

- [SV merging within bubbles or bubble clusters](#sv-merging-within-bubbles-or-bubble-clusters)
- [Phased haplotype checking](#phased-haplotype-checking)
- [REFLEN vs SVLEN](#reflen-vs-svlen)

## SV merging within bubbles or bubble clusters

For a pangenome graph constructed from phased genome assemblies (e.g., `minigraph-cactus`), all variants are identified by genome alignments and represented in bubbles (or more generally [snarls](https://github.com/vgteam/vg/wiki/Snarls-and-chains)) in the graph. Since different bubbles do not overlap, most similar SVs should reside within the same bubble. However, it is possible that bubbles in repeat regions overlap with each other after variant normalization. Therefore, we further group those overlapping bubbles into bubble clusters. The SV merging procedure is as follows:

1. Identify all bubbles and their variants by checking `INFO/BUBBLE_ID`. For concatenated variants, the bubble IDs of original variants are retrieved from `INFO/CONCAT`
2. If any variants (include both small variants and SVs) from different bubbles are resided at the same position, the bubbles are grouped into a bubble cluster
3. Variants within the same bubble are pairwise compared and merged. The most common SV is selected as the representative SV.
4. After within bubble SV merging, SVs within the same bubble cluster are further compared and merged. The most common SV is selected as the representative SV.

## Phased haplotype checking

Since pangenome VCF is fully phased, the haplotype information can be leverage to avoid over-merging SVs. This approach is based on the fact that one haplotype cannot harbor more than 1 non-redundant SV at the same locus. When merging two SVs, the following rules are applied:

| SV1  | SV2  | Compatibility | Merged genotype                     |
| ---- | ---- | ------------- | ----------------------------------- |
| 0\|0 | 1\|0 | compatible    | 1\|0                                |
| 0\|1 | 1\|0 | compatible    | 1\|1                                |
| 1\|0 | 1\|0 | incompatible  | not merge if found in any haplotype |


## REFLEN vs SVLEN

`SVLEN` (i.e., length difference between REF and ALT alleles) is usually used to determine the size of SVs in most SV merging tools. This approach works well for simple insertions and deletions, which are the majority of SVs detected by read alignment-based callers. However, variants deconstructed from pangenome graph include many complex variants with both REF and ALT alleles > 1bp. `SVLEN` is not the best way to represent the size of such variants. For example:

```
COMPLEX1    ATCG      ACCTAAAAAAAAAA
COMPLEX2    ATCGTCG   ACCTAAAAAAAAAAAAA

```

In this case, the `SVLEN` of both variants are 10. However, COMPLEX2 is actually larger than COMPLEX1. Therefore, we use `REFLEN` (i.e., length of REF) instead to represent the size of such SVs.