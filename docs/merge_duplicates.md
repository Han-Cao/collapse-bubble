# Merge duplicates

`merge_duplicates.py` performs a two-step process to merge duplicated VCF records and overlapping alleles at the same position:

**Step 1**: Concatenate the alleles at the same position if they are on the same haplotype. This step is similar as the haplotype reconstruction that has been disccused in [vt](https://github.com/atks/vt/issues/16).

```
C AA   0|1
C AA   0|1
C AAAA 1|0

to:   
C AA   0|0
C AAAA 1|0
C AAAA 0|1
```

**Step 2**: Merge the genotypes of duplicated variants: 
```
C AA   0|0
C AAAA 1|0
C AAAA 0|1

to:  
C AA   0|0
C AAAA 1|1
```
Here, we get a non-overlapping, deduplicated VCF.

## Table of contents

- [Input requirements](#input-requirements)
- [How it works](#how-it-works)
    - [Proof](#proof)
    - [Demonstration](#demonstration)
- [Notes](#notes)
- [Limitations](#limitations)

## Input requirements

When concatenating phased alleles, it requires that:

1. The input VCF should not have any overlapping alleles on the same haplotype before `bcftools norm`
2. The input VCF should be sorted by only chr and position after `bcftools norm`.

These ensure when 2 alleles are found at the same position, all variants except the first one must be left aligned to this position. Therefore, we can concatenate indels to the first variant to reconstruct the original haplotype. This makes it possbible to concatenate alleles without querying reference genome.

### Important note:
Sorting by `chr` and  `pos` only can be done using `bcftools v1.7` or earlier (see bcftools issue #756). For newer `bcftools`, variants are sorted by `chr`, `pos`, `ref`, `alt`. The following bash script can sort  by `chr` and `pos` only:

```
(bcftools view -h input.vcf.gz ; bcftools view -H input.vcf.gz | sort -s -k1,1d -k2,2n) | bgzip > output.vcf.gz
```

## How it works

**TL;DR**: To concatenate variant A and variant B, we can right shift B to the end of A's reference allele, and then concatenate their alleles, the concatenated allele is:

$$ Allele_{concat} =  Allele_A + S_B[n:m] + S_B[0:n] $$

where $Allele_A$ is the REF or ALT allele of variant A to be concatenated, $S_B$ is the left-trimed indel sequence of variant B, $m = len(S_B)$, and $n = [len(REF_A) - 1] \mod m$.

Particularly, for deletions, the above is equivalent to insert the left-trimed $REF_B$ after the first base of $REF_A$ :

$$ REF_{concat} = REF_A[0] + S_B + REF_A[1:len(REF_A)] $$

Example:
```
ref  GGCTAGCTA      (span from 1-9)
hap1 AAA-A----      (GGCT to AAA + GCTA del)
hap2 AAA-AGCTAGCTA  (GGCT to AAA + GCTA ins) 

VCF:
                           hap1   hap2
var   1   GGCT    AAA      1      1
del   1   GGCTA   G        1      0
ins   1   G       GGCTA    0      1

For hap1:
ref = var_ref[0] + S_del + var_ref[1:4]
    = G + GCTA + GCT
    = GGCTAGCT
alt = var_alt 
    = AAA

For hap2:
ref = var_ref
    = GGCT
alt = var_alt + S_ins[3:4] + S_ins[0:3]
    = AAA + A + GCT
    = AAAAGCT
Right trim to normalize:
ref = G
alt = AAAA
```


### Proof:

Definition:
- Allele: We assume all alleles are already normalized, which have at least 1 base in both REF and ALT, and have been left/right trimmed if possible
- Indel sequence: The sequence of an indel only include the inserted / deleted sequence, for example, given variant with alleles of GCCC G, the indel sequence is CCC

According to the [left-align algorithm](https://genome.sph.umich.edu/wiki/Variant_Normalization), we prove:

1. **Only indels can be left aligned**: 

   Variants are left aligned only if right trim generate empty ref or alt, so they are indels

2. **If an indel of m length is left aligned from A to B, we can right shift it to any position in [A,B]**:

   This can be done by revert left align

3. **If right shift an indel x bases, the new indel sequence is: s[n:m] + s[0:n], where n = x % m, s is the original indel sequence**

   Right shift 1 base: left trim 1 allele, extend the next base of reference genome to the right end of ref and alt. As this variant has been left aligned, the extended allele is always the same as the first base of the indel sequence Therefore, right shift an allele is the same as right shift a string if the reference genome allow right shift. For example:
```
ref GGCTAGCTA
alt GGCTA

             GGCTAGCTA
left aligned G----GCTA  (GGCTA G in VCF, del sequence is GCTA)
right shift  GG----CTA  (GCTAG G in VCF, del sequence is  CTAG)

-> right shift 1 base   =  s -> s[1:m] + s[0:1]
-> right shift n bases  =  s -> s[n:m] + s[0:n], if n <= m
-> right shift m bases  =  s -> s
-> right shift x bases  =  s -> s[n:m] + s[0:n], where n = x % m
```

4. **If an indel of length m can right shift longer than its length, then the indel is a tandem repeat, with repeat motif of length M, and m % M == 0**

   Let the indel right shift m bases, the new sequence is the same as the original one, then it is a tandem repeat

Moreover, in a non-overlapping VCF (requirement 1):

5. **For variant with reference allele spanning from position A to B, only indels can be left aligned to A, and we can right shift the indels back to B + 1**

   If not, there will be overlapping alleles on the same haplotype before left alignment.

Therefore, to concatenate variant A with an indel left aligned to the same position, we can always right shift the indel to the end of A's reference allele and concatenate the indel sequences to reconstruct the haplotype.


### Demonstration:
Given reference gneome G(GCTA)n, there is an complex variant at the left end: `GGCT` to `AAA` replacement. Moreover, there are also GCTA indels at the right end of GCTA repeats. Because the left end replacement and right end indels are far from each other, they were called as separated variants in VCF. Therefore, after variant normalization, the output VCF include 3 variants at the same position:

```
VCF representation:
                           hap1   hap2
var   1   GGCT    AAA      1      1
del   1   GGCTA   G        1      0
ins   1   G       GGCTA    0      1

For visualization, we only consider the sequence of first 2 repeats G(GCTA)2
ref  GGCTAGCTA      (span from 1-9)
hap1 AAA-A----      (GGCT to AAA + GCTA del)
hap2 AAA-AGCTAGCTA  (GGCT to AAA + GCTA ins) 
```

**Concatenate deletions**:

```
ref  GGCTAGCTA      (span from 1-9)
hap1 AAA-A----      (GGCT to AAA + GCTA del)

VCF representation:
var   1   GGCT     AAA
del   1   GGCTA    G
hap1  1   GGCTAGCT AAA   (expected output)

Concatenate procedure:
1. GGCTA G -> TAGCT T  (Right shift from pos 2 to pos 5)
2. GGCTAGCT  AAA       (Concat alleles, same as expected)
3. GGCTTCCTA AAAA      (Add the remaining alleles in reference)
-> AAAA is the reconstructed hap1

In general:
Let x = len(var_ref) - 1               (base to right shift)
    m = len(del),
    n = x % m,
    s = left-trimmed deletion sequence (i.e., GCTA)

Then, concatednated alleles are:
ref = var_ref + ( s[n:m] + s[0:n] )
    = GGCT + A + GCT 
    = GGCTAGCT

alt = var_alt 
    = AAA

Morever, the concatenated ref allele is also equivalent to inserting the indel sequence s after position 1:
ref = G + GCTA + GCT
    = GGCTAGCT

Proof:
if x < m, then 
n = x % m = x, and var_ref[1:] = s[0:n]
ref = var_ref + ( s[n:m] + s[0:n])
    = var_ref[0] + s[0:n] + s[n:m] + s[0:n]
    = var_ref[0] + s + var_ref[1:]
if x > m, then s is tandem repeat and x has at least 1 copy of full repeat motif. As the copy number of tandem repeat in ref will not affect the results, this is equivalent to x < m.
```

**Concatenate insertions**:

```
ref  GGCTAGCTA      (span from 1-9)
hap2 AAA-AGCTAGCTA  (GGCT to AAA + GCTA ins) 

VCF representation:
var   1   GGCT  AAA    
ins   1   G     GGCTA
hap2  1   G     AAAA   (expected output)

Concatenate procedure:
1. G GGCTA -> T TAGCT       (Right shift from pos 2 to pos 5)
2. GGCT  AAAAGCT            (Concat alleles)
3. GGCTAGCTA AAAAGCTAGCTA   (Add the remaining alleles in reference)
-> AAAAGCTAGCTA is the reconstructed hap1

Usually we need to further normalize the concatenated alleles by right trimming:
GGCT    AAAAGCT     (right trim GCT)
G       AAAA        (same as expected output)

In general:
Let x = len(var_ref) - 1               (base to right shift)
    m = len(ins),
    n = x % m,
    s = left-trimmed ins sequence      (i.e., GCTA)

Then, concatednated alleles are:
ref = var_ref 
    = GGCT

alt = var_alt + ( s[n:m] + s[0:n] )
    = AAA + A + GCT 
    = AAAAGCT

Right trim:
ref = G
alt = AAAA
```

Q.E.D

## Notes

After allele concatenation, the output VCF it possible to be further normalized by `bcftools norm`. For example:

```
Input
chr20   195745  >40785190>40785193      C       T
chr20   195745  >40785193>40785203_4    CGTGT   C

Output
chr20   195745  chr20:195745_0          CGTGT   T
```

In the new variant, both `REF` and `ALT` end with `T`, so it can be further left align to `chr20 195744 TCGTG   T`. This convert the original SNP + repeat deletion at chr20:195745 to a single deletion at chr20:195744, which could be misleading for downstream analysis. Therefore, it is not recommended to use `bcftools norm` after `merge_duplicates.py`.

## Limitations:

To concatenate alleles, the script only consider the combination of alternative alleles seen on existing haplotypes. This can correctly identify alternative alleles but may not be perfect to distinguish between the reference allele and missing genotypes. To perfectly determine whether one haplotypes has reference allele or missing genotype, we need a greedy search for all possible combinations:

```
              sample1   sample2
1  C AA       1|0       0|0
2  C AAA      0|0       .|0
3  C AAA      0|0       1|0
4  C AAAA     1|0       0|0

-> C AAAAAAA  1|0       ?|0
```

- For sample1, we concat 1,4 and get the final 6A insertion. 
- For sample2, if also consider 1,4, its genotype is `0|0`, but if we consider 2,3, it should be `.|0`.
 
However, a greedy search is too expensive when there are many variants to be considered. Since most combinations should not affect the final result, we currently only check the conbinations that are seen (i.e., 1,4 not 3,4) in existing haplotypes. 

If you already specify `--merge-mis-as-ref`, this would not be a problem. Because the genotypes of concatenated alleles are missing if and only if all existing genotypes are missing:

```
              sample1   sample2
1  C AA       1|0       0|.
2  C AAA      0|.       .|.
3  C AAA      0|.       1|.
4  C AAAA     1|.       0|.

-> C AAAAAAA  1|0       0|.
```
