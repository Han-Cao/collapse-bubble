#!/usr/bin/env python3

# Merge duplicated variants in biallelic PanGenie results
# Last update: 23-Jul-2025
# Author: Han Cao

import argparse
import sys
from collections import defaultdict

import pysam
import numpy as np

# AC <-> GT conversion
# biallelic PanGenie VCF should only have these genotypes
AC_DICT = {
    (None, None): -1,
    (0, 0): 0,
    (0, 1): 1,
    (1, 0): 1,
    (1, 1): 2
}

GT_DICT = {
    -1: (None, None),
    0: (0, 0),
    1: (0, 1),
    2: (1, 1)
}


def read_sample_list(fname: str) -> list:
    """ Read sample list """

    with open(fname, 'r') as f:
        samples = [x.strip() for x in f.readlines()]
    
    return samples


def ac_array(variant: pysam.VariantRecord) -> np.ndarray:
    """ Convert genotypes to AC array """

    ac_lst = [AC_DICT[sample['GT']] for sample in variant.samples.values()]

    return np.array(ac_lst)


def get_af(variant: pysam.VariantRecord, sample_list: list=None) -> float:
    """ Get allele frequency from list of samples """

    if sample_list is None:
        sample_list = variant.samples.keys()

    ac = 0
    an = 0

    for sample in sample_list:
        for gt in variant.samples[sample]['GT']:
            if gt is not None:
                an += 1
            if gt == 1:
                ac += 1

    af = ac / an if an != 0 else 0

    return af


def count_gt(ac_arr: np.ndarray) -> tuple[int, int, int]:
    """ Count genotypes """

    ac_hom1 = np.sum(ac_arr == 0)
    ac_hets = np.sum(ac_arr == 1)
    ac_hom2 = np.sum(ac_arr == 2)

    return ac_hom1, ac_hets, ac_hom2


# equivalent to https://csg.sph.umich.edu/abecasis/Exact/snp_hwe.r
def exact_hwe(obs_hom1: int, obs_hets: int, obs_hom2: int) -> float:

    # Prepare genotype counts
    N = obs_hom1 + obs_hom2 + obs_hets
    obs_homr = min(obs_hom1, obs_hom2)
    # obs_homc = max(obs_hom1, obs_hom2)
    rare = obs_homr * 2 + obs_hets

    # Initialize probability array
    probs = np.zeros(rare + 1)

    # Find midpoint of the distribution
    mid = int(rare * (2 * N - rare) / (2 * N))
    if (mid % 2) != (rare % 2):
        mid += 1

    probs[mid] = 1.0
    mysum = 1.0

    # Calculate probabilities from midpoint down
    curr_hets = mid
    curr_homr = (rare - mid) / 2
    curr_homc = N - curr_hets - curr_homr

    while curr_hets >= 2:
        probs[curr_hets - 2] = probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0))
        mysum += probs[curr_hets - 2]

        curr_hets -= 2
        curr_homr += 1
        curr_homc += 1

    # Calculate probabilities from midpoint up
    curr_hets = mid
    curr_homr = (rare - mid) / 2
    curr_homc = N - curr_hets - curr_homr
    
    while curr_hets <= rare - 2:
        probs[curr_hets + 2] = probs[curr_hets] * 4.0 * curr_homr * curr_homc / ((curr_hets + 2.0) * (curr_hets + 1.0))
        mysum += probs[curr_hets + 2]
        
        curr_hets += 2
        curr_homr -= 1
        curr_homc -= 1

    # P-value calculation
    target = probs[obs_hets]
    p = min(1.0, np.sum(probs[probs <= target]) / mysum)
    
    return p


def conflict_ratio(ac_arr1: np.ndarray, ac_arr2: np.ndarray) -> bool:
    """ Calculate the ratio of conflicting genotypes (1/1 vs 1/-) """

    ac_arr_sum = ac_arr1 + ac_arr2
    n_hom_alt = np.sum((ac_arr1 == 2) | (ac_arr2 == 2))
    if n_hom_alt == 0:
        return 0
    else:
        return np.sum(ac_arr_sum > 2) / n_hom_alt


def merge_ac(ac_arr1: np.ndarray, ac_arr2: np.ndarray, hwe_ratio: float) -> np.ndarray:
    """ Merge AC arrays by sum or max """

    # convert -1 to 0 to allow merging with missing
    ac_arr1_fix = np.where(ac_arr1 == -1, 0, ac_arr1)
    ac_arr2_fix = np.where(ac_arr2 == -1, 0, ac_arr2)
    
    # sum and max of two AC arrays
    ac_arr_sum = ac_arr1_fix + ac_arr2_fix
    ac_arr_sum[ac_arr_sum > 2] = 2
    ac_arr_max = np.maximum(ac_arr1, ac_arr2)

    # set missing to -1
    mis_mask = (ac_arr1 == -1) & (ac_arr2 == -1)
    ac_arr_sum[mis_mask] = -1
    ac_arr_max[mis_mask] = -1

    # early return if two method are the same
    if np.array_equal(ac_arr_sum, ac_arr_max):
        return ac_arr_sum

    # compare HWE between sum and max using non-missing AC
    sum_hom1, sum_hets, sum_hom2 = count_gt(ac_arr_sum)
    max_hom1, max_hets, max_hom2 = count_gt(ac_arr_max)

    # add small number to avoid 0/0
    sum_p = exact_hwe(sum_hom1, sum_hets, sum_hom2) + 1e-10
    max_p = exact_hwe(max_hom1, max_hets, max_hom2) + 1e-10

    # if two genotypes pass conflict cutoff, we prefer sum
    # so we do max only when max_p / sum_p > hwe_ratio
    if (max_p / sum_p) > hwe_ratio:
        return ac_arr_max
    else:
        return ac_arr_sum
    

def merge_variant(variant: pysam.VariantRecord, ac_arr: np.ndarray) -> pysam.VariantRecord:
    """ Update genotypes from AC array """
    
    # Note: By default, PanGenie biallelic VCF does not have AC related INFO
    # so we don't need to update those INFO fields, manually add them by bcftools later
    
    new_var = variant.copy()
    for i in range(len(new_var.samples)):
        new_var.samples[i]['GT'] = GT_DICT[ac_arr[i]]
    
    return new_var


def merge_duplicates(var_lst: list[pysam.VariantRecord], hwe_ratio: float, max_conflict: float) -> pysam.VariantRecord:
    """ 
    Merge duplicated variants as much as possible, can return more than 1 if conflict

    To speed up, we first check compatibility between the most common variant and the rest
    This procedure will loop over all remaing variants until no more variants can be merged
    """

    ret_lst = []
    ac_arr_lst = []

    # vectorize genotypes
    for var in var_lst:
        ac_arr = ac_array(var)
        ac_arr_lst.append(ac_arr)
    
    # rank variant index by AC
    ac_lst = [np.sum(x[x != -1]) for x in ac_arr_lst]
    remain_var_idx = list(np.argsort(ac_lst)[::-1])

    # merge from the most common one
    while len(remain_var_idx) > 0:
        merge_ac_arr = ac_arr_lst[remain_var_idx[0]]
        drop_lst = [0]

        for i in range(1, len(remain_var_idx)):
            other_ac_arr = ac_arr_lst[remain_var_idx[i]]

            if conflict_ratio(merge_ac_arr, other_ac_arr) > max_conflict:
                continue
            
            merge_ac_arr = merge_ac(merge_ac_arr, other_ac_arr, hwe_ratio)
            drop_lst.append(i)

        # update variant after merging all compatible variants
        if len(drop_lst) > 1:
            merge_var = merge_variant(var_lst[remain_var_idx[0]], merge_ac_arr)
            ret_lst.append(merge_var)
        # if unable to merge any variants, no need to update genotypes
        else:
            ret_lst.append(var_lst[remain_var_idx[0]])
        
        # remove merged variants
        drop_set = set(drop_lst)
        remain_var_idx = [x for i, x in enumerate(remain_var_idx) if i not in drop_set]

    return ret_lst
    
 
def pos_workder(var_lst: list[pysam.VariantRecord], 
                outvcf: pysam.VariantFile,
                refvcf: pysam.VariantFile,
                in_samples: list, 
                ref_samples: list,
                hwe_ratio: float,
                max_conflict: float) -> int:
    """ Worker to dedup and write variants at the same position """

    dup_dict = defaultdict(list) # id -> list of VariantRecord
    n_remove = 0

    # group variants by id
    for var in var_lst:
        dup_dict[var.id].append(var)
    
    # dedup and write
    for var_dups in dup_dict.values():
        if len(var_dups) > 1:
            # merge duplicates as much as possible
            merge_var_lst = merge_duplicates(var_dups, hwe_ratio, max_conflict)
            # if any conflict, we choose the one with AF close to reference
            if len(merge_var_lst) > 1:
                merge_var_af_lst = [get_af(x, in_samples) for x in merge_var_lst]
                ref_var = fetch_variant(refvcf, var.chrom, var.pos, var.id)
                ref_var_af = get_af(ref_var, ref_samples)
                af_diff = np.array(merge_var_af_lst) - ref_var_af
                merge_var_idx = np.argmin(np.abs(af_diff))
            else:
                merge_var_idx = 0
            
            var_out = merge_var_lst[merge_var_idx]
            n_remove += len(var_dups) - 1
        else:
            var_out = var_dups[0]

        outvcf.write(var_out)

    return n_remove


def fetch_variant(vcf: pysam.VariantFile, chr:str, pos: int, id: str) -> pysam.VariantRecord:
    """ Fetch variant from VCF by id """

    for var in vcf.fetch(chr, pos - 1, pos):
        if var.id == id:
            return var
    
    raise ValueError(f'Variant {chr}:{pos} not found')


def main(args: argparse.Namespace) -> None:

    # read VCF
    invcf = pysam.VariantFile(args.input_vcf, 'r')
    refvcf = pysam.VariantFile(args.ref_vcf, 'r')
    outvcf = pysam.VariantFile(args.out_vcf, 'w', header=invcf.header)

    # parse sample lists
    if args.in_samples is not None:
        in_samples = read_sample_list(args.in_samples)
    else:
        in_samples = list(invcf.header.samples)
    
    if args.ref_samples is not None:
        ref_samples = read_sample_list(args.ref_samples)
    else:
        ref_samples = list(refvcf.header.samples)
    
    # main worker
    n_remove = 0
    invcf_iter = invcf.fetch()
    prev_var = next(invcf_iter)
    pos_var_lst = [prev_var]
    
    for cur_var in invcf_iter:
        # same chromosome
        if cur_var.chrom == prev_var.chrom:
            if cur_var.pos == prev_var.pos:
                pos_var_lst.append(cur_var)
            else:
                assert cur_var.pos > prev_var.pos, "Input VCF is not sorted"
                n_remove += pos_workder(pos_var_lst, outvcf, refvcf, in_samples, ref_samples, 
                                        args.hwe_ratio, args.max_conflict)
                pos_var_lst = [cur_var]

        # different chromosome
        else:
            n_remove += pos_workder(pos_var_lst, outvcf, refvcf, in_samples, ref_samples, 
                                    args.hwe_ratio, args.max_conflict)
            pos_var_lst = [cur_var]

        prev_var = cur_var.copy()

    # process the last position
    n_remove += pos_workder(pos_var_lst, outvcf, refvcf, in_samples, ref_samples, args.hwe_ratio, args.max_conflict)

    print(f'Removed {n_remove} duplicated variants', file=sys.stderr)

    invcf.close()
    refvcf.close()
    outvcf.close() 


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='dedub_biallelic.py', description='Merge duplicated variants in biallelic PanGenie results')
    parser.add_argument('-i', '--input-vcf', metavar='VCF', required=True, 
                        help='Biallelic Pangenie genotyping result VCF')
    parser.add_argument('-r', '--ref-vcf', metavar='VCF', required=True, 
                        help='Reference biallelic VCF')
    parser.add_argument('-o', '--out-vcf', metavar='VCF', required=True, 
                        help='Output VCF')
    parser.add_argument('--in-samples', metavar='TXT', default=None, 
                        help='Input sample list for AF calculation')
    parser.add_argument('--ref-samples', metavar='TXT', default=None, 
                        help='Reference sample list for AF calculation')
    parser.add_argument('--max-conflict', metavar='0.05', type=float, default=0.05,
                        help='Maximum ratio of conflicting genotypes (i.e., 1/1 vs 1/-) to allow merging two variants. Default: 0.05')
    parser.add_argument('--hwe-ratio', metavar='10', type=float, default=10, 
                        help='HWE P value ratio to determine whether merge two 0/1 genotypes into 0/1 or 1/1.' + \
                             'If HWE_P(0/1) / HWE_P(1/1) > hwe_ratio, it keeps as 0/1. Default: 10')
    
    args = parser.parse_args()
    main(args)