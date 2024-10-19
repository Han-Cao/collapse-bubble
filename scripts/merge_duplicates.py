#!/usr/bin/env python3

# Merge duplicated variants in phased VCF
# Last update: 19-Oct-2024
# Author: Han Cao

import argparse
import logging
from collections import defaultdict
from dataclasses import dataclass

import pysam


@dataclass
class VariantCounter:
    """ Counter for different types of variants """
    input: int = 0              # Number of input variants
    output: int = 0             # Number of output variants
    merge: int = 0              # Number of varaints merged by genotype
    conflict_hap: int = 0       # Number of varaints with haplotype conflict
    conflict_missing: int = 0   # Number of varaints with missing genotype conflict


def merge_two_hap(hap1: int, hap2: int) -> int:
    """ Merge two phased genotypes """
    
    if hap1 is None:
        return hap2
    if hap2 is None:
        return hap1
    hap_sum = hap1 + hap2
    # 1 haplotype should not have 2 variants at the same position
    if hap_sum > 1:
        return -1
    return hap_sum


def check_missing_conflict(gt1: tuple, gt2: tuple) -> bool:
    """ Check if two haplotypes are co-missing """

    for hap1, hap2 in zip(gt1, gt2):
        if (hap1 is None) ^ (hap2 is None):
            return True
    
    return False

def merge_two_gt(gt1: tuple, gt2: tuple) -> tuple:
    """ Merge two genotypes """

    gt_merge = []
    ploidy = len(gt1)

    for i in range(ploidy):
        gt_merge.append(merge_two_hap(gt1[i], gt2[i]))

    return tuple(gt_merge)


def update_ac(variant: pysam.VariantRecord) -> pysam.VariantRecord:
    """ Update AC, AN, AF with new genotypes """

    ac = 0
    an = 0
    
    for sample in variant.samples.values():
        gt = sample['GT']
        for x in gt:
            if x is not None:
                an += 1
            if x == 1:
                ac += 1

    # for AC and AF, Number=A
    variant.info['AC'] = (ac,)
    variant.info['AN'] = an
    variant.info['AF'] = (ac / an,)
   
    return variant


def merge_genotypes(var_lst: list, counter: VariantCounter, keep: str, missing: str) -> list:
    """ Merge genotypes of list of varaints, return merged variant """

    res_lst = []
    remain_var_lst = var_lst.copy()

    while len(remain_var_lst) > 0:
        merge_var = remain_var_lst.pop(0)
        # the last one, no need to merge
        if len(remain_var_lst) == 0:
            res_lst.append(merge_var)
            break
        
        finish_var = []
        dup_id = []

        if keep == 'first':
            drop_id = []
        # pairwise genotype merging
        for other_var in remain_var_lst:
            flag_conflict_hap = False
            flag_conflict_missing = False
            new_gt_lst = []
            # merge genotype and save to list
            for i in range(len(merge_var.samples)):
                gt1 = merge_var.samples[i]['GT']
                gt2 = other_var.samples[i]['GT']
                new_gt = merge_two_gt(gt1, gt2)

                # check conflicts
                if not flag_conflict_missing:
                    flag_conflict_missing = check_missing_conflict(gt1, gt2)
                # -1 indicate haplotype conflict, don't merge
                if -1 in new_gt:
                    flag_conflict_hap = True
                    break
                new_gt_lst.append(new_gt)
            if flag_conflict_hap:
                counter.conflict_hap += 1
                if keep == 'first':
                    drop_id.append(other_var.id)
                continue

            if flag_conflict_missing:
                counter.conflict_missing += 1
                if missing == 'warn':
                    logger = logging.getLogger(__name__)
                    logger.warning(f'missing genotype conflict at {merge_var.chrom}:{merge_var.pos}: {merge_var.id} and {other_var.id}')
                elif missing == 'conflict':
                    counter.conflict_hap += 1
                    if keep == 'first':
                        drop_id.append(other_var.id)
                    continue
            
            # update merged genotypes
            for sample, new_gt in zip(merge_var.samples.values(), new_gt_lst):
                sample['GT'] = new_gt
                sample.phased = True
            
            finish_var.append(other_var)
            dup_id.append(other_var.id)
            counter.merge += 1
            
        # update list and results
        remain_var_lst = [x for x in remain_var_lst if x not in finish_var]
        merge_var = update_ac(merge_var)
        if len(dup_id) > 0:
            merge_var.info['DUP_ID'] = ":".join(dup_id)
        if keep == 'first' and len(drop_id) > 0:
            merge_var.info['DROP_ID'] = ":".join(drop_id)
        res_lst.append(merge_var)

        # if only output the first one, no need to loop
        if keep == 'first':
            break

    return res_lst
            

def write_vcf(outvcf: pysam.VariantFile, working_var: dict, prev_var: pysam.VariantRecord, 
              counter: VariantCounter, keep: str, missing: str) -> None:
    """ Merge duplicated records and write to output VCF """

    # only 1 variant at the previous position
    if len(working_var) == 0:
        outvcf.write(prev_var)
        counter.output += 1
        return

    # more than 1 varaints at the previous position
    out_lst = []
    for _, var_lst in working_var.items():
        if len(var_lst) > 1:
            # merge genotypes SVs
            merge_var_lst = merge_genotypes(var_lst, counter, keep, missing)
            out_lst += merge_var_lst
        else:
            out_lst.append(var_lst[0])
        
    for variant in out_lst:
        outvcf.write(variant)
    
    counter.output += len(out_lst)


def add_header(header: pysam.VariantHeader) -> pysam.VariantHeader:
    """ Add required headers """
    
    header.add_line('##INFO=<ID=DUP_ID,Number=A,Type=String,Description="Colon-separated duplicated variants merged into this variant">')
    header.add_line('##INFO=<ID=DROP_ID,Number=A,Type=String,Description="Colon-separated duplicated varinats dropped due to haplotype conflict">')

    return header

def main(args: argparse.Namespace):

    # setup logger
    logging.basicConfig(level=logging.INFO,
                        format='[%(asctime)s] - [%(levelname)s]: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger(__name__)
    
    # initialize
    invcf = pysam.VariantFile(args.invcf, 'rb')
    outvcf = pysam.VariantFile(args.outvcf, 'w', header=add_header(invcf.header))
    prev_var = None
    counter = VariantCounter()
    # (ref, alt) -> list of variants
    # this store variants (if more than 1) at the position of the previous variant
    working_var = defaultdict(list)


    logger.info(f'Merge duplicated variants from: {args.invcf}')
    # to check duplicates, we always write variants after reading its next one
    # if cur_var.pos = prev_var.pos, save them to working_var
    # if cur_var.pos > prev_var.pos, write all variants in prev_var.pos
    for cur_var in invcf.fetch():
        counter.input += 1
        # get the first variant
        if prev_var is None:
            prev_var = cur_var.copy()
            continue

        # process all variants at the previous position
        if cur_var.pos > prev_var.pos:
            write_vcf(outvcf, working_var, prev_var, counter, args.keep, args.missing)
            working_var = defaultdict(list)
        # store variants at the same position
        elif cur_var.pos == prev_var.pos:
            # for the first 2 variants, store both
            if len(working_var) == 0:
                working_var[(prev_var.ref, prev_var.alts[0])].append(prev_var)
            working_var[(cur_var.ref, cur_var.alts[0])].append(cur_var)
        else:
            # end of current chromosome, write all variants on the previous chromosome
            if prev_var.chrom != cur_var.chrom:
                write_vcf(outvcf, working_var, prev_var, counter, args.keep, args.missing)
                working_var = defaultdict(list)
            # only unsorted VCF will have cur_pos < prev_pos, report error and exit
            else:
                raise ValueError('Input VCF is not sorted')

        prev_var = cur_var.copy()
                
    # process at the end of the VCF
    write_vcf(outvcf, working_var, prev_var, counter, args.keep, args.missing)

    logger.info(f'Read {counter.input} variants')
    logger.info(f'{counter.merge} variants are merged')
    if args.missing == 'warn' and counter.conflict_missing > 0:
        logger.warning(f'{counter.conflict_missing} variants have non-missing genotypes merged with missing genotypes')

    if args.keep == 'all':
        logger.info(f'{counter.conflict_hap} variants are not merged due to haplotype conflict')
    elif args.keep == 'first':
        logger.info(f'{counter.conflict_hap} variants are dropped due to haplotype conflict')
    logger.info(f'Write {counter.output} variants to {args.outvcf}')

    invcf.close()
    outvcf.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='merge_duplicates.py', description='Merge duplicated variants in phased VCF')
    parser.add_argument('-i', '--invcf', metavar='VCF', help='Input VCF, sorted and phased', required=True)
    parser.add_argument('-o', '--outvcf', metavar='VCF', help='Output VCF', required=True)
    parser.add_argument('--keep', metavar='all|first', default='all', help='For duplicated variants with haplotype conflict  (e.g., 1|0 vs 1|0), output all of them (all) or only the first one (first).')
    parser.add_argument('--missing', metavar='warn|merge|conflict', default='warn',
                        help='For duplicated variants with missing conflict (e.g., .|0 vs 1|0), merge them with warning (warn), without warning (merge), or treat as haplotype conflict (conflict).')

    args = parser.parse_args()

    # check args
    if args.keep not in {'all', 'first'}:
        raise ValueError('--keep must be "all" or "first"')
    
    if args.missing not in {'warn', 'merge', 'conflict'}:
        raise ValueError('--missing must be "warn", "merge", or "conflict"')

    main(args)
