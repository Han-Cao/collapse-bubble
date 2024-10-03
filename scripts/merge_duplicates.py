#!/usr/bin/env python3

# Merge duplicated variants in phased VCF
# Last update: 03-Oct-2024
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


def check_missing_conflict(hap1: int, hap2: int) -> bool:
    """ Check if two haplotypes are co-missing """
    
    return (hap1 is None) ^ (hap2 is None)


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


def merge_genotypes(var_lst: list, counter: VariantCounter) -> list:
    """ Merge genotypes of list of varaints, return merged variant or the first one (haplotype conflict) """

    res_lst = []
    remain_var_lst = var_lst.copy()

    while len(remain_var_lst) > 0:
        merge_var = remain_var_lst.pop(0)
        # the last one, no need to merge
        if len(remain_var_lst) == 0:
            res_lst.append(merge_var)
            break
        
        finish_var = []
        merge_id = []
        # TODO: for conflicting variants, should we only keep the first one and save dropped varaints to INFO/DROP_ID?
        # drop_id = []
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
                flag_conflict_missing = check_missing_conflict(gt1, gt2)
                # -1 indicate haplotype conflict, don't merge
                if -1 in new_gt:
                    flag_conflict_hap = True
                    break
                new_gt_lst.append(new_gt)
            if flag_conflict_hap:
                counter.conflict_hap += 1
                # drop_id.append(other_var.id)
                continue
            # TODO: we just report missing conflicts, should we also avoid merging?
            if flag_conflict_missing:
                logger = logging.getLogger(__name__)
                logger.warning(f'missing genotype conflict at {merge_var.chrom}:{merge_var.pos}: {merge_var.id} and {other_var.id}')
                counter.conflict_missing += 1
            
            # update merged genotypes
            for sample, new_gt in zip(merge_var.samples.values(), new_gt_lst):
                sample['GT'] = new_gt
                sample.phased = True
            
            finish_var.append(other_var)
            merge_id.append(other_var.id)
            counter.merge += 1
            
        # update list and results
        remain_var_lst = [x for x in remain_var_lst if x not in finish_var]
        merge_var = update_ac(merge_var)
        if len(merge_id) > 0:
            merge_var.info['MERGE_ID'] = ":".join(merge_id)
        # if len(drop_id) > 0:
        #     merge_var.info['DROP_ID'] = ":".join(drop_id)
        res_lst.append(merge_var)

    # update AC, AN, AF and return
    return res_lst
            

def write_vcf(outvcf: pysam.VariantFile, working_var: dict, prev_var: pysam.VariantRecord, counter: VariantCounter) -> None:
    """ Merge duplicated records and write to output VCF, return number of output variants """

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
            merge_var_lst = merge_genotypes(var_lst, counter)
            out_lst += merge_var_lst
        else:
            out_lst.append(var_lst[0])
        
    for variant in out_lst:
        outvcf.write(variant)
    
    counter.output += len(out_lst)


def add_header(header: pysam.VariantHeader) -> pysam.VariantHeader:
    """ Add required headers """
    
    header.add_line('##INFO=<ID=MERGE_ID,Number=A,Type=String,Description="Colon-separated duplicated variants merged into this variant">')
    # header.add_line('##INFO=<ID=DROP_ID,Number=A,Type=String,Description="Colon-separated duplicated varinats dropped due to haplotype conflict">')

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
    working_var = defaultdict(list) # (ref, alt) -> list of variants
    counter = VariantCounter()

    logger.info(f'Merge duplicated variants from: {args.invcf}')
    for cur_var in invcf.fetch():
        counter.input += 1
        # get the first variant
        if prev_var is None:
            prev_var = cur_var.copy()
            continue

        # process all variants at the previous position
        if cur_var.pos > prev_var.pos:
            write_vcf(outvcf, working_var, prev_var, counter)
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
                write_vcf(outvcf, working_var, prev_var, counter)
                working_var = defaultdict(list)
            # only unsorted VCF will have cur_pos < prev_pos, report error and exit
            else:
                raise ValueError('Input VCF is not sorted')

        prev_var = cur_var.copy()
                
    # process at the end of the VCF
    write_vcf(outvcf, working_var, prev_var, counter)

    logger.info(f'Read {counter.input} variants')
    logger.info(f'{counter.merge} variants are merged')
    if counter.conflict_missing > 0:
        logger.warning(f'{counter.conflict_missing} variants have non-missing genotypes merged with missing genotypes')
    logger.info(f'{counter.conflict_hap} variants are not merged due to haplotype conflict')
    logger.info(f'Write {counter.output} variants to {args.outvcf}')

    invcf.close()
    outvcf.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='graph_to_collapse.py', description='Map variants in PanGenie graph VCF to collapsed variants')
    parser.add_argument('-i', '--invcf', metavar='VCF', help='Input VCF, sorted and phased', required=True)
    parser.add_argument('-o', '--outvcf', metavar='VCF', help='Output VCF', required=True)

    args = parser.parse_args()
    main(args)
