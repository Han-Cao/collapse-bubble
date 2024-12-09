#!/usr/bin/env python3

# Merge duplicated variants in phased VCF
# Last update: 08-Dec-2024
# Author: Han Cao

import argparse
import logging
from collections import defaultdict
from dataclasses import dataclass

import pysam


logger = logging.getLogger(__name__)

@dataclass
class VariantCounter:
    """ Counter for different types of variants """
    input: int = 0              # Number of input variants
    output: int = 0             # Number of output variants
    merge: int = 0              # Number of varaints merged by genotype
    concat: int = 0             # Number of concatenated varaints
    conflict_missing: int = 0   # Number of varaints with missing genotype conflict


class VariantDict:
    """ (ref, alt) -> list of variants """

    def __init__(self):
        self.var_dict = defaultdict(list)


    def __getitem__(self, key):
        return self.var_dict[key]

    
    def __len__(self):
        return len(self.var_dict)

    
    def add_variant(self, variant: pysam.VariantRecord):
        """ Add variant to dictionary when reading from VCF """
        assert len(variant.alts) == 1
        self.var_dict[variant.alleles].append(variant)

    
    def get_duplicates(self) -> tuple[pysam.VariantRecord, pysam.VariantRecord]:
        """ Get a pair of duplicated variants """

        for _, var_lst in self.var_dict.items():
            if len(var_lst) > 1:
                return var_lst[0], var_lst[1]
        
        return None, None


    def update_merging(self, merged_var: pysam.VariantRecord, concat_var: pysam.VariantRecord):
        """ Update the dict after merging 2 duplicates """

        # always update the merged variant
        self.var_dict[merged_var.alleles][0] = merged_var
        self.var_dict[merged_var.alleles].pop(1)

        # if a concat variant is created, update
        if concat_var is not None:
            self.var_dict[concat_var.alleles].append(concat_var) 


    def get_sorted_keys(self):
        return sorted(self.var_dict.keys())



# enumerate possible haplotype merging results
HAP_MERGE_DICT = {
    (None, None): (None, None),
    (None, 0): (None, 0),
    (0, None): (None, 0),
    (None, 1): (1, None),
    (1, None): (1, None),
    (0, 0): (0, 0),
    (0, 1): (1, 0),
    (1, 0): (1, 0),                     
    (1, 1): (0, 1)
}

HAP_MERGE_DICT_MIS_REF = {
    (None, None): (None, None),
    (None, 0): (0, 0),
    (0, None): (0, 0),
    (None, 1): (1, 0),
    (1, None): (1, 0),
    (0, 0): (0, 0),
    (0, 1): (1, 0),
    (1, 0): (1, 0),                     
    (1, 1): (0, 1)
}


def merge_two_gt(gt1: tuple, gt2: tuple, mis_as_ref: bool) -> tuple:
    """ Merge two genotypes """

    gt_merge = []
    gt_concat = []
    ploidy = len(gt1)

    for i in range(ploidy):
        if mis_as_ref:
            gt_merge_i, gt_concat_i = HAP_MERGE_DICT_MIS_REF[(gt1[i], gt2[i])]
        else:
            gt_merge_i, gt_concat_i = HAP_MERGE_DICT[(gt1[i], gt2[i])]
        gt_merge.append(gt_merge_i)
        gt_concat.append(gt_concat_i)    
    
    return tuple(gt_merge), tuple(gt_concat)


def merge_two_variant(var1: pysam.VariantRecord, var2: pysam.VariantRecord, 
                      counter: VariantCounter, mis_as_ref: bool) -> list:
    """ 
    Merge two varaints
    Return:
        1. original variant with merged genotypes
        2. variant alleles will be concat if see 2 alleles on the same haplotype, otherwise None
    """

    merged_var = var1.copy()
    concat_var = var2.copy()
    flag_concat = False
    flag_missing_conflict = False

    for i in range(len(var1.samples)):
        gt1 = var1.samples[i]['GT']
        gt2 = var2.samples[i]['GT']

        if not flag_missing_conflict and check_missing_conflict(gt1, gt2):
            flag_missing_conflict = True

        gt_merge, gt_concat = merge_two_gt(gt1, gt2, mis_as_ref)
        merged_var.samples[i]['GT'] = gt_merge
        merged_var.samples[i].phased = True
        concat_var.samples[i]['GT'] = gt_concat
        concat_var.samples[i].phased = True

        if not flag_concat and 1 in gt_concat:
            flag_concat = True

    if flag_concat:
        concat_var.alleles = concat_alleles(concat_var.alleles)
        counter.concat += 1
    else:
        concat_var = None
        counter.merge += 1

    if flag_missing_conflict:
        # logger.warning(f'Sample {var1.samples[i].name} has missing genotype conflict at {var1.chrom}:{var1.pos} {var1.ref} {var1.alts[0]}')
        counter.conflict_missing += 1

    return merged_var, concat_var


def merge_var_dict(var_dict: VariantDict, counter: VariantCounter, mis_as_ref: bool, no_id: bool) -> VariantDict:
    """ Perform variant mergeing per position """

    while True:
        dup_var1, dup_var2 = var_dict.get_duplicates()
        if dup_var1 is None:
            break

        merged_var, concat_var = merge_two_variant(dup_var1, dup_var2, counter, mis_as_ref)
        if not no_id:
            merged_var = add_info_id(merged_var, dup_var2, 'DUP_ID')
            if concat_var is not None:
                concat_var = add_info_id(concat_var, dup_var1, 'CONCAT_ID')
        var_dict.update_merging(merged_var, concat_var)
    
    return var_dict


def concat_alleles(alleles: tuple[str, str]) -> tuple:
    """ Concatenate duplicated alleles """

    # left-trim same bases to find the unique sequence
    ref_uniq = alleles[0]
    alt_uniq = alleles[1]

    while len(ref_uniq) > 0 and len(alt_uniq) > 0 and ref_uniq[0] == alt_uniq[0]:
        ref_uniq = ref_uniq[1:]
        alt_uniq = alt_uniq[1:]
       
    return alleles[0] + ref_uniq, alleles[1] + alt_uniq


def check_missing_conflict(gt1: tuple, gt2: tuple) -> bool:
    """ Check if two haplotypes are co-missing """

    for hap1, hap2 in zip(gt1, gt2):
        if (hap1 is None) ^ (hap2 is None):
            return True
    
    return False


def add_info_id(var_keep: pysam.VariantRecord, var_add: pysam.VariantRecord, key: str) -> pysam.VariantRecord:
    """ Add DUP_ID or CONCAT_ID to INFO """

    id_lst = [var_keep.info[key][0]] if key in var_keep.info else []
    id_lst += [var_keep.info[key][0]] if key in var_add.info else []
    id_lst += [var_add.id]

    var_keep.info[key] = ':'.join(id_lst)

    return var_keep


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
    if an == 0:
        variant.info['AF'] = (0,)
    else:
        variant.info['AF'] = (ac / an,)
   
    return variant


def write_vcf(outvcf: pysam.VariantFile, var_dict: VariantDict, prev_var: pysam.VariantRecord, 
              counter: VariantCounter, mis_as_ref: bool, no_id: bool) -> None:
    """ Merge duplicated records and write to output VCF """

    # only 1 variant at the previous position
    if len(var_dict) == 0:
        outvcf.write(prev_var)
        counter.output += 1
        return None

    # more than 1 varaints at the previous position
    var_dict = merge_var_dict(var_dict, counter, mis_as_ref, no_id)

    for key in var_dict.get_sorted_keys():
        # assert len(var_dict[key]) == 1          # TEST ONLY, should always be true
        out_var = update_ac(var_dict[key][0])
        outvcf.write(out_var)
        counter.output += 1


def add_header(header: pysam.VariantHeader) -> pysam.VariantHeader:
    """ Add required headers """
    
    header.add_line('##INFO=<ID=DUP_ID,Number=A,Type=String,Description="Colon-separated duplicated variants merged into this variant">')
    header.add_line('##INFO=<ID=CONCAT_ID,Number=A,Type=String,Description="Colon-separated variants concatenated into this variant">')
    # header.add_line('##INFO=<ID=DROP_ID,Number=A,Type=String,Description="Colon-separated duplicated varinats dropped due to haplotype conflict">')

    return header


def check_vcf(vcf: pysam.VariantFile, n: int) -> None:
    """ Check VCF format """

    for i, var in enumerate(vcf.fetch()):
        if i > n:
            return

        for sample in var.samples.values():
            ploidy = len(sample['GT'])
            if ploidy > 2:
                raise ValueError(f'Ploidy {ploidy} is not supported yet')
            if ploidy == 2 and not sample.phased:
                raise ValueError(f'Unphased diploid sample is not supported')


def main(args: argparse.Namespace):

    # setup logger
    logging.basicConfig(level=logging.INFO,
                        format='[%(asctime)s] - [%(levelname)s]: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    
    # initialize
    invcf = pysam.VariantFile(args.invcf, 'rb')
    check_vcf(invcf, n=1000)
    if args.no_info_id:
        header = invcf.header
    else:
        header = add_header(invcf.header)
    outvcf = pysam.VariantFile(args.outvcf, 'w', header=header)
    counter = VariantCounter()
    # (ref, alt) -> list of variants
    # this store variants (if more than 1) at the position of the previous variant
    working_var_dict = VariantDict()


    logger.info(f'Merge duplicated variants from: {args.invcf}')
    # to check duplicates, we always write variants after reading its next one
    # if cur_var.pos = prev_var.pos, save them to working_var
    # if cur_var.pos > prev_var.pos, write all variants in prev_var.pos
    invcf_iter = invcf.fetch()
    prev_var = next(invcf_iter)
    counter.input += 1

    for cur_var in invcf_iter:
        counter.input += 1

        # process all variants at the previous position
        if cur_var.pos > prev_var.pos:
            write_vcf(outvcf, working_var_dict, prev_var, counter, args.merge_mis_as_ref, args.no_info_id)
            working_var_dict = VariantDict()
        # store variants at the same position
        elif cur_var.pos == prev_var.pos:
            # for the first 2 variants, store both
            if len(working_var_dict) == 0:
                working_var_dict.add_variant(prev_var)
            working_var_dict.add_variant(cur_var)
        else:
            # end of current chromosome, write all variants on the previous chromosome
            if prev_var.chrom != cur_var.chrom:
                write_vcf(outvcf, working_var_dict, prev_var, counter, args.merge_mis_as_ref, args.no_info_id)
                working_var_dict = VariantDict()
            # only unsorted VCF will have cur_pos < prev_pos, report error and exit
            else:
                raise ValueError('Input VCF is not sorted')

        prev_var = cur_var.copy()
                
    # process at the end of the VCF
    write_vcf(outvcf, working_var_dict, prev_var, counter, args.merge_mis_as_ref, args.no_info_id)

    logger.info(f'Read {counter.input} variants')
    logger.info(f'{counter.merge} variants are merged')
    logger.info(f'{counter.concat} variants are concatenated')
    if counter.conflict_missing > 0:
        logger.info(f'{counter.conflict_missing} variants have non-missing genotypes merged with missing genotypes')

    logger.info(f'Write {counter.output} variants to {args.outvcf}')

    invcf.close()
    outvcf.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='merge_duplicates.py', description='Merge duplicated variants in phased VCF')
    parser.add_argument('-i', '--invcf', metavar='VCF', help='Input VCF, sorted and phased', required=True)
    parser.add_argument('-o', '--outvcf', metavar='VCF', help='Output VCF', required=True)
    parser.add_argument('--merge-mis-as-ref', action='store_true', 
                        help='Convert missing to ref when merging missing genotypes with non-missing genotypes')
    parser.add_argument('--no-info-id', action='store_true', 
                        help='Do not add INFO/DUP_ID and INFO/CONCAT_ID to output VCF')

    args = parser.parse_args()

    main(args)
