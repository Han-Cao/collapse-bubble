#!/usr/bin/env python3

# Merge duplicated variants in phased VCF
# Last update: 14-Dec-2024
# Author: Han Cao

import argparse
import logging
from collections import defaultdict
from dataclasses import dataclass

import pysam
import numpy as np


logger = logging.getLogger(__name__)

GT_DICT = {0: 0,
           1: 1,
           -1: None}

@dataclass
class VariantCounter:
    """ Counter for different types of variants """
    input: int = 0              # Number of input variants
    output: int = 0             # Number of output variants
    merge: int = 0              # Number of varaints merged by genotype
    concat: int = 0             # Number of concatenated varaints
    # conflict_missing: int = 0   # Number of varaints with missing genotype conflict


class HapIndexer:
    """ Map haplotype index to sample index """

    def __init__(self, var: pysam.VariantRecord):

        self.idx_dict = {}      # VCF haplotype index -> sample index, ploidy, sample nth hap
        hap_idx = 0
        for samp_idx, sample in enumerate(var.samples.values()):
            ploidy = len(sample['GT'])
            for n_hap in range(ploidy):
                self.idx_dict[hap_idx] = (samp_idx, ploidy, n_hap)
                hap_idx += 1
    
    def __getitem__(self, key):
        return self.idx_dict[key]


class DuplicatesDict:
    """ (ref, alt) -> list of variants """

    def __init__(self) -> None:
        self.var_lst = []                   # no. -> pysam.VariantRecord
        self.n_var = 0                      # number of variants in var_lst
        self.n_merge = 0                    # number of variants merged
        self.dup_dict = defaultdict(list)   # alleles -> list -> var idx
        self.gt_mat = np.array([])          # haplotype genotype matrix


    def __getitem__(self, key) -> pysam.VariantRecord:
        return self.var_lst[key]

    
    def __len__(self) -> int:
        return self.n_var

    
    def add_variant(self, variant: pysam.VariantRecord) -> None:
        """ Add variant to dictionary when reading from VCF """

        if variant.ref is None or variant.alts is None:
            raise ValueError(f"REF and ALT alleles must be provided at {variant.chrom}:{variant.pos}")

        if len(variant.alts) != 1:
            raise ValueError(f"Only biallelic variants are supported, found {len(variant.alts)} alts at {variant.chrom}:{variant.pos}")
        

        self.var_lst.append(variant)
        self.dup_dict[variant.alleles].append(self.n_var)
        self.n_var += 1


    def concatenate(self, hap_indexer: HapIndexer, method: str, 
                    mis_as_ref: bool, track: str, max_repeat: int) -> int:
        """ Concatenate alleles """

        gt_mat = get_gt_mat(self.var_lst)

        if method == "none":
            self.gt_mat = gt_mat
            return 0
        
        new_vars = []
        if method == "position":
            # here, we pass the original gt_mat, so it is updated in-place
            new_vars, new_vars_gt_mat, _ = concat_variants(self.var_lst, 
                                                           gt_mat, 
                                                           hap_indexer, 
                                                           mis_as_ref, 
                                                           track)
        elif method == "repeat":
            motif_dict = defaultdict(list)
            new_vars_gt_lst = []
            # find repeat motif for all variants
            for i, var in enumerate(self.var_lst):
                ref, alt = var.alleles
                if not is_indel(ref, alt):
                    continue
                indel_seq = get_indel_seq(ref, alt)
                if max_repeat is not None and len(indel_seq) > max_repeat:
                    continue
                motif = find_rep_motif(indel_seq)
                motif_dict[motif].append(i)
            
            for _, concat_var_idx in motif_dict.items():
                if len(concat_var_idx) < 2:
                    continue
                # here, we pass a copy of gt_mat, so we need to update the gt_mat
                tmp_vars, tmp_gt_mat, update_gt_mat = concat_variants([self.var_lst[i] for i in concat_var_idx],
                                                                      gt_mat[concat_var_idx],
                                                                      hap_indexer,
                                                                      mis_as_ref,
                                                                      track)
                gt_mat[concat_var_idx] = update_gt_mat
                # concat 2 indels mean result in no variant, make sure append real variant
                if len(tmp_vars) > 0:
                    new_vars.extend(tmp_vars)
                    new_vars_gt_lst.append(tmp_gt_mat)
            
            # combine all genotypes
            if len(new_vars_gt_lst) > 0:
                new_vars_gt_mat = np.vstack(new_vars_gt_lst)

        else:
            raise ValueError(f"Unknown concatenation method: {method}")
            

        # no variants are concatenated
        if len(new_vars) == 0:
            self.gt_mat = gt_mat
            return 0
        else:
            self.gt_mat = np.vstack([gt_mat, new_vars_gt_mat]) # type: ignore

        for var in new_vars:
            self.add_variant(var)
        assert np.all(self.gt_mat == get_gt_mat(self.var_lst)) # TEST ONLY, should be true

        return len(new_vars)

    
    def merge(self, mis_as_ref: bool, track: str) -> int:
        """ Merge genotypes """

        dup_var_idx = get_duplicates(self.dup_dict)
        while len(dup_var_idx) > 1:
            merge_variants([self.var_lst[i] for i in dup_var_idx],
                           self.gt_mat[dup_var_idx],
                           mis_as_ref,
                           track)
            self.update_merging((self.var_lst[dup_var_idx[0]].alleles))
            dup_var_idx = get_duplicates(self.dup_dict)
        
        return self.n_merge
    

    def update_merging(self, dup_key: tuple):
        """ Update the dup_dict after merging duplicates """

        self.n_merge += len(self.dup_dict[dup_key]) - 1
        self.dup_dict[dup_key] = [self.dup_dict[dup_key][0]]

    
    def get_sorted_idx(self, keep_order: bool) -> list:
        """ Get index of sorted alleles """

        # keep input order
        if keep_order:
            return sorted([x[0] for _, x in self.dup_dict.items()])
        # sort by alleles
        else:
            sorted_keys = sorted(self.dup_dict.keys())
            return [self.dup_dict[k][0] for k in sorted_keys]



def get_duplicates(dup_dict: dict[str, list]) -> list:
    """ Get the index of first duplicated variants """

    for _, dup_var_idx in dup_dict.items():
        if len(dup_var_idx) > 1:
            return dup_var_idx
    
    return []


## functions to concatenate variants
def concat_variants(var_lst: list[pysam.VariantRecord], 
                    gt_mat: np.ndarray, 
                    hap_indexer: HapIndexer,
                    mis_as_ref: bool, 
                    track: str) -> tuple[list[pysam.VariantRecord], np.ndarray, np.ndarray]:
    """
    Concatenate selected variants and update genotypes
    Genotypes of input variants are updated in-place if input if a view not a copy

    Return: 
        1. list of new variant
        2. genotype matrix of new variants
        3. genotype matrix of updated input variants
    """
    
    # gt_mat: n_var x n_hap
    alt_mat = gt_mat == 1

    if mis_as_ref:
        # any non-missing -> non-missing
        non_mis_hap_idx = np.any(gt_mat != -1, axis=0)
    else:
        non_mis_hap_idx = None   # if not mis_as_ref, this should not be used

    alt_sum_arr = alt_mat.sum(axis=0)
    ret_gt_lst = []
    ret_var_lst = []
    concat_n = 0
    while alt_sum_arr.max() > 1:
        # select the one with the most alt alleles to process
        # this ensure that, if other variants have all 1, no additional alt exist
        hap_idx = alt_sum_arr.argmax()
        # find variants with alt
        var_idx = np.where(alt_mat[:, hap_idx])[0]

        # concat genotypes
        new_gt = concat_gt_mat(gt_mat[var_idx], mis_as_ref)
        if mis_as_ref:
            force_ref_idx = (new_gt == -1) & non_mis_hap_idx
            new_gt[force_ref_idx] = 0

        # update genotypes of finished hap as they have been consumed
        reset_hap_idx = np.where(new_gt == 1)[0]
        alt_sum_arr[reset_hap_idx] = 0
        gt_mat[np.ix_(var_idx, reset_hap_idx)] = 0
        # in-place update VariantRecord
        for i in var_idx:
            update_gt_at(var_lst[i], reset_hap_idx.tolist(), hap_indexer, 
                         [0] * reset_hap_idx.shape[0])

        # create new variant
        new_var = create_concat_var([var_lst[i] for i in var_idx], new_gt.tolist(),
                                    concat_n, track)
        if new_var is not None:
            ret_var_lst.append(new_var)
            ret_gt_lst.append(new_gt)
            concat_n += 1
    
    return ret_var_lst, np.array(ret_gt_lst), gt_mat


def get_gt_mat(var_lst: list[pysam.VariantRecord]) -> np.ndarray:
    """ Save genotypes to a matrix """

    gt_mat_lst = []
    for variant in var_lst:
        gt_var_lst = []
        for sample in variant.samples.values():
            for x in sample['GT']:
                if x is None:
                    gt_var_lst.append(-1)
                else:
                    gt_var_lst.append(x)
        gt_mat_lst.append(gt_var_lst)

    return np.array(gt_mat_lst, dtype=np.int8)


def concat_gt_mat(gt_mat: np.ndarray, mis_as_ref: bool) -> np.ndarray:
    """ Get concatenated genotypes from a matrix """

    # all 1 -> 1, any 0 -> 0, else -1
    gt_arr = np.empty(gt_mat.shape[1], dtype=np.int8)
    gt_arr.fill(-1)

    is_alt = np.all(gt_mat == 1, axis=0)
    if not mis_as_ref:
        is_ref = np.any(gt_mat == 0, axis=0)
    # if mis_as_ref, only all -1 -> -1
    else:
        is_ref = np.any(gt_mat != -1, axis=0)
        is_ref[is_alt] = False
    gt_arr[is_alt] = 1
    gt_arr[is_ref] = 0

    return gt_arr


def create_concat_var(var_lst: list[pysam.VariantRecord], gt_lst: list[int],
                      var_id_no: int, track: str) -> pysam.VariantRecord | None:
    """ Concatenate a list of variants to create a new one """

    concat_ref, concat_alt = concat_alleles(var_lst)
    # no variants after merging alleles
    if concat_ref == concat_alt:
        return None
    
    new_var = var_lst[-1].copy()
    new_var.id = f'{new_var.chrom}:{new_var.pos}_{var_id_no}'
    new_var.alleles = (concat_ref, concat_alt)

    new_var.info.clear()
    if track == "ID":
        new_var.info['CONCAT'] = ':'.join([var.id for var in var_lst if var.id is not None])
    elif track == "AT":
        ref_at = ':'.join([var.info['AT'][0] for var in var_lst])
        alt_at = ':'.join([var.info['AT'][1] for var in var_lst])
        new_var.info['CONCAT'] = ref_at + '_' + alt_at

    update_gt_all(new_var, gt_lst)

    return new_var


def concat_alleles(var_lst: list[pysam.VariantRecord]) -> tuple[str, str]:
    """ Concatenate alleles """

    assert len(var_lst) > 1    # TEST ONLY, should always be true

    ref, alt = var_lst[0].alleles   # type: ignore
    ref0 = ref[0]
    ref_trim = ref[1:]

    # # TEST code
    # if len(ref) > 1 and len(alt) > 1:
    #     print(f'Found complex base variant at {var_lst[0].chrom}:{var_lst[0].pos}')

    for i in range(1, len(var_lst)):
        add_ref, add_alt = var_lst[i].alleles  # type: ignore
        # process indels
        if is_indel(add_ref, add_alt):
            add_indel = get_indel_seq(add_ref, add_alt)
            check_indel(ref_trim, add_indel, var_lst[i])

            # deletion
            if len(add_ref) > 1:
                ref_trim = add_indel + ref_trim
            # insertion
            else:
                n_shift = len(ref_trim)
                alt = alt + right_shift_indel(add_indel, n_shift)
        else:
            check_replacement(ref, alt, var_lst[i])
    
    # right-trim
    ref = ref0 + ref_trim
    ref, alt = right_trim_alleles((ref, alt))

    return ref, alt


def check_replacement(ref: str, alt: str, var_add: pysam.VariantRecord) -> None:
    """ Check SNP, MNP, and complex replacement variants """

    ref_add, alt_add = var_add.alleles # type: ignore
    if ref[:len(ref_add)] != ref_add or alt[:len(alt_add)] != alt_add:
        logger.warning(f"Identify haplotype conflict between:\n" + 
                       f"{var_add.chrom}:{var_add.pos}:{ref}_{alt}\n" +
                       f"{var_add.chrom}:{var_add.pos}:{var_add.ref}_{var_add.alts[0]} (force this to reference)") # type: ignore
    else:
        # here we assume all compatible replacements are just redundantly called
        logger.warning(f"Ignore redundant non-indel overlapping:\n" +
                       f"{var_add.chrom}:{var_add.pos}:{ref}_{alt}\n" +
                       f"{var_add.chrom}:{var_add.pos}:{var_add.ref}_{var_add.alts[0]}") # type: ignore


def check_indel(ref_trim: str, indel_seq: str, var_add: pysam.VariantRecord, ) -> None:
    """ Check if a variant can be correctly concatenated with an INDEL """
    
    # check if the second can be right-shifted to the end of base's ref
    indel_len = len(indel_seq)
    ref_trim_len = len(ref_trim)

    if ref_trim_len <= indel_len:
        expect_ref_trim = indel_seq[0:ref_trim_len]
    else:
        n_copy = ref_trim_len // indel_len
        n_shift = ref_trim_len % indel_len
        if n_shift == 0:
            expect_ref_trim = indel_seq * n_copy
        else:
            expect_ref_trim = indel_seq * n_copy + indel_seq[:n_shift]

    if expect_ref_trim != ref_trim:
        raise ValueError(f"{var_add.chrom}:{var_add.pos}:{var_add.ref}_{var_add.alts[0]} " +  # type: ignore
                         f"cannot be right shifted to concatenate with {ref_trim}.")


def right_shift_indel(seq: str, n: int) -> str:
    """ Right shift a sequence n bases """

    n = n % len(seq)
    return seq[n:] + seq[:n]


def right_trim_alleles(alleles: tuple[str, str]) -> tuple:
    """ Left-trim same bases """

    ref_trim = alleles[0]
    alt_trim = alleles[1]

    while len(ref_trim) > 1 and len(alt_trim) > 1 and ref_trim[-1] == alt_trim[-1]:
        ref_trim = ref_trim[:-1]
        alt_trim = alt_trim[:-1]
       
    return ref_trim, alt_trim


def get_indel_seq(ref: str, alt: str) -> str:
    """ Get indel sequence """

    if len(ref) == 1:
        return alt[1:]
    else:
        return ref[1:]


def update_gt_all(var: pysam.VariantRecord, gt_lst: list[int], phased: bool=True):
    """ In-place update genotypes from a list """

    # update genotypes
    hap_idx = 0
    for sample in var.samples.values():
        ploidy = len(sample['GT'])
        new_gt = gt_lst[hap_idx: hap_idx + ploidy]
        new_gt = [GT_DICT[x] for x in new_gt]
        sample['GT'] = tuple(new_gt)
        sample.phased = phased
        hap_idx += ploidy

        
def update_gt_at(var: pysam.VariantRecord, idx_lst: list[int], hap_indexer: HapIndexer,
                 gt_lst: list[int], phased: bool=True):
    """ In-place update genotypes at a specific position """

    for i, gt in zip(idx_lst, gt_lst):  

        samp_idx, ploidy, n_hap = hap_indexer[i]
        
        old_gt = var.samples[samp_idx]['GT']
        new_gt = list(old_gt)
        new_gt[n_hap] = GT_DICT[gt]
        var.samples[samp_idx]['GT'] = tuple(new_gt)
        var.samples[samp_idx].phased = phased


def is_indel(ref: str, alt: str) -> bool:
    """ Check if a variant is an indel """

    if ref[0] != alt[0]:
        return False

    return (len(ref) != 1) != (len(alt) != 1)
    

## end of concat functions




## functions to merge identical variants

def merge_variants(var_lst: list[pysam.VariantRecord], gt_mat: np.ndarray, 
                   mis_as_ref: bool, merge_id: str) -> None:
    """ In-place merge identical variants """

    merged_var = var_lst[0]
    gt_merge = merge_gt_mat(gt_mat, mis_as_ref)
    update_gt_all(merged_var, gt_merge.tolist())

    dup_lst = []
    concat_lst = []

    for var in var_lst:
        if 'CONCAT' in var.info:
            concat_lst.append(var.info['CONCAT'][0])
        elif merge_id == "ID":
            dup_lst.append(var.id)
        elif merge_id == "AT":
            dup_lst.append(var.info['AT'][0] + '_' + var.info['AT'][1])

    if len(dup_lst) > 1:
        merged_var.info['DUP'] = dup_lst[1:]

    if len(concat_lst) > 0:
        merged_var.info['CONCAT'] = concat_lst

    # check if multiple alleles on the same haplotype:
    if np.any(np.sum(gt_mat==1, axis=0) > 1):
        logger.warning(f'More than 1 alleles on the same haplotype for variant: ' +
                       f'{merged_var.chrom}:{merged_var.pos}:{merged_var.ref}_{merged_var.alts[0]}') # type: ignore
            

def merge_gt_mat(gt_mat: np.ndarray, mis_as_ref: bool) -> np.ndarray:
    """ Merge genotype matrix """

    # initialize with all missing
    ret_gt = np.empty(gt_mat.shape[1], dtype=np.int8)
    ret_gt.fill(-1)

    # any 1 -> 1
    is_alt = np.any(gt_mat == 1, axis=0)
    # if mis_as_ref, any 0 -> 0 if not is_alt
    if mis_as_ref:
        is_ref = np.any(gt_mat != -1, axis=0)
        is_ref[is_alt] = False
    # otherwise, all 0 -> 0
    else:
        is_ref = np.all(gt_mat == 0, axis=0)
    ret_gt[is_alt] = 1
    ret_gt[is_ref] = 0

    return ret_gt

## end of merge functions


# O(n) algorithm from https://stackoverflow.com/questions/6021274
def find_rep_motif(allele: str) -> str:
    if not allele:
        return allele

    nxt = [0]*len(allele)
    for i in range(1, len(nxt)):
        k = nxt[i - 1]
        while True:
            if allele[i] == allele[k]:
                nxt[i] = k + 1
                break
            elif k == 0:
                nxt[i] = 0
                break
            else:
                k = nxt[k - 1]

    rep_len = len(allele) - nxt[-1]
    if len(allele) % rep_len != 0:
        return allele

    return allele[0:rep_len]


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


def process_dups(outvcf: pysam.VariantFile, var_dict: DuplicatesDict, prev_var: pysam.VariantRecord, 
                 hap_indexer: HapIndexer, counter: VariantCounter, concat_method: str, 
                 mis_as_ref: bool, track: str, max_repeat: int, keep_order: bool) -> None:
    """ Merge duplicated records and write to output VCF """

    # only 1 variant at the previous position
    n_input = len(var_dict)
    if n_input == 0:
        outvcf.write(prev_var)
        counter.output += 1
        return None

    # merge if more than 1 varaints at the previous position
    n_concat = var_dict.concatenate(hap_indexer=hap_indexer, 
                                    method=concat_method,
                                    mis_as_ref=mis_as_ref,
                                    track=track,
                                    max_repeat=max_repeat)
    n_merge =var_dict.merge(mis_as_ref=mis_as_ref, 
                            track=track)
    logger.debug(f'{prev_var.chrom}:{prev_var.pos}: read {n_input}, concatenate {n_concat}, merge {n_merge}')

    counter.concat += n_concat
    counter.merge += n_merge


    for k in var_dict.get_sorted_idx(keep_order):
        out_var = var_dict.var_lst[k]
        out_var = update_ac(out_var)
        outvcf.write(out_var)
        counter.output += 1


def add_header(header: pysam.VariantHeader) -> pysam.VariantHeader:
    """ Add required headers """
    
    header.add_line('##INFO=<ID=DUP,Number=.,Type=String,Description="List of variant ID or allele traversal merged into this variant">')
    header.add_line('##INFO=<ID=CONCAT,Number=.,Type=String,Description="List of variant ID or allele traversal concatenated into this variant">')

    return header


def check_vcf(vcf: pysam.VariantFile, n: int) -> None:
    """ Check VCF format """

    for i, var in enumerate(vcf.fetch()):
        if i > n:
            return

        for sample in var.samples.values():
            ploidy = len(sample['GT'])
            if ploidy >= 2 and not sample.phased:
                raise ValueError(f'Unphased sample is not supported')


def main(args: argparse.Namespace):

    # setup logger
    if args.debug:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO
    logging.basicConfig(level=log_level,
                        format='[%(asctime)s] - [%(levelname)s]: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    
    # initialize
    invcf = pysam.VariantFile(args.invcf, 'rb')
    check_vcf(invcf, n=1000)
    if args.track is None:
        header = invcf.header
    else:
        header = add_header(invcf.header)
    outvcf = pysam.VariantFile(args.outvcf, 'w', header=header)
    counter = VariantCounter()
    hap_indexer = HapIndexer(next(invcf.fetch()))

    logger.info(f'Merge duplicated variants from: {args.invcf}')
    if args.concat == "none":
        logger.info("Skip variant concatenation")
    logger.info(f'Concatenate variants with same {args.concat}')

    # to check duplicates, we always write variants after reading its next one
    # if cur_var.pos = prev_var.pos, save them to working_var
    # if cur_var.pos > prev_var.pos, write all variants in prev_var.pos
    working_var_dict = DuplicatesDict()
    invcf_iter = invcf.fetch()
    prev_var = next(invcf_iter)
    counter.input += 1

    for cur_var in invcf_iter:
        counter.input += 1

        # process all variants at the previous position
        if cur_var.pos > prev_var.pos:
            process_dups(outvcf, working_var_dict, prev_var, hap_indexer, counter,
                         args.concat, args.merge_mis_as_ref, args.track, 
                         args.max_repeat, args.keep_order)
            working_var_dict = DuplicatesDict()
        # store variants at the same position
        elif cur_var.pos == prev_var.pos:
            # for the first 2 variants, store both
            if len(working_var_dict) == 0:
                working_var_dict.add_variant(prev_var)
            working_var_dict.add_variant(cur_var)
        else:
            # end of current chromosome, write all variants on the previous chromosome
            if prev_var.chrom != cur_var.chrom:
                process_dups(outvcf, working_var_dict, prev_var, hap_indexer, counter, 
                             args.concat, args.merge_mis_as_ref, args.track, 
                             args.max_repeat, args.keep_order)
                working_var_dict = DuplicatesDict()
            # only unsorted VCF will have cur_pos < prev_pos, report error and exit
            else:
                raise ValueError('Input VCF is not sorted')

        prev_var = cur_var.copy()
                
    # process at the end of the VCF
    process_dups(outvcf, working_var_dict, prev_var, hap_indexer, counter, 
                 args.concat, args.merge_mis_as_ref, args.track, 
                 args.max_repeat, args.keep_order)

    logger.info(f'Read {counter.input} variants')
    logger.info(f'{counter.merge} variants are merged')
    logger.info(f'{counter.concat} variants are concatenated')
    # if counter.conflict_missing > 0:
    #     logger.info(f'{counter.conflict_missing} variants have non-missing genotypes merged with missing genotypes')

    logger.info(f'Write {counter.output} variants to {args.outvcf}')

    invcf.close()
    outvcf.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='merge_duplicates.py', description='Merge duplicated variants in phased VCF')
    parser.add_argument('-i', '--invcf', metavar='VCF', help='Input VCF, sorted and phased', required=True)
    parser.add_argument('-o', '--outvcf', metavar='VCF', help='Output VCF', required=True)
    parser.add_argument('-c', '--concat', choices=['position', 'repeat', 'none'], default='position',
                        help='Concatenate variants when they have identical "position" (default) or "repeat" motif, "none" to skip')
    parser.add_argument('-m', '--max-repeat', default=None, type=int,
                        help='Maximum size a variant to search for repeat motif (default: None)')
    parser.add_argument('-t', '--track', choices=['ID', 'AT'], default=None,
                        help='Track how variants are merged by "ID", "AT", or disable (default)')
    parser.add_argument('--merge-mis-as-ref', action='store_true', 
                        help='Convert missing to ref when merging missing genotypes with non-missing genotypes')
    parser.add_argument('--keep-order', action='store_true',
                        help='keep the order of variants in the input VCF (default: sort by chr, pos, alleles)')
    parser.add_argument('--debug', action='store_true', 
                        help='Debug mode')

    args = parser.parse_args()

    main(args)
