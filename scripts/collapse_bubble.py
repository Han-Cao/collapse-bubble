#!/usr/bin/env python3

# SV merging within the same bubble
# Last update: 10-Jun-2025
# Author: Han Cao

import logging
import argparse
from collections import defaultdict

import pysam
import truvari
import numpy as np
import pandas as pd


def get_variant_type(variant: pysam.VariantRecord, min_sv_len: int, sv_only: bool=True) -> str:
    """ Return variant type """

    len_ref = len(variant.ref)
    len_alt = len(variant.alts[0])

    if max(len_ref, len_alt) < min_sv_len:
        if sv_only:
            return None
        if len_ref == 1 and len_alt == 1:
            return 'SNP'
        elif len_ref == len_alt:
            return 'MNP'
    
    # check INV annotation from vcfwave
    if 'INV' in variant.info:
        return 'INV'
    # get SVTYPE from allele length
    if len_ref == 1 and len_alt > 1:
        return 'INS'
    elif len_ref > 1 and len_alt == 1:
        return 'DEL'
    else:
        return 'COMPLEX'
    

def clean_bubble_dict(bubble_dict: dict) -> None:
    """ Clean and annotate the bubble dict in place """

    drop_bubbles = []
    for k, v in bubble_dict.items():
        # count SVTYPE
        sv_types, sv_counts = np.unique(v['type'], return_counts=True)
        # remove bubble with only one SV
        if sv_counts.max() == 1:
            drop_bubbles.append(k)
            continue
        # save counts
        bubble_dict[k]['counts'] = {t: c for t, c in zip(sv_types, sv_counts)}
    
    for k in drop_bubbles:
        del bubble_dict[k]
    
    return bubble_dict


def get_vcf_iter(vcf: pysam.VariantFile, chr: str) -> pysam.VariantFile:
    if chr is not None:
        vcf_iter = vcf.fetch(contig=chr)
    else:
        vcf_iter = vcf.fetch()
    
    return vcf_iter


def parse_vcf(vcf: pysam.VariantFile, chr: str, min_len: int) -> dict:
    """ Read through VCF and save bubbles need to be collapsed """

    logger = logging.getLogger(__name__)
    # dict to store SVs
    # bubble_id -> {'id': [], 'type': []}
    bubble_dict = {}

    vcf_iter = get_vcf_iter(vcf, chr)
    n_var = 0
    n_sv = 0
    for variant in vcf_iter:
        n_var += 1
        # check if allele length is >= min_len
        svtype = get_variant_type(variant, min_len)
        if svtype is None:
            continue
        n_sv += 1
        
        bubble_id = variant.info['BUBBLE_ID']
        if bubble_id not in bubble_dict:
            bubble_dict[bubble_id] = {'id': [], 'type': []}

        # make sure id is unique
        if variant.id in bubble_dict[bubble_id]['id']:
            raise ValueError(f'Duplicate ID found: {variant.id}, please annotate variant with unique ID first')
        # store SV into bubble_dict
        bubble_dict[bubble_id]['id'].append(variant.id)
        bubble_dict[bubble_id]['type'].append(svtype)
    
    logger.info(f'Read {n_var} variants')
    logger.info(f'Found {n_sv} SVs >= {min_len} bp')
    # clean and annotate
    clean_bubble_dict(bubble_dict)
    logger.info(f'Found {len(bubble_dict)} bubbles with more than 1 SVs')

    return bubble_dict


def retrive_svtype(variant: pysam.VariantRecord, bubble_dict: dict) -> str:
    """ Retrive SVTYPE from bubble dict """

    bubble = bubble_dict[variant.info['BUBBLE_ID']]
    idx = bubble['id'].index(variant.id) if variant.id in bubble['id'] else None

    if idx is None:
        return None

    return bubble['type'][idx]


def annotate_sv(variant: pysam.VariantRecord, svtype: str) -> pysam.VariantRecord:
    """ Annotate SV based on SVTYPE """

    new_var = variant.copy()

    new_var.info['SVTYPE'] = svtype
    # truvari always use INFO/SVLEN to compare size, so we set it differently
    # this will not affect SVLEN in output VCF
    # INS, DEL: len(ref) - len(alt)
    # INV, COMPLEX: len(ref)
    if svtype == 'COMPLEX' or svtype == 'INV':
        new_var.info['SVLEN'] = len(variant.ref)
    else:
        new_var.info['SVLEN'] = len(variant.alts[0]) - len(variant.ref)

    return new_var


# functions to convert genotypes to binary code
def is_none(x) -> bool:
    return x is None
def has_var(x) -> bool:
    return x == 1


def hap_conflict(var1: pysam.VariantRecord, var2: pysam.VariantRecord) -> bool:
    """ Check whether the haplotypes are consistent if two SVs are merged """

    # # test only
    # # check if genotype always co-missing on the same haplotype
    # # this is expected for varaints decomposed from the same bubble
    # var1_gt = [is_none(x) for sample in var1.samples.values() for x in sample['GT']]
    # var1_gt = np.array(var1_gt, dtype=bool)
    # var2_gt = [is_none(x) for sample in var2.samples.values() for x in sample['GT']]
    # var2_gt = np.array(var2_gt, dtype=bool)
    # assert np.array_equal(var1_gt, var2_gt)

    # check haplotype consistency
    var1_gt = [has_var(x) for sample in var1.samples.values() for x in sample['GT']]
    var1_gt = np.array(var1_gt, dtype=bool)
    var2_gt = [has_var(x) for sample in var2.samples.values() for x in sample['GT']]
    var2_gt = np.array(var2_gt, dtype=bool)

    conflict = (var1_gt & var2_gt).any()

    return conflict


def mac(variant: pysam.VariantRecord) -> int:
    """ Return minor allele count based on AC, AN """

    alt_count = variant.info['AC'][0] if isinstance(variant.info['AC'], tuple) else variant.info['AC']
    ref_count = variant.info['AN'] - alt_count

    return min(alt_count, ref_count)


def collapse_bubble(var_lst: list, matcher: truvari.Matcher) -> dict:
    """
    Collapse SVs within the same bubble
    Return the map from original SV ID to collapsed SV ID
    """

    id_map = {} # SV ID -> Collapsed SV ID
    collapsed_sv = defaultdict(list) # Collapsed SV ID -> list of SV records

    # start from the most frequent SVs
    var_remain = sorted(var_lst, key=mac, reverse=True)
    while len(var_remain) > 0:
        collapse_var = var_remain.pop(0)
        # always include itself
        id_map[collapse_var.id] = collapse_var.id
        collapsed_sv[collapse_var.id].append(collapse_var)

        drop_var = []
        # SV comparison
        for other_var in var_remain:
            res_match = matcher.build_match(collapse_var, other_var, skip_gt=True, short_circuit=True)
            if res_match.state:
                # check haplotype consistency against all SVs already collapsed
                flag_conflict = False
                for check_var in collapsed_sv[collapse_var.id]:
                    if hap_conflict(check_var, other_var):
                        flag_conflict = True
                        break
                if flag_conflict:
                    continue

                # merge SVs
                id_map[other_var.id] = collapse_var.id
                collapsed_sv[collapse_var.id].append(other_var)
                drop_var.append(other_var)
        
        var_remain = [x for x in var_remain if x not in drop_var]
           
    return id_map


def get_variant_info(variant: pysam.VariantRecord, collapse_id: str, info_lst: list) -> dict:
    """ Save variant info into dict """

    output = {'CHROM': variant.chrom,
              'POS': variant.pos,
              'Bubble_ID': variant.info['BUBBLE_ID'],
              'Variant_ID': variant.id,
              'Collapse_ID': collapse_id
            }

    if len(info_lst) > 0:
        for tag in info_lst:
            tag_value = variant.info[tag]
            # extract value for single element tuple
            if isinstance(tag_value, tuple) and len(tag_value) == 1:
                tag_value = tag_value[0]
            output[f'INFO_{tag}'] = tag_value

    return output


def collapse_vcf(vcf: pysam.TabixIterator, 
                 matcher: truvari.Matcher, 
                 bubble_dict: dict, 
                 info_lst: list) -> pd.DataFrame:
    """ Collapse bubbles in the VCF """
    
    working_bubbles = defaultdict(list) # bubble_id -> sv_id
    variant_map = []

    for variant in vcf:
        bubble_id = variant.info['BUBBLE_ID']
        # bubble to skip as it does not has SV to collapse
        if bubble_id not in bubble_dict:
            continue

        svtype = retrive_svtype(variant, bubble_dict)
        # not a SV
        if svtype is None:
            continue
        
        working_bubbles[bubble_id].append(annotate_sv(variant, svtype))

        # collapse if all SVs in this bubble have been processed
        n_working = len(working_bubbles[bubble_id])
        n_total = len(bubble_dict[bubble_id]['id'])

        if n_working == n_total:
            bubble_map = collapse_bubble(working_bubbles[bubble_id], matcher)
            for sv in working_bubbles[bubble_id]:
                variant_map.append(get_variant_info(sv, bubble_map[sv.id], info_lst))
            # free memory
            del working_bubbles[bubble_id]
    
    # If everthing works well, all bubbles should be removed from working_bubble
    assert len(working_bubbles) == 0

    # convert to dataframe, remove SVs will not be merged
    df_collapse = pd.DataFrame.from_dict(variant_map)
    sv_counts = df_collapse.value_counts('Collapse_ID')
    keep_sv = sv_counts[sv_counts > 1].index
    df_collapse = df_collapse.loc[df_collapse['Collapse_ID'].isin(keep_sv)].reset_index(drop=True)

    # create dataframe
    return df_collapse


def add_header(header: pysam.VariantHeader) -> pysam.VariantHeader:
    """ Add required INFO to header  """

    if 'TYPE' not in header.info:
        header.add_line('##INFO=<ID=TYPE,Number=A,Type=String,Description="Type of variant">')
    
    header.add_line('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">')
    header.add_line('##INFO=<ID=REFLEN,Number=1,Type=Integer,Description="Length of REF allele">')
    header.add_line('##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">')
    header.add_line('##INFO=<ID=ID_LIST,Number=.,Type=String,Description="List of SVs collapsed into this SV">')

    return header


GT_DICT = {0: 0,
           1: 1,
           -1: None}


def collapse_genotype(var_lst: list, var_collapse: pysam.VariantRecord) -> pysam.VariantRecord:
    """ Merge genotypes from list of SVs """

    missing_lst = []
    has_var_lst = []

    for variant in var_lst:
        has_var_lst.append([has_var(x) for sample in variant.samples.values() for x in sample['GT']])
        missing_lst.append([is_none(x) for sample in variant.samples.values() for x in sample['GT']])
    
    has_var_arr = np.array(has_var_lst)
    missing_arr = np.array(missing_lst)

    # merge genotypes
    # has any SV -> has collapsed SV
    merge_has_var = has_var_arr.any(axis=0)
    # all missing -> missing
    # actually for MC output VCF, 1 missing GT always = all missing
    merge_missing = missing_arr.all(axis=0)
    # if missing -> -1, else, 0 or 1
    merge_gt = merge_has_var.astype(int)
    merge_gt[merge_missing] = -1

    # check only: 
    # for one haplotype, maximum one SV
    if has_var_arr.sum(axis=0).max() > 1:
        logger = logging.getLogger(__name__)
        logger.warning(f'More than 1 SV on the same haplotype, collapse SV ID: {var_collapse.id}')
    # . and 1 should not exist on the same haplotype
    if (merge_has_var & merge_missing).any():
        logger = logging.getLogger(__name__)
        logger.warning(f'Inconsistent genotypes (. vs 1) on the same haplotype, collapse SV ID: {var_collapse.id}')
    
    # update genotypes
    # TODO: use update_gt_all, can support any ploidy
    hap_idx = 0
    merge_gt = merge_gt.tolist()
    for sample in var_collapse.samples.values():
        ploidy = len(sample['GT'])
        if ploidy == 1:
            sample['GT'] = (GT_DICT[merge_gt[hap_idx]],)
            hap_idx += 1
        elif ploidy == 2:
            sample['GT'] = (GT_DICT[merge_gt[hap_idx]], GT_DICT[merge_gt[hap_idx+1]])
            sample.phased = True # keep the phase
            hap_idx += 2
        else:
            raise ValueError('Ploidy > 2 not supported')
    
    return var_collapse


# TODO: we should update AC, AN, AF of merged VCF
def write_outvcf(invcf: pysam.VariantFile, outvcf: pysam.VariantFile,
                 df_collapse: pd.DataFrame, chr: str, min_len: int) -> None:
    """ Convert SVs and write to output VCF """

    invcf_iter = get_vcf_iter(invcf, chr)

    # extract ID mapping dict from df
    id_map = df_collapse[['Variant_ID', 'Collapse_ID']].set_index('Variant_ID')['Collapse_ID'].to_dict()
    # get summary of collapsed SV
    df_summary = df_collapse.groupby('Collapse_ID').agg({'CHROM': 'first', 'POS': 'first', 'Variant_ID': 'count'})
    df_summary = df_summary.rename(columns={'Variant_ID': 'count'})

    # cache variants
    working_collapse = defaultdict(list) # collapse SV ID -> list of SVs to be merged

    for variant in invcf_iter:
        var_out = variant.copy()
        var_type = get_variant_type(variant, min_len, sv_only=False)
        var_out.info['TYPE'] = var_type

        # nothing to do with SNP
        if var_type == 'SNP':
            outvcf.write(var_out)
            continue
        # MNP only need REFLEN
        var_out.info['REFLEN'] = len(variant.ref)
        if var_type == 'MNP':
            outvcf.write(var_out)
            continue
        # For SVs, add SVLEN and SVTYPE, check if need to collapse
        var_out.info['SVTYPE'] = var_type
        var_out.info['SVLEN'] = len(variant.alts[0]) - len(variant.ref)
        if variant.id not in id_map:
            outvcf.write(var_out)
            continue

        # collapse SVs
        collapse_id = id_map[variant.id]
        working_collapse[collapse_id].append(var_out)

        # process if all SVs have been collected
        collapse_info = df_summary.loc[collapse_id]
        if len(working_collapse[collapse_id]) == collapse_info['count']:
            # find target SV
            for var_collapse in working_collapse[collapse_id]:
                if var_collapse.id == collapse_id:
                    break
            # modify collapsed SV genotypes and INFO
            var_collapse = collapse_genotype(working_collapse[collapse_id], var_collapse)
            var_collapse.info['ID_LIST'] = ','.join([x.id for x in working_collapse[collapse_id]])
            outvcf.write(var_collapse)
            # free memory
            del working_collapse[collapse_id]
    
    # if all SVs are processed, working_collapse should be empty
    assert len(working_collapse) == 0


def main(args: argparse.Namespace):

    # setup logger
    logging.basicConfig(level=logging.INFO,
                        format='[%(asctime)s] - [%(levelname)s]: %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    logger = logging.getLogger(__name__)
    # set up matcher
    matcher = truvari.Matcher()
    matcher.params.typeignore = False
    matcher.params.refdist = args.refdist
    matcher.params.pctseq = args.pctseq
    matcher.params.pctsize = args.pctsize
    matcher.params.pctovl = args.pctovl
    logger.info(f'Setup Truvari Matcher: refdist={args.refdist}, pctseq={args.pctseq}, pctsize={args.pctsize}, pctovl={args.pctovl}')

    # parse input VCF
    logger.info(f'Read input VCF: {args.invcf}')
    invcf = pysam.VariantFile(args.invcf, 'rb')
    bubble_dict = parse_vcf(invcf, args.chr, args.min_len)
    new_header = add_header(invcf.header)

    # collapse VCF
    invcf_iter = get_vcf_iter(invcf, args.chr)
    info_lst = [] if args.info is None else args.info.split(',')

    df_collapse = collapse_vcf(invcf_iter, matcher, bubble_dict, info_lst)
    logger.info(f'Collapse {len(df_collapse)} SVs into {df_collapse["Collapse_ID"].unique().shape[0]} SVs')

    # write collapse map
    df_collapse.to_csv(args.map, sep='\t', index=False)
    logger.info(f'Write collapse map to {args.map}')

    # write output VCF
    outvcf = pysam.VariantFile(args.outvcf, 'w', header=new_header)

    write_outvcf(invcf, outvcf, df_collapse, args.chr, args.min_len)
    logger.info(f'Write output VCF to {args.outvcf}')

    invcf.close()
    outvcf.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='collapse_bubble.py', description='Collapse biallelic SVs within the same bubble in VCF')
    io_arg = parser.add_argument_group('Input / Output arguments')
    io_arg.add_argument('-i', '--invcf', metavar='VCF', help='Input VCF', required=True)
    io_arg.add_argument('-o', '--outvcf', metavar='VCF', help='Output VCF', required=True)
    io_arg.add_argument('-m', '--map', metavar='TSV', help='Write SV mapping table to this file. Default: None', type=str, required=True)
    io_arg.add_argument('--chr', metavar='CHR', help='chromosome to work on. Default: all', type=str, default=None)
    io_arg.add_argument('--info', metavar='TAG', help='Comma-separated INFO/TAG list to include in the output map. Default: None', type=str, default=None)

    collapse_arg = parser.add_argument_group('Collapse arguments')
    collapse_arg.add_argument('-l', '--min-len', metavar='50', help='Minimum allele length of variants to be included, defined as max(len(alt), len(ref)). Default: 50', type=int, default=50)
    collapse_arg.add_argument('-r', '--refdist', metavar='100', help='Max reference location distance. Default: 100', type=int, default=100)
    collapse_arg.add_argument('-p', '--pctseq', metavar='0.9', help='Min percent sequence similarity (REF for DEL, ALT for other SVs). Default: 0.9', type=float, default=0.9)
    collapse_arg.add_argument('-P', '--pctsize', metavar='0.9', help='Min percent size similarity (SVLEN for INS, DEL; REFLEN for INV, COMPLEX). Default: 0.9', type=float, default=0.9)
    collapse_arg.add_argument('-O', '--pctovl', metavar='0.9', help='Min pct reciprocal overlap. Default: 0.9', type=float, default=0.9)

    args = parser.parse_args()
    main(args)
