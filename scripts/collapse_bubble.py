#!/usr/bin/env python3

# SV merging for pangenome VCF
# Last update: 30-Jun-2025
# Author: Han Cao

import logging
import argparse
from collections import defaultdict

import pysam
import truvari
import numpy as np
import pandas as pd


class BubbleClusters:
    """ Clusters of bubbles overlapping with the same tandem repeat """

    def __init__(self) -> None:

        self.bubbles = {} # bubble_id -> {'id': [var_ids], 'svtype': [var_svtypes], 'cluster': cluster_id}
        self.clusters = defaultdict(set) # cluster_id -> set of bubble_ids
        self.cluster_vars = defaultdict(list) # cluster_id -> list of var_ids
        self.overlaps = defaultdict(set) # bubble_id -> set of overlapping bubbles (self excluded)
        
        self._chr = None # current chr
        self._pos = 0 # current position
        self._pos_bubbles = set() # all bubbles in the current position


    def get_bubble_cluster(self, bubble_id: str) -> int:
        """ Get cluster id for bubble """

        if bubble_id in self.bubbles:
            return self.bubbles[bubble_id]['cluster']
        else:
            return None
    

    def get_cluster_variants(self, cluster_id: int) -> list:
        """ Get variants in a cluster """

        return self.cluster_vars[cluster_id]
     

    def get_svtype(self, variant: pysam.VariantRecord) -> str:
        """ Get SVTYPE of a variant """

        bubble_id = get_bubble_ids(variant)[0]
        if variant.id in self.bubbles[bubble_id]['id']:
            idx = self.bubbles[bubble_id]['id'].index(variant.id)
        # non-SV are not included in self.bubbles, return None
        else:
            return None
        
        return self.bubbles[bubble_id]['svtype'][idx]
    

    def _update_overlap(self) -> None:
        """ Update overlaps after position change """

        for b1 in self._pos_bubbles:
            for b2 in self._pos_bubbles:
                if b1 != b2:
                    self.overlaps[b1].add(b2)
    
    
    def add_variant(self, variant: pysam.VariantRecord, svtype: str) -> None:
        """ Add variant to bubble and check for overlaps """

        # get bubble id
        bubble_ids = get_bubble_ids(variant)
        
        # collect bubbles at the same position
        if variant.pos == self._pos and variant.chrom == self._chr:
            self._pos_bubbles.update(bubble_ids)

        # all variants at the previous position has been processed
        else:
            # update overlaps 
            self._update_overlap()
            # update position
            self._chr = variant.chrom
            self._pos = variant.pos
            self._pos_bubbles = set(bubble_ids)
        
        # skip non-SV ensures bubbles without any SV not in self.bubbles
        if svtype is None:
            return
        
        # add SV to bubble
        # Note: concatenated variants belong to all the original bubbles
        for b_id in bubble_ids:
            if b_id not in self.bubbles:
                self.bubbles[b_id] = {'id': [variant.id], 'svtype': [svtype], 'cluster': None}
            else:
                # check unique id
                if variant.id in self.bubbles[b_id]['id']:
                    raise ValueError(f'Duplicate variant ID {variant.id} in bubble {b_id}')
                self.bubbles[b_id]['id'].append(variant.id)
                self.bubbles[b_id]['svtype'].append(svtype)

    
    def cluster_bubbles(self) -> None:
        """ Cluster bubbles """

        cluster_id = 0
        for b_id in self.bubbles.keys():

            # find overlap bubbles with SVs
            overlap_bubbles = [x for x in self.overlaps[b_id] if x in self.bubbles]
            
            # singleton bubble
            if len(overlap_bubbles) == 0:
                self.bubbles[b_id]['cluster'] = cluster_id
                self.clusters[cluster_id].add(b_id)
                cluster_id += 1
            
            # multiple bubble cluster
            else:
                # find if any bubble already clustered
                exist_cluster = []
                cluster_bubbles = [b_id] + overlap_bubbles
                for id in cluster_bubbles:
                    if self.bubbles[id]['cluster'] is not None:
                        exist_cluster.append(self.bubbles[id]['cluster'])

                # create a new cluster
                unique_cluster = list(set(exist_cluster))
                if len(unique_cluster) == 0:
                    for id in cluster_bubbles:
                        self.bubbles[id]['cluster'] = cluster_id
                    self.clusters[cluster_id].update(cluster_bubbles)
                    cluster_id += 1
                
                # add to existing cluster
                else:
                    assert len(unique_cluster) == 1
                    for id in cluster_bubbles:
                        self.bubbles[id]['cluster'] = unique_cluster[0]
                    self.clusters[unique_cluster[0]].update(cluster_bubbles)
        
        # test only: check all bubbles have been clustered
        for b_id in self.bubbles.keys():
            assert self.bubbles[b_id]['cluster'] is not None


    def finalize_clusters(self) -> None:
        """ Remove clusters with only one SV per type, prepare variant list per cluster """

        drop_clusters = []
        # identify clusters to drop
        for c_id in self.clusters:
            svtypes = []
            var_ids = []
            for b_id in self.clusters[c_id]:
                svtypes += self.bubbles[b_id]['svtype']
                var_ids += self.bubbles[b_id]['id']
            _, sv_counts = np.unique(svtypes, return_counts=True)

            # if not remove, save var_ids
            if sv_counts.max() == 1:
                drop_clusters.append(c_id)
            else:
                # make sure we only save unique var_ids
                # keep the order to make the output more deterministic
                self.cluster_vars[c_id] = list(dict.fromkeys(var_ids))
        
        # remove cluster and related bubbles
        for c_id in drop_clusters:
            for b_id in self.clusters[c_id]:
                del self.bubbles[b_id]
            del self.clusters[c_id]


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
    

def get_bubble_ids(variant: pysam.VariantRecord) -> list[str]:
    """ Get list bubble ID of a variant """

    if 'BUBBLE_ID' in variant.info:
        bubble_ids = [variant.info['BUBBLE_ID']]
    elif 'CONCAT' in variant.info:
        concat_var_ids = variant.info['CONCAT'][0].split(':')
        # concat id in format of bubble_id.type.number
        bubble_ids = [x.rsplit('.', 2)[0] for x in concat_var_ids]
        # remove duplicates, this happens when variants from the same bubble are concat
        bubble_ids = list(set(bubble_ids))
    else:
        raise ValueError(f'Cannot find INFO/BUBBLE_ID or INFO/CONCAT in variant {variant.id}')
    
    return bubble_ids
    

def get_vcf_iter(vcf: truvari.VariantFile, chr: str) -> truvari.VariantFile:
    if chr is not None:
        vcf_iter = vcf.fetch(contig=chr)
    else:
        vcf_iter = vcf.fetch()
    
    return vcf_iter


def parse_vcf(vcf: truvari.VariantFile, chr: str, min_len: int) -> BubbleClusters:
    """ Read through VCF and save bubbles need to be collapsed """

    logger = logging.getLogger(__name__)

    bubble_clusters = BubbleClusters()

    vcf_iter = get_vcf_iter(vcf, chr)
    n_var = 0
    n_sv = 0

    for variant in vcf_iter:
        n_var += 1
        # find sv by allele length >= min_len
        svtype = get_variant_type(variant, min_len)
        bubble_clusters.add_variant(variant, svtype)
        if svtype is not None:
            n_sv += 1
    
    # update overlap for the last position
    bubble_clusters._update_overlap()
        
    logger.info(f'Read {n_var} variants')
    logger.info(f'Found {n_sv} SVs >= {min_len} bp')

    # find cluster and clean up
    bubble_clusters.cluster_bubbles()
    bubble_clusters.finalize_clusters()

    logger.info(f'Found {len(bubble_clusters.clusters)} bubble clusters with more than 1 SVs per SVTYPE')

    return bubble_clusters


def annotate_sv(variant: pysam.VariantRecord, svtype: str) -> pysam.VariantRecord:
    """ Annotate SV based on SVTYPE """

    new_var = variant.copy()

    new_var.info['SVTYPE'] = svtype
    # truvari always use INFO/SVLEN to compare size, so we set it differently:
    # INS, DEL: len(ref) - len(alt)
    # INV, COMPLEX: len(ref)
    # this will not affect SVLEN in output VCF
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


def collapse_bubble(var_lst: list[truvari.VariantRecord]) -> tuple[dict, list]:
    """ Collapse SVs within the same cluster, return the matching results """

    match_map = {} # SV ID -> {Collapsed SV ID, Matching stats}
    conflict_lst = [] # list of conflicting SVs {'Variant_ID', 'Conflict_ID', 'Collapse_ID'}
    collapsed_sv = defaultdict(list) # Collapsed SV ID -> list of SV records

    # start from the most frequent SVs
    var_remain = sorted(var_lst, key=mac, reverse=True)
    while len(var_remain) > 0:
        collapse_var = var_remain.pop(0)
        # always include itself
        match_map[collapse_var.id] = {'Collapse_ID': collapse_var.id}
        collapsed_sv[collapse_var.id].append(collapse_var)

        drop_var = []
        # SV comparison
        # TODO: If there are too many SVs, first filter by distance can save time
        # e.g., var_remain_near = [var for var in var_remain if near()]
        for other_var in var_remain:
            res_match = collapse_var.match(other_var)
            if res_match.state:
                # check haplotype consistency against all SVs already collapsed
                # TODO: this can speed up by first merging genotypes and check all genotypes at once
                # To achieve this, we have to do SV merging and genotype merging in one loop
                # Will work on this if run time is an issue
                flag_conflict = False
                for check_var in collapsed_sv[collapse_var.id]:
                    if hap_conflict(check_var, other_var):
                        flag_conflict = True
                        conflict_lst.append({'Variant_ID': other_var.id, 'Conflict_ID': check_var.id, 'Collapse_ID': collapse_var.id})
                        break
                if flag_conflict:
                    continue

                # merge SVs
                match_map[other_var.id] = {'Collapse_ID': collapse_var.id,
                                           'PctSeqSimilarity': res_match.seqsim,
                                           'PctSizeSimilarity': res_match.sizesim,
                                           'PctRecOverlap': res_match.ovlpct,
                                           'SizeDiff': res_match.sizediff,
                                           'StartDistance': res_match.st_dist,
                                           'EndDistance': res_match.ed_dist,
                                           'TruScore': res_match.score
                                           }
                collapsed_sv[collapse_var.id].append(other_var)
                drop_var.append(other_var)
        
        var_remain = [x for x in var_remain if x not in drop_var]
           
    return match_map, conflict_lst


def get_variant_info(variant: pysam.VariantRecord, res_match_dict: dict, info_lst: list) -> dict:
    """ Save variant info into dict """

    bubble_ids = get_bubble_ids(variant)
    if len(bubble_ids) == 1:
        bubble_id = bubble_ids[0]
    else:
        bubble_id = ','.join(bubble_ids)

    output = {'CHROM': variant.chrom,
              'POS': variant.pos,
              'Bubble_ID': bubble_id,
              'Variant_ID': variant.id,
            }
    output.update(res_match_dict)

    if len(info_lst) > 0:
        for tag in info_lst:
            tag_value = variant.info[tag]
            # extract value for single element tuple
            if isinstance(tag_value, tuple) and len(tag_value) == 1:
                tag_value = tag_value[0]
            output[f'INFO_{tag}'] = tag_value

    return output


def collapse_vcf(vcf: truvari.VariantFile, 
                 bubble_clusters: BubbleClusters, 
                 info_lst: list) -> tuple[pd.DataFrame, pd.DataFrame]:
    """ Collapse bubbles in the VCF """
    
    working_clusters = defaultdict(list) # cluster_id -> list of VariantRecord
    collapse_map = []
    conflict_map = []

    for variant in vcf:
        bubble_ids = get_bubble_ids(variant)
        cluster_id = bubble_clusters.get_bubble_cluster(bubble_ids[0])
        # bubble without cluster means it doesn't need to be collapsed
        if cluster_id is None:
            continue

        svtype = bubble_clusters.get_svtype(variant)
        # not a SV
        if svtype is None:
            continue
        
        # cache annotated SVTYPE and SVLEN for comparison
        working_clusters[cluster_id].append(annotate_sv(variant, svtype))

        # collapse if all SVs in this cluster have been processed
        n_working = len(working_clusters[cluster_id])
        cluster_var_ids = bubble_clusters.get_cluster_variants(cluster_id)
        n_total = len(cluster_var_ids)

        if n_working == n_total:
            cluster_map, conflict_lst = collapse_bubble(working_clusters[cluster_id])
            conflict_map += conflict_lst
            
            for sv in working_clusters[cluster_id]:
                assert sv.id in cluster_var_ids # test only: comapre working varinats and those in the cluster
                collapse_map.append(get_variant_info(sv, cluster_map[sv.id], info_lst))

            # free memory
            del working_clusters[cluster_id]
    
    # If everthing works well, all clusters should be removed from working_clusters
    assert len(working_clusters) == 0

    # convert to dataframe, remove SVs will not be merged
    df_conflict = pd.DataFrame(conflict_map)
    df_collapse = pd.DataFrame.from_dict(collapse_map)
    sv_counts = df_collapse.value_counts('Collapse_ID')
    keep_sv = sv_counts[sv_counts > 1].index
    df_collapse = df_collapse.loc[df_collapse['Collapse_ID'].isin(keep_sv)].reset_index(drop=True)

    # create dataframe
    return df_collapse, df_conflict


def add_header(header: pysam.VariantHeader) -> pysam.VariantHeader:
    """ Add required INFO to header  """

    if 'TYPE' not in header.info:
        header.add_line('##INFO=<ID=TYPE,Number=.,Type=String,Description="Type of variant">')
    if 'SVTYPE' not in header.info:
        header.add_line('##INFO=<ID=SVTYPE,Number=.,Type=String,Description="Type of structural variant">')
    if 'SVLEN' not in header.info:
        header.add_line('##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">')

    header.add_line('##INFO=<ID=REFLEN,Number=1,Type=Integer,Description="Length of REF allele">')
    header.add_line('##INFO=<ID=ID_LIST,Number=.,Type=String,Description="List of SVs collapsed into this SV">')

    return header


GT_DICT = {0: 0,
           1: 1,
           -1: None}


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

    # test only: 
    # for one haplotype, maximum one SV
    if has_var_arr.sum(axis=0).max() > 1:
        logger = logging.getLogger(__name__)
        logger.error(f'More than 1 SV on the same haplotype, collapse SV ID: {var_collapse.id}')
        raise ValueError(f'More than 1 SV on the same haplotype, collapse SV ID: {var_collapse.id}')
    # . and 1 should not exist on the same haplotype
    # TODO: when extending to overlapping bubbles, this may not true, need testing
    if (merge_has_var & merge_missing).any():
        logger = logging.getLogger(__name__)
        logger.warning(f'{var_collapse.id} has inconsistent genotypes (. vs 1) on the same haplotype after collapsing')
    
    # update genotypes
    update_gt_all(var_collapse, merge_gt)

    # update AC, AN, AF
    ac = int(np.sum(merge_has_var))
    an = int(np.sum(~merge_missing))
    af = float(ac / an) if an > 0 else None

    var_collapse.info['AC'] = (ac, )
    var_collapse.info['AN'] = an
    var_collapse.info['AF'] = (af, )
    
    return var_collapse


def write_outvcf(invcf: truvari.VariantFile, outvcf: truvari.VariantFile,
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
    p = truvari.VariantParams(refdist=args.refdist,
                              pctseq=args.pctseq,
                              pctsize=args.pctsize,
                              pctovl=args.pctovl,
                              no_roll=True, # In Truvari 5.3.0, no_roll=True enables rolling...
                              typeignore=False,
                              short_circuit=True,
                              skip_gt=True)
    logger.info(f'Setup Truvari Matcher: refdist={args.refdist}, pctseq={args.pctseq}, pctsize={args.pctsize}, pctovl={args.pctovl}')

    # parse input VCF
    logger.info(f'Read input VCF: {args.invcf}')
    invcf = truvari.VariantFile(args.invcf, 'rb', params=p)
    new_header = add_header(invcf.header)
    bubble_clusters = parse_vcf(invcf, args.chr, args.min_len)

    # collapse VCF   
    invcf_iter = get_vcf_iter(invcf, args.chr)
    info_lst = [] if args.info is None else args.info.split(',')

    df_collapse, df_conflict = collapse_vcf(invcf_iter, bubble_clusters, info_lst)
    logger.info(f'Collapse {len(df_collapse)} SVs into {df_collapse["Collapse_ID"].unique().shape[0]} SVs')
    logger.info(f'{len(df_conflict)} SV pairs are not collapsed due to haplotype conflict')

    # write collapse and conflict map
    file_out_collapse = f'{args.map}.collapse.txt'
    file_out_conflict = f'{args.map}.conflict.txt'

    df_collapse.to_csv(file_out_collapse, sep='\t', index=False)
    logger.info(f'Write collapse map to {file_out_collapse}')
    df_conflict.to_csv(file_out_conflict, sep='\t', index=False)
    logger.info(f'Write conflict map to {file_out_conflict}')

    # write output VCF
    outvcf = truvari.VariantFile(args.outvcf, 'w', header=new_header)

    write_outvcf(invcf, outvcf, df_collapse, args.chr, args.min_len)
    logger.info(f'Write output VCF to {args.outvcf}')

    invcf.close()
    outvcf.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='collapse_bubble.py', 
                                   description='Collapse biallelic SVs within the same bubble in VCF')

    io_arg = parser.add_argument_group('Input / Output arguments')
    io_arg.add_argument('-i', '--invcf', metavar='VCF', required=True, 
                       help='Input VCF')
    io_arg.add_argument('-o', '--outvcf', metavar='VCF', required=True, 
                       help='Output VCF')
    io_arg.add_argument('-m', '--map', metavar='PREFIX', type=str, required=True,
                       help='Write collapsed and conflicting SV tables to PREFIX.collapse.txt and PREFIX.conflict.txt.')
    io_arg.add_argument('--chr', metavar='CHR', type=str, default=None,
                       help='chromosome to work on. Default: all')
    io_arg.add_argument('--info', metavar='TAG', type=str, default=None,
                       help='Comma-separated INFO/TAG list to include in the output map. Default: None')

    collapse_arg = parser.add_argument_group('Collapse arguments')
    collapse_arg.add_argument('-l', '--min-len', metavar='50', type=int, default=50,
                             help='Minimum allele length of variants to be included, defined as max(len(alt), len(ref)). Default: 50')
    collapse_arg.add_argument('-r', '--refdist', metavar='100', type=int, default=100,
                             help='Max reference location distance. Default: 100')
    collapse_arg.add_argument('-p', '--pctseq', metavar='0.9', type=float, default=0.9,
                             help='Min percent sequence similarity (REF for DEL, ALT for other SVs). Default: 0.9')
    collapse_arg.add_argument('-P', '--pctsize', metavar='0.9', type=float, default=0.9,
                             help='Min percent size similarity (SVLEN for INS, DEL; REFLEN for INV, COMPLEX). Default: 0.9')
    collapse_arg.add_argument('-O', '--pctovl', metavar='0.9', type=float, default=0.9,
                             help='Min pct reciprocal overlap. Default: 0.9')
    
    args = parser.parse_args()
    main(args)
