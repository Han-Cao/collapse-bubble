#!/usr/bin/env python3

# SV merging for pangenome VCF
# Last update: 29-Dec-2025
# Author: Han Cao
# Contributor: Quanyu Chen

import logging
import argparse
from collections import defaultdict, deque

import pysam
import truvari
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)
logger.propagate = False

class BubbleClusters:
    """ Clusters of bubbles overlapping with the same tandem repeat """

    def __init__(self) -> None:

        self.bubbles = {} # bubble_id -> {'id': [var_ids], 'svtype': [var_svtypes], 'cluster': cluster_id}
        self.bubbles_id = [] # list of bubble_ids, this keep the bubble order in the VCF
        self.cluster = defaultdict(set) # cluster_id -> set of bubble_ids
        self.cluster_vars = {} # cluster_id -> {'bubble_id': [var_ids], '_MULTI': [var_ids]}
        self.cluster_nvar = {} # cluster_id -> number of variants
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
    

    def validate_cluster_variants(self, cluster_id: int, bubble_id: str, var_lst: list[pysam.VariantRecord]) -> None:
        """ Validate if all variants in a bubble have been collected (Test only) """

        assert cluster_id in self.cluster_vars
        assert bubble_id in self.cluster_vars[cluster_id]
        check_id_lst = [x.id for x in var_lst]
        check_id_set = set(check_id_lst)
        assert len(check_id_set) == len(check_id_lst)
        assert set(self.cluster_vars[cluster_id][bubble_id]) == check_id_set
    

    def get_cluster_nvar(self, cluster_id: int) -> int:
        """ Get number of variants in a cluster """

        return self.cluster_nvar[cluster_id]
      

    def retrieve_svtype(self, variant: pysam.VariantRecord) -> str:
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
                    logger.debug(f'Overlap between {b1} and {b2} at {self._chr}:{self._pos}')
    

    def _get_all_overlap(self, b_id: str) -> set:
        """ Find all overlapping bubbles by BFS """

        visited = {b_id}
        queue = deque([b_id])
        while queue:
            current = queue.popleft()
            for neighbor in self.overlaps.get(current, set()):
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)

        return visited
    
    
    def add_variant(self, variant: pysam.VariantRecord, svtype: str) -> None:
        """ Add variant to bubble and check for overlaps """

        # get bubble id
        bubble_ids = get_bubble_ids(variant)
        
        # collect bubbles at the same position
        if variant.pos == self._pos and variant.chrom == self._chr:
            self._pos_bubbles.update(bubble_ids)

        # check if VCF is sorted
        elif variant.pos < self._pos and variant.chrom == self._chr:
            raise ValueError('Input VCF is not sorted')

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
                self.bubbles_id.append(b_id)
            else:
                # check unique id
                if variant.id in self.bubbles[b_id]['id']:
                    raise ValueError(f'Duplicate variant ID {variant.id} in bubble {b_id}')
                self.bubbles[b_id]['id'].append(variant.id)
                self.bubbles[b_id]['svtype'].append(svtype)

    
    def make_clusters(self) -> None:
        """ Cluster bubbles """

        logger.debug('Making clusters...')
        cluster_id = 0
        for b_id in self.bubbles_id:
            # skip if this bubble was already clustered
            if self.bubbles[b_id]['cluster'] is not None:
                continue
            # find overlap bubbles (only bubbles with SVs are included)
            overlap_bubbles = [x for x in self._get_all_overlap(b_id) if x in self.bubbles]
            
            # single bubble cluster
            if len(overlap_bubbles) == 0:
                logger.debug(f'Single-bubble cluster {cluster_id}: {b_id}')
                self.bubbles[b_id]['cluster'] = cluster_id
                self.cluster[cluster_id].add(b_id)
                cluster_id += 1
            
            # multiple bubble cluster
            else:
                logger.debug(f'Bubble {b_id} overlaps with {overlap_bubbles}')
                # find if any bubble already clustered
                # TODO: this should be useless when using BFS to cluster bubbles, consider removing
                exist_cluster = []
                cluster_bubbles = [b_id] + overlap_bubbles
                
                for id in cluster_bubbles:
                    if self.bubbles[id]['cluster'] is not None:
                        exist_cluster.append(self.bubbles[id]['cluster'])
                        logger.debug(f'Cluster {self.bubbles[id]["cluster"]} already exists for bubble {id}')

                # create a new cluster
                unique_cluster = list(set(exist_cluster))
                
                if len(unique_cluster) == 0:
                    logger.debug(f'Multi-bubble cluster {cluster_id}: {cluster_bubbles}')
                    for id in cluster_bubbles:
                        self.bubbles[id]['cluster'] = cluster_id
                    self.cluster[cluster_id].update(cluster_bubbles)
                    cluster_id += 1
                
                # add to existing cluster
                # TODO: this should be useless when using BFS to cluster bubbles, consider removing
                else:
                    logger.warning(f'Bubble {b_id} is clustered 2 times, this should not happen!')
                    # already clustered bubbles must in the same cluster
                    logger.debug(f'Assigning {b_id} to existing cluster {unique_cluster}')
                    assert len(unique_cluster) == 1
                    for id in cluster_bubbles:
                        self.bubbles[id]['cluster'] = unique_cluster[0]
                    self.cluster[unique_cluster[0]].update(cluster_bubbles)
        
        # test only: check all bubbles have been clustered
        for b_id in self.bubbles.keys():
            assert self.bubbles[b_id]['cluster'] is not None

        # clean cluster
        self._clean_clusters()


    def _clean_clusters(self) -> None:
        """ Remove clusters with only one SV per type, prepare variant list per cluster """

        drop_clusters = []
        # identify clusters to drop
        for c_id in self.cluster:
            svtypes = []
            var_ids = []
            for b_id in self.cluster[c_id]:
                svtypes += self.bubbles[b_id]['svtype']
                var_ids += self.bubbles[b_id]['id']

            var_ids = np.array(var_ids)
            svtypes = np.array(svtypes)
            # find concat vars by checking duplicates
            uniq_id, uniq_idx, uniq_counts = np.unique(var_ids, return_index=True, return_counts=True)
            # count svtypes after deduplicate by var_ids
            _, sv_counts = np.unique(svtypes[uniq_idx], return_counts=True)
            
            # if only one SV per type, drop
            if sv_counts.max() == 1:
                drop_clusters.append(c_id)
            # if more than one SV per type, save var_ids for collapsing
            # Note: in cluster_vars, we separate multi bubble variants
            # this allow we can first collapse variants within the same bubble, then across bubbles
            else:
                self.cluster_nvar[c_id] = len(uniq_id)
                concat_ids = uniq_id[uniq_counts > 1].tolist()
                self.cluster_vars[c_id] = {'_MULTI': concat_ids}
                for b_id in self.cluster[c_id]:
                    # exclude concatenated variants
                    self.cluster_vars[c_id][b_id] = [x for x in self.bubbles[b_id]['id'] if x not in concat_ids]
        
        # remove single SV cluster and bubbles
        for c_id in drop_clusters:
            for b_id in self.cluster[c_id]:
                del self.bubbles[b_id]
            del self.cluster[c_id]


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
    """ Get bubble ID list/ of a variant """

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

    # find clusters
    bubble_clusters.make_clusters()

    logger.info(f'Found {len(bubble_clusters.cluster)} bubble clusters with more than 1 SVs per SVTYPE')

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

    # TODO: consider only generate gt_array for collapse_var once to speed up
    # given the current sample size of pangenomes, this may not be very important
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


def match_summary(collapse_id: str, match_result: truvari.MatchResult) -> dict:
    """ Collect matching stats """

    return {
        'Collapse_ID': collapse_id,
        'PctSeqSimilarity': match_result.seqsim,
        'PctSizeSimilarity': match_result.sizesim,
        'PctRecOverlap': match_result.ovlpct,
        'SizeDiff': match_result.sizediff,
        'StartDistance': match_result.st_dist,
        'EndDistance': match_result.ed_dist,
        'TruScore': match_result.score
    }


def collapse_bubble(var_lst: list[truvari.VariantRecord], collapse_chain: dict) -> tuple[list, dict, dict, list]:
    """ 
    Collapse SVs within the same bubble or cluster, this will be run in 2 passes:
    1. collapse SVs within the same bubble
    2. collapse SVs within across bubbles within the same cluster

    Input:
    var_lst: list of SVs to be collapsed
    collapse_chain: dict of SV ID -> {Collapsed SV ID} in previous run
    
    Return:
    keep_vars: list of non-redundant SVs after collapse
    collapse_dict: dict of collapase id -> [collapsed VariantRecord]
    match_map: dict of SV matching results
    conflict_map: dict of conflicting SVs {'Variant_ID': 'Collapse_ID'}

    Note:
    For conflicting SVs, if variants within same bubble, they are compared by raw genotypes
    if variants across bubbles, they are compared after within bubble collapse
    """
    
    keep_vars = [] # list of collapsed SV
    collapse_dict = defaultdict(list) # Collapsed SV ID -> list of SV records
    match_map = {} # Variant ID -> {Collapsed SV ID, Matching stats}
    conflict_map = {} # Variant ID -> Collapsed SV ID


    # start from the most frequent SVs
    var_remain = sorted(var_lst, key=mac, reverse=True)
    while len(var_remain) > 0:
        collapse_var = var_remain[0]
        drop_idx = [0]

        # SV comparison
        for i in range(1, len(var_remain)):
            candidate_var = var_remain[i]
            res_match = collapse_var.match(candidate_var)
            logger.debug(f'Directly compare {candidate_var.id} with {collapse_var.id}: {res_match.state}. ' + 
                         f'(seqsim:{res_match.seqsim}, sizesim:{res_match.sizesim}, ovlpct:{res_match.ovlpct})')
            if res_match.state:
                # check haplotype consistency against collapsed SV
                if hap_conflict(candidate_var, collapse_var):
                    conflict_map[candidate_var.id] = collapse_var.id
                    continue
                
                # if candidate var is already collased from multiple varinats
                # check all the chained variants to avoid over-merging
                # we don't need to check haplotype conflict as they are all merged into candidate_var
                # TODO: consider add --chain like truvari?
                if candidate_var.id in collapse_chain:
                    flag_chain_break = False
                    chain_match_map = {} # cached match_lst for chained variants
                    for chain_var in collapse_chain[candidate_var.id]:
                        chain_match = collapse_var.match(chain_var)
                        logger.debug(f'Chained compare {chain_var.id} with {collapse_var.id}: {chain_match.state}. ' + 
                                     f'(seqsim:{chain_match.seqsim}, sizesim:{chain_match.sizesim}, ovlpct:{chain_match.ovlpct})')
                        # collect match summary
                        if chain_match.state:
                            chain_match_map[chain_var.id] = match_summary(collapse_var.id, chain_match)
                        else:
                            flag_chain_break = True
                            break
                    
                    if flag_chain_break:
                        continue
                    # only update match_map if this variant can be merged (i.e., no chain break)
                    else:
                        match_map.update(chain_match_map)

                # merge and update collapsed SV genotype
                collapse_var = collapse_genotype([collapse_var, candidate_var], collapse_var, update_info=False)

                # update SV merging results
                match_map[candidate_var.id] = match_summary(collapse_var.id, res_match)
                collapse_dict[collapse_var.id].append(candidate_var)

                drop_idx.append(i)
        
        # all variants have been compared, update variant list
        keep_vars.append(collapse_var)
        drop_set = set(drop_idx)
        var_remain = [x for i, x in enumerate(var_remain) if i not in drop_set]
            
    return keep_vars, collapse_dict, match_map, conflict_map


def get_variant_info(variant: pysam.VariantRecord, match_summary_dict: dict, info_lst: list) -> dict:
    """ Save variant info into dict """

    bubble_ids = get_bubble_ids(variant)
    if len(bubble_ids) == 1:
        bubble_id = bubble_ids[0]
    else:
        bubble_id = ','.join(bubble_ids)

    output = {'Chrom': variant.chrom,
              'Position': variant.pos,
              'Bubble_ID': bubble_id,
              'Variant_ID': variant.id,
            }
    output.update(match_summary_dict)

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
    
    working_clusters = {} # cluster_id -> {'bubble_id': [var_ids], '_MULTI': [var_ids]}
    working_clusters_nvar = defaultdict(int) # cluster_id -> number of variants
    collapse_map_lst = [] # list of matching results
    conflict_map = {} # Variant_ID -> Collapse_ID

    for variant in vcf:
        bubble_ids = get_bubble_ids(variant)
        cluster_id = bubble_clusters.get_bubble_cluster(bubble_ids[0])
        # bubble without cluster means it doesn't need to be collapsed
        if cluster_id is None:
            continue

        svtype = bubble_clusters.retrieve_svtype(variant)
        # not a SV or singleton SV that not need collapse
        if svtype is None:
            continue
        
        if cluster_id not in working_clusters:
            working_clusters[cluster_id] = defaultdict(list)
        # cache annotated SVTYPE and SVLEN for comparison
        bubble_id = bubble_ids[0] if len(bubble_ids) == 1 else '_MULTI'
        working_clusters[cluster_id][bubble_id].append(annotate_sv(variant, svtype))
        working_clusters_nvar[cluster_id] += 1

        # collapse if all SVs in this cluster have been processed
        n_working = working_clusters_nvar[cluster_id]
        n_total = bubble_clusters.get_cluster_nvar(cluster_id)

        if n_working == n_total:
            cluster_vars = [] # variants for pass 2
            cluster_match = {} # variant_id -> match_summary
            collapse_chain = {} # chained variants for pass 2
            # pass 1: within bubble collapse
            for b_id, b_vars in working_clusters[cluster_id].items():
                bubble_clusters.validate_cluster_variants(cluster_id, b_id, b_vars) # test only

                if b_id == '_MULTI':
                    continue
                
                b_remain_vars, b_collapse, b_match, b_conflict = collapse_bubble(b_vars, collapse_chain={})
                # prepare for pass 2
                assert set(b_collapse.keys()).isdisjoint(set(collapse_chain.keys())) # test only
                cluster_vars += b_remain_vars
                collapse_chain.update(b_collapse)
                cluster_match.update(b_match)
                conflict_map.update(b_conflict)
            
            # include multi-bubble variants
            if '_MULTI' in working_clusters[cluster_id]:
                cluster_vars += working_clusters[cluster_id]['_MULTI']
            
            # pass 2: within cluster collapse
            _, _, c_match, c_conflict = collapse_bubble(cluster_vars, collapse_chain)

            cluster_match.update(c_match)
            conflict_map.update(c_conflict)

            # extract variant annotations and prepare for output
            for _, b_vars in working_clusters[cluster_id].items():
                for var in b_vars:
                    if var.id in cluster_match:
                        collapse_map_lst.append(get_variant_info(var, cluster_match[var.id], info_lst))

            # free memory
            del working_clusters[cluster_id]
    
    # If everthing works well, all clusters should be removed from working_clusters
    assert len(working_clusters) == 0

    # convert to dataframe
    df_conflict = pd.DataFrame(list(conflict_map.items()), columns=['Variant_ID', 'Collapse_ID'])
    df_collapse = pd.DataFrame.from_dict(collapse_map_lst)

    # no merging
    if df_collapse.empty:
        return df_collapse, df_conflict

    # check Collapse_ID no longer collapse to any other variant
    assert df_collapse['Collapse_ID'].isin(df_collapse['Variant_ID']).sum() == 0
    # remove duplicated conflict, this is due to 

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


def collapse_genotype(var_lst: list, var_collapse: pysam.VariantRecord, update_info: bool=True) -> pysam.VariantRecord:
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
        logger.error(f'More than 1 SV on the same haplotype, collapse SV ID: {var_collapse.id}')
        raise ValueError(f'More than 1 SV on the same haplotype, collapse SV ID: {var_collapse.id}')
    # . and 1 should not exist on the same haplotype
    # TODO: when extending to overlapping bubbles, this may not true, need testing
    if (merge_has_var & merge_missing).any():
        logger.warning(f'{var_collapse.id} has inconsistent genotypes (. vs 1) on the same haplotype after collapsing')
    
    # update genotypes
    update_gt_all(var_collapse, merge_gt)

    # update AC, AN, AF
    if update_info:
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

    # handle with empty collapse df
    if df_collapse.empty:
        id_map = {}
    else:
        # extract ID mapping dict from df
        df_id_map = df_collapse[['Variant_ID', 'Collapse_ID']]
        uniq_collapse_ids = df_id_map['Collapse_ID'].unique()
        # include self mapping, i.e., Collapse_ID -> Collapse_ID
        df_id_map = pd.concat([df_id_map, pd.DataFrame({'Variant_ID': uniq_collapse_ids, 
                                                        'Collapse_ID': uniq_collapse_ids})])
        id_map = df_id_map.set_index('Variant_ID')['Collapse_ID'].to_dict()

        # get summary of collapsed SV
        df_summary = df_id_map.groupby('Collapse_ID').agg({'Variant_ID': 'count'})
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


def setup_logger(log_level) -> None:
    """ Setup logger """
    
    # Clear and setup
    logger.handlers.clear()
    logger.setLevel(log_level)
    
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        '[%(asctime)s] - [%(levelname)s]: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    
    return None


def main() -> None:

    args = parse_args()

    # setup logger
    if args.debug:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO
    setup_logger(log_level)

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
    if df_collapse.empty:
        logger.warning('No SVs are collapsed, please check input VCF')
    else:
        logger.info(f'Collapse {len(df_collapse)} SVs into {df_collapse["Collapse_ID"].unique().shape[0]} SVs')
    logger.info(f'{len(df_conflict)} SV pairs are not collapsed due to haplotype conflict')
    # write output VCF
    outvcf = truvari.VariantFile(args.outvcf, 'w', header=new_header)

    write_outvcf(invcf, outvcf, df_collapse, args.chr, args.min_len)
    logger.info(f'Write output VCF to {args.outvcf}')

    # write collapse and conflict map
    collapse_cols = ['Chrom', 'Position', 'Bubble_ID', 'Variant_ID', 'Collapse_ID']
    collapse_cols += [f'INFO_{x}' for x in info_lst]
    collapse_cols += ['PctSeqSimilarity', 'PctSizeSimilarity', 'PctRecOverlap', 
                      'SizeDiff', 'StartDistance', 'EndDistance', 'TruScore']
    df_collapse = df_collapse[collapse_cols] if not df_collapse.empty else pd.DataFrame(columns=collapse_cols)
    file_out_collapse = f'{args.map}.collapse.txt'
    file_out_conflict = f'{args.map}.conflict.txt'

    df_collapse.to_csv(file_out_collapse, sep='\t', index=False)
    logger.info(f'Write collapse map to {file_out_collapse}')
    df_conflict.to_csv(file_out_conflict, sep='\t', index=False)
    logger.info(f'Write conflict map to {file_out_conflict}')

    invcf.close()
    outvcf.close()


def parse_args() -> argparse.Namespace:
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
                        help='chromosome to work on, all if not specified. Default: %(default)s')
    io_arg.add_argument('--info', metavar='TAG', type=str, default=None,
                        help='Comma-separated INFO/TAG list to include in the output map. Default: %(default)s')

    collapse_arg = parser.add_argument_group('Collapse arguments')
    collapse_arg.add_argument('-l', '--min-len', metavar='50', type=int, default=50,
                              help='Minimum allele length of variants to be included, defined as max(len(alt), len(ref)). Default: %(default)s')
    collapse_arg.add_argument('-r', '--refdist', metavar='500', type=int, default=500,
                              help='Max reference location distance. Default: %(default)s')
    collapse_arg.add_argument('-p', '--pctseq', metavar='0.9', type=float, default=0.9,
                              help='Min percent sequence similarity (REF for DEL, ALT for other SVs). Default: %(default)s')
    collapse_arg.add_argument('-P', '--pctsize', metavar='0.9', type=float, default=0.9,
                              help='Min percent size similarity (SVLEN for INS, DEL; REFLEN for INV, COMPLEX). Default: %(default)s')
    collapse_arg.add_argument('-O', '--pctovl', metavar='0', type=float, default=0,
                              help='Min pct reciprocal overlap. Default: %(default)s')
    
    other_arg = parser.add_argument_group('Other arguments')
    other_arg.add_argument('--debug', action='store_true', 
                           help='Debug mode')
    
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    main()
