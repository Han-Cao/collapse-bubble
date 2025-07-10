import pytest
import subprocess
import os

import pysam
import pandas as pd
import numpy as np

# Constants
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_NAME = 'collapse_bubble.py'
SCRIPT = os.path.join(TEST_DIR, '..', 'scripts', SCRIPT_NAME)
INPUT_DIR = os.path.join(TEST_DIR, 'collapse_bubble', 'input')
TRUTH_DIR = os.path.join(TEST_DIR, 'collapse_bubble', 'truth')
OUTPUT_DIR = os.path.join(TEST_DIR, 'collapse_bubble', 'output')
# TODO: overlap, merge_repeat, merge_position
TYPE = ['disjoint', 'overlap', 'chain', 'merge_repeat', 'no_collapse']


# Run command
def run_script(vcf_type: str) -> None:

    invcf = os.path.join(INPUT_DIR, vcf_type + '.input.vcf.gz')
    outvcf = os.path.join(OUTPUT_DIR, vcf_type + '.output.vcf.gz')
    outmap = os.path.join(OUTPUT_DIR, vcf_type + '.output.mapping')

    # create output dir if not exist
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    # clean previous test run if exist
    if os.path.exists(outvcf):
        os.remove(outvcf)

    command = [SCRIPT, 
               '-i', invcf, 
               '-o', outvcf,
               '-m', outmap,
               '-r', '100', '-p', '0.9', '-P', '0.9', '-O', '0.9',
               '--info', 'SVTYPE,SVLEN']
    
    subprocess.run(command, check=True)


# Test if the command runs successfully
@pytest.mark.order(1)
@pytest.mark.parametrize("vcf_type", TYPE)
def test_script_execution(vcf_type: str) -> None:
    try:
        run_script(vcf_type)  # This will raise CalledProcessError on failure
        assert True
    except subprocess.CalledProcessError as e:
        pytest.fail(f"{SCRIPT_NAME} failed with error: {e}")


# Validate output VCF
def is_none(x) -> bool:
    return x is None

def has_var(x) -> bool:
    return x == 1

def vcf2np(vcf: pysam.VariantFile) -> tuple:
    
    id_lst = []
    has_var_lst = []
    has_gt_lst = []

    for variant in vcf.fetch():
        id_lst.append(variant.id)
        has_var_lst.append([has_var(x) for sample in variant.samples.values() for x in sample['GT']])
        has_gt_lst.append([not is_none(x) for sample in variant.samples.values() for x in sample['GT']])
    
    return np.array(id_lst), np.array(has_var_lst), np.array(has_gt_lst)


def is_conflict(has_var1: np.ndarray, has_var2: np.ndarray) -> bool:
    return (has_var1 & has_var2).any()


@pytest.mark.order(2)
@pytest.mark.parametrize("vcf_type", TYPE)
def test_validate_output(vcf_type: str) -> None:
    file_invcf = os.path.join(INPUT_DIR, vcf_type + '.input.vcf.gz')
    file_outvcf = os.path.join(OUTPUT_DIR, vcf_type + '.output.vcf.gz')
    file_mapping = os.path.join(OUTPUT_DIR, vcf_type + '.output.mapping.collapse.txt')
    file_conflict = os.path.join(OUTPUT_DIR, vcf_type + '.output.mapping.conflict.txt')

    invcf = pysam.VariantFile(file_invcf, 'rb')
    outvcf = pysam.VariantFile(file_outvcf, 'rb')
    df_map = pd.read_csv(file_mapping, sep='\t')
    df_conflict = pd.read_csv(file_conflict, sep='\t')

    # generate genotype matrix
    in_idx, in_has_var, in_has_gt = vcf2np(invcf)
    out_idx, out_has_var, out_has_gt = vcf2np(outvcf)

    # variant ID to np.array index mapping
    df_in_idx = pd.DataFrame(data = {'input': list(range(len(in_idx)))}, index = in_idx)
    df_out_idx = pd.DataFrame(data = {'output': list(range(len(out_idx)))}, index = out_idx)
    df_idx = pd.merge(df_in_idx, df_out_idx, how='outer', left_index=True, right_index=True)
    df_idx = df_idx.fillna(-1)
    df_idx['output'] = df_idx['output'].astype(int)

    # 0. check if mapping files meet requirements
    uniq_collapse_id = df_map['Collapse_ID'].unique()
    uniq_variants_id = df_map['Variant_ID'].unique()

    # collapse id should not exist in variant_id column
    assert df_map['Variant_ID'].isin(uniq_collapse_id).sum() == 0, \
        f'Error: Collapse ID column has overlapping with Variant ID column for {vcf_type}'
    # matching and conflict should not overlap
    df_merge_map_conflict = pd.merge(df_map, df_conflict, how='inner', on=['Variant_ID', 'Collapse_ID'])
    assert len(df_merge_map_conflict) == 0, \
        f'Error: Collapse mapping file include conflicting pairs for {vcf_type}'
    # no duplicated records in mapping files
    assert len(df_map) == len(df_map.drop_duplicates(['Variant_ID', 'Collapse_ID'])), \
        f'Error: Collapse mapping file has duplicated records for {vcf_type}'
    assert len(df_conflict) == len(df_conflict.drop_duplicates(['Variant_ID', 'Collapse_ID'])), \
        f'Error: Conflict mapping file has duplicated records for {vcf_type}'
    # no duplicated merging in mapping files
    assert len(uniq_variants_id) == len(df_map), \
        f'Error: A single variant collapse to multiple variants in {vcf_type}'

    # 1. check number of varaints
    in_has_var_sum = in_has_var.sum(axis=0)
    out_has_var_sum = out_has_var.sum(axis=0)
    assert np.array_equal(in_has_var_sum, out_has_var_sum), \
        f'Error: Number of variants in the input VCF ({in_has_var_sum}) and output VCF ({out_has_var_sum}) are different'

    # 2. check genotypes for variants without collapse
    id_no_collapse = df_idx.index[~(df_idx.index.isin(uniq_collapse_id) | df_idx.index.isin(uniq_variants_id))]
    id_missing = id_no_collapse[df_idx.loc[id_no_collapse, 'output'] == -1].tolist()

    assert len(id_missing) == 0, \
        f'Error: Variants {id_missing} without collapse are missing in output VCF'

    assert np.array_equal(in_has_var[df_idx.loc[id_no_collapse, 'input']], 
                        out_has_var[df_idx.loc[id_no_collapse, 'output']]), \
        'Error: Some variants without collapse have different genotypes'

    assert np.array_equal(in_has_gt[df_idx.loc[id_no_collapse, 'input']], 
                        out_has_gt[df_idx.loc[id_no_collapse, 'output']]), \
        'Error: Some variants without collapse have different missing genotypes'

    # 3. check collapsed variants

    # check if variant number matched
    assert len(in_idx) == (len(out_idx) + df_map.shape[0]), \
        f'Error: Input variants no. ({len(in_idx)}) != output variants no. ({len(out_idx)}) + collapsed variants no. ({df_map.shape[0]})'

    # check genotypes
    for collapse_id in uniq_collapse_id:
        df_collapse = df_map.loc[df_map['Collapse_ID'] == collapse_id]
        id_lst = df_collapse['Variant_ID'].to_list()

        # make sure the collapsed_id not in id_lst
        assert collapse_id not in id_lst

        id_lst += [collapse_id]
        # check consistentcy on missing genotypes
        assert np.array_equal(in_has_gt[df_idx.loc[id_lst, 'input']].all(axis=0),
                              out_has_gt[df_idx.loc[collapse_id, 'output']]), \
            f'Error: Collapsed variants in {id_lst} have inconsistent missing genotypes'

        # check consistency on genotypes
        assert np.array_equal(in_has_var[df_idx.loc[id_lst, 'input']].sum(axis=0),
                              out_has_var[df_idx.loc[collapse_id, 'output']]), \
            f'Error: Collapsed variants in {id_lst} have inconsistent genotypes'
    
    # 4. check conflicting haplotypes
    for row in df_conflict.itertuples():
        if row.Variant_ID in df_map['Variant_ID'].values:
            target_id = df_map.loc[df_map['Variant_ID'] == row.Variant_ID, 'Collapse_ID'].values[0]
        else:
            target_id = row.Variant_ID
        if row.Collapse_ID in df_map['Variant_ID'].values:
            conflict_id = df_map.loc[df_map['Variant_ID'] == row.Collapse_ID, 'Collapse_ID'].values[0]
        else:
            conflict_id = row.Collapse_ID
        
        target_has_var = out_has_var[df_idx.loc[target_id, 'output']]
        conflict_has_var = out_has_var[df_idx.loc[conflict_id, 'output']]
        assert is_conflict(target_has_var, conflict_has_var), \
            f'Error: {row.Variant_ID} and {row.Collapse_ID} are not conflicting'
        
    invcf.close()
    outvcf.close()


# Compare output files with truth
@pytest.mark.order(3)
@pytest.mark.parametrize("vcf_type", TYPE)
def test_output(vcf_type: str) -> None:

    # Get truth and test files for this suffix
    file_truth_vcf = os.path.join(TRUTH_DIR, vcf_type + '.output.vcf.gz')
    file_truth_collapse = os.path.join(TRUTH_DIR, vcf_type + '.output.mapping.collapse.txt')
    file_truth_conflict = os.path.join(TRUTH_DIR, vcf_type + '.output.mapping.conflict.txt')
    file_output_vcf = os.path.join(OUTPUT_DIR, vcf_type + '.output.vcf.gz')
    file_output_collapse = os.path.join(OUTPUT_DIR, vcf_type + '.output.mapping.collapse.txt')
    file_output_confilict = os.path.join(OUTPUT_DIR, vcf_type + '.output.mapping.conflict.txt')

    # Read files
    truth_vcf = pysam.VariantFile(file_truth_vcf, 'rb')
    df_collapse_truth = pd.read_csv(file_truth_collapse, sep='\t')
    df_confilict_truth = pd.read_csv(file_truth_conflict, sep='\t')
    output_vcf = pysam.VariantFile(file_output_vcf, 'rb')
    df_collapse_output = pd.read_csv(file_output_collapse, sep='\t')
    df_confilict_output = pd.read_csv(file_output_confilict, sep='\t')

    # compare mapping files
    # for multie bubble variant, we need to first sort Bubble_ID to ensure identical output
    df_collapse_truth['Bubble_ID'] = df_collapse_truth['Bubble_ID'].str.split(',').apply(sorted).str.join(',')
    df_collapse_output['Bubble_ID'] = df_collapse_output['Bubble_ID'].str.split(',').apply(sorted).str.join(',')
    assert df_collapse_truth.equals(df_collapse_output), f"Collapse mapping mismatch for VCF: {vcf_type}"
    assert df_confilict_truth.equals(df_confilict_output), f"Conflict mapping mismatch for VCF: {vcf_type}"

    # Compare header
    for truth_header, output_header in zip(truth_vcf.header.records, output_vcf.header.records):
        assert str(truth_header) == str(output_header), f"Header mismatch for VCF: {vcf_type}"
    
    # Compare samples
    assert list(truth_vcf.header.samples) == list(output_vcf.header.samples), f"Sample mismatch for VCF: {vcf_type}"

    # Compare records
    for truth_record, output_record in zip(truth_vcf, output_vcf):
        assert str(truth_record) == str(output_record), f"Record mismatch for VCF: {vcf_type}"

    truth_vcf.close()
    output_vcf.close()

    # Clean up
    os.remove(file_output_vcf)
    os.remove(file_output_collapse)
    os.remove(file_output_confilict)