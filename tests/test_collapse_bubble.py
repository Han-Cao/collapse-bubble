import pytest
import subprocess
import os
import shutil

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
TYPE = ['disjoint']


# Clean previous test run and create output directory
def prepare_outdir() -> None:

    if os.path.exists(OUTPUT_DIR):
        shutil.rmtree(OUTPUT_DIR)
    os.makedirs(OUTPUT_DIR)

# Run command
def run_script(vcf_type: str) -> None:
    command = ['python', SCRIPT, 
               '-i', os.path.join(INPUT_DIR, vcf_type + '.input.vcf.gz'), 
               '-o',  os.path.join(OUTPUT_DIR, vcf_type + '.output.vcf.gz'),
               '-m',  os.path.join(OUTPUT_DIR, vcf_type + '.output.mapping.txt'),
               '--info', 'SVTYPE,SVLEN']
    subprocess.run(command, check=True)


# Test if the command runs successfully
@pytest.mark.order(1)
@pytest.mark.parametrize("vcf_type", TYPE)
def test_script_execution(vcf_type: str) -> None:
    prepare_outdir()
    try:
        run_script(vcf_type)  # This will raise CalledProcessError on failure
        assert True
    except subprocess.CalledProcessError as e:
        pytest.fail(f"{SCRIPT_NAME} failed with error: {e}")


# Validate output VCF
def has_gt(x) -> bool:
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
        has_gt_lst.append([not has_gt(x) for sample in variant.samples.values() for x in sample['GT']])
    
    return np.array(id_lst), np.array(has_var_lst), np.array(has_gt_lst)


@pytest.mark.order(2)
@pytest.mark.parametrize("vcf_type", TYPE)
def test_validate_output(vcf_type: str) -> None:
    file_invcf = os.path.join(INPUT_DIR, vcf_type + '.input.vcf.gz')
    file_outvcf = os.path.join(OUTPUT_DIR, vcf_type + '.output.vcf.gz')
    file_mapping = os.path.join(OUTPUT_DIR, vcf_type + '.output.mapping.txt')

    invcf = pysam.VariantFile(file_invcf, 'rb')
    outvcf = pysam.VariantFile(file_outvcf, 'rb')
    df_map = pd.read_csv(file_mapping, sep='\t')

    # generate genotype matrix
    in_idx, in_has_var, in_has_gt = vcf2np(invcf)
    out_idx, out_has_var, out_has_gt = vcf2np(outvcf)

    # variant ID to np.array index mapping
    df_in_idx = pd.DataFrame(data = {'input': list(range(len(in_idx)))}, index = in_idx)
    df_out_idx = pd.DataFrame(data = {'output': list(range(len(out_idx)))}, index = out_idx)
    df_idx = pd.merge(df_in_idx, df_out_idx, how='outer', left_index=True, right_index=True)
    df_idx = df_idx.fillna(-1)
    df_idx['output'] = df_idx['output'].astype(int)

    # 1. check number of varaints
    in_has_var_sum = in_has_var.sum(axis=0)
    out_has_var_sum = out_has_var.sum(axis=0)
    assert np.array_equal(in_has_var_sum, out_has_var_sum), \
        f'Error: Number of variants ({in_has_var_sum}) in input VCF and output VCF ({out_has_var_sum}) are different'

    # 2. check genotypes for variants without collapse
    id_no_collapse = df_idx.index[~df_idx.index.isin(df_map['Variant_ID'])]
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

    collapse_lst = df_map['Collapse_ID'].unique()
    # check if variant number matched
    assert len(in_idx) == (len(out_idx) + df_map.shape[0] - len(collapse_lst)), \
        f'Error: Input variants no. ({len(in_idx)}) != output variants no. ({len(out_idx)}) + collapsed variants no. ({df_map.shape[0] - len(collapse_lst)})'

    # check genotypes
    for collapse_id in collapse_lst:
        df_collapse = df_map.loc[df_map['Collapse_ID'] == collapse_id]
        id_lst = df_collapse['Variant_ID'].to_list()
        
        # check consistentcy on missing genotypes
        assert np.array_equal(in_has_gt[df_idx.loc[id_lst, 'input']].all(axis=0),
                              out_has_gt[df_idx.loc[collapse_id, 'output']]), \
            f'Error: Collapsed variants in {id_lst} have inconsistent missing genotypes'

        # check consistency on genotypes
        assert np.array_equal(in_has_var[df_idx.loc[id_lst, 'input']].sum(axis=0),
                              out_has_var[df_idx.loc[collapse_id, 'output']]), \
            f'Error: Collapsed variants in {id_lst} have inconsistent genotypes'
        
    invcf.close()
    outvcf.close()

# Compare output files with truth
@pytest.mark.order(3)
@pytest.mark.parametrize("vcf_type", TYPE)
def test_output(vcf_type: str) -> None:

    # Get truth and test files for this suffix
    file_truth_vcf = os.path.join(TRUTH_DIR, vcf_type + '.output.vcf.gz')
    file_truth_mapping = os.path.join(TRUTH_DIR, vcf_type + '.output.mapping.txt')
    file_output_vcf = os.path.join(OUTPUT_DIR, vcf_type + '.output.vcf.gz')
    file_output_mapping = os.path.join(OUTPUT_DIR, vcf_type + '.output.mapping.txt')

    # Read files
    truth_vcf = pysam.VariantFile(file_truth_vcf, 'rb')
    df_mapping_truth = pd.read_csv(file_truth_mapping, sep='\t')
    output_vcf = pysam.VariantFile(file_output_vcf, 'rb')
    df_mapping_output = pd.read_csv(file_output_mapping, sep='\t')

    # compare mapping files
    assert df_mapping_truth.equals(df_mapping_output), f"Mapping mismatch for VCF: {vcf_type}"

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
    os.remove(file_output_mapping)