import pytest
import subprocess
import os
 
import pysam

# Constants
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_NAME = 'annotate_var_id.py'
SCRIPT = os.path.join(TEST_DIR, '..', 'scripts', SCRIPT_NAME)
INPUT_DIR = os.path.join(TEST_DIR, 'annotate_var_id', 'input')
TRUTH_DIR = os.path.join(TEST_DIR, 'annotate_var_id', 'truth')
OUTPUT_DIR = os.path.join(TEST_DIR, 'annotate_var_id', 'output')
TYPE = ['raw', 'wave']

# Run command
def run_script(vcf_type: str) -> None:

    invcf = os.path.join(INPUT_DIR, vcf_type + '.input.vcf.gz')
    outvcf = os.path.join(OUTPUT_DIR, vcf_type + '.output.vcf.gz')

    # create output dir if not exist
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    # clean previous test run if exist
    if os.path.exists(outvcf):
        os.remove(outvcf)

    command = [SCRIPT, '-i', invcf, '-o',  outvcf]
    
    if vcf_type == 'wave':
        command += ['--suffix-sep', '_']

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


# Compare output files with truth
@pytest.mark.order(2)
@pytest.mark.parametrize("vcf_type", TYPE)
def test_output(vcf_type: str) -> None:

    # Get truth and test files for this suffix
    file_truth_vcf = os.path.join(TRUTH_DIR, vcf_type + '.output.vcf.gz')
    file_output_vcf = os.path.join(OUTPUT_DIR, vcf_type + '.output.vcf.gz')

    # Read files
    truth_vcf = pysam.VariantFile(file_truth_vcf, 'rb')
    output_vcf = pysam.VariantFile(file_output_vcf, 'rb')

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

    os.remove(file_output_vcf)
