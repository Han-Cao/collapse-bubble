#!/usr/bin/env python3

# Annotate and assign unique variant ID for pangenome VCF
# Last update: 30-Sep-2024
# Author: Han Cao

import argparse
from collections import defaultdict

import pysam

parser = argparse.ArgumentParser(prog='annotate_var_id.py', description='Annotate and assign unique variant ID for pangenome VCF')
parser.add_argument('-i', '--input', metavar='VCF', help='Input VCF', required=True)
parser.add_argument('-o', '--output', metavar='VCF', help='Output VCF', required=True)


def get_var_type(record: pysam.VariantRecord) -> str:
    len_ref = len(record.ref)
    len_alt = len(record.alts[0])

    if len_ref == 1:
        if len_alt > 1:
            return 'INS'
        else:
            return 'SNP'
    else:
        if len_ref == len_alt:
            return 'MNP'
        elif len_alt == 1:
            return 'DEL'
        else:
            return 'COMPLEX'


def main(input: str, output: str):

    invcf = pysam.VariantFile(input, 'rb')
    header = invcf.header
    header.add_line('##INFO=<ID=BUBBLE_ID,Number=1,Type=String,Description="ID of pangenome bubble">')
    outvcf = pysam.VariantFile(output, 'w', header=header)

    id_dict = defaultdict(int)

    for record in invcf:
        record.info['BUBBLE_ID'] = record.id
        var_type = get_var_type(record)

        # generate new unique ID
        # Format: <Bubble_ID>.<SVTYPE>.<No.>
        # start from 1
        id_dict[record.id] += 1
        new_id = record.info['BUBBLE_ID'] + '.' + var_type + '.' + str(id_dict[record.id])

        record.id = new_id
        outvcf.write(record)
    
    invcf.close()
    outvcf.close()

if __name__ == '__main__':
    args = parser.parse_args()
    main(args.input, args.output)
