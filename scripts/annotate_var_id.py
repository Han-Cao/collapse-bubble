#!/usr/bin/env python3

# Annotate and assign unique variant ID for pangenome VCF
# Last update: 25-Jun-2025
# Author: Han Cao

import argparse
from collections import defaultdict

import pysam

parser = argparse.ArgumentParser(prog='annotate_var_id.py', description='Annotate and assign unique variant ID for pangenome VCF')
parser.add_argument('-i', '--input', metavar='VCF', help='Input VCF', required=True)
parser.add_argument('-o', '--output', metavar='VCF', help='Output VCF', required=True)
parser.add_argument('--suffix-sep', default=None, type=str,
                    help='Separator between bubble ID and suffix, e.g., "_" for vcfwave processed VCF (default: None)')


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


def main(input: str, output: str, suffix_sep) -> None:

    invcf = pysam.VariantFile(input, 'rb')
    header = invcf.header
    header.add_line('##INFO=<ID=BUBBLE_ID,Number=1,Type=String,Description="ID of pangenome bubble">')
    outvcf = pysam.VariantFile(output, 'w', header=header)

    id_dict = defaultdict(int)

    for record in invcf:
        if suffix_sep is not None:
            bubble_id = record.id.rsplit(suffix_sep, 1)[0]
        else:
            bubble_id = record.id

        record.info['BUBBLE_ID'] = bubble_id
        var_type = get_var_type(record)

        # generate new unique ID
        # Format: <Bubble_ID>.<SVTYPE>.<No.>
        # start from 1
        id_dict[bubble_id] += 1
        new_id = bubble_id + '.' + var_type + '.' + str(id_dict[bubble_id])

        record.id = new_id
        outvcf.write(record)
    
    invcf.close()
    outvcf.close()

if __name__ == '__main__':
    args = parser.parse_args()
    main(args.input, args.output, args.suffix_sep)
