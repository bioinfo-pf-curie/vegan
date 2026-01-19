#!/usr/bin/env python

import cyvcf2
import pandas as pd
import argparse
import sys

"""
Parse inputs
"""


def argsParse():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--vcf", help="Input file (.vcf)")
    parser.add_argument("--sample", help="Specify the sample ID to focus on", type=str, default=None)
    parser.add_argument("-o", "--output", help="Output file (.txt)")

    # Others
    parser.add_argument("--verbose", help="Active verbose mode", action="store_true")

    args = parser.parse_args()
    return (args)

if __name__ == "__main__":

    args = argsParse()

    # Load Data
    if args.sample is not None:
        vcf = cyvcf2.VCF(args.vcf)
        count = 0
        for sample in vcf.samples:
            count = count + 1
            if str(sample) == str(args.sample):
                vcf = cyvcf2.VCF(args.vcf, samples=args.sample)
            elif count == len(vcf.samples):
                print("Error: Name of the sample incorrect")
                sys.exit(-1)
    else:
        vcf = cyvcf2.VCF(args.vcf)

    # Sample name
    if len(vcf.samples) > 1:
        sys.stderr.write("Error: " + len(vcf.samples) + " sample detected. This version is designed for a single sample !")
        sys.exit(-1)


    with open (args.output, "w") as fileout:
        header="ID\tchr\tpos\tref\talt\ttype\n"
        fileout.write("{}".format(header))
        for variant in vcf:
            type = variant.var_type
            chr=variant.CHROM
            position=variant.POS
            ref=variant.REF
            alt=variant.ALT[0]
            fileout.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(vcf.samples[0],chr,position,ref,alt,type))
