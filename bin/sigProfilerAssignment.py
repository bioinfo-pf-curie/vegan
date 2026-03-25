#!/usr/bin/env python

import os
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--in_dir", help="Input directory", type=str, default='')
parser.add_argument("-o", "--out_dir", help="Output directory", type=str, default='')
parser.add_argument("-e", "--exome", help="Enable exome mode", action="store_true")
parser.add_argument("-a", "--assembly", help="Assembly (GRCh38 or GRCh37)", type=str, default='GRCh38')
parser.add_argument("-c", "--cosmic_version", help="Cosmic Version (3.3)", type=float, default=3.3)
parser.add_argument("-v", "--volume", help="TSB reference files directory", type=str, default=None)

args = parser.parse_args()

import SigProfilerAssignment as spa

from SigProfilerAssignment import Analyzer as Analyze


Analyze.cosmic_fit(samples=args.in_dir, output=args.out_dir, input_type="vcf", context_type="96",
                   collapse_to_SBS96=True, cosmic_version=args.cosmic_version, exome=args.exome,
                   genome_build=args.assembly, signature_database=None,
                   exclude_signature_subgroups=None, export_probabilities=False,
                   export_probabilities_per_mutation=False, make_plots=True,
                   sample_reconstruction_plots=False, verbose=False, volume=args.volume)

version = spa.__version__
output_text = f"SigProfilerAssignment {version}"
file_path = os.path.join(os.getcwd(), "versions.txt")
with open(file_path, 'w') as file:
    file.write(output_text) 
