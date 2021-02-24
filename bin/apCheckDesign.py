#!/usr/bin/env python

#############################################################################################
# Copyright Institut Curie 2020                                                             #
#                                                                                           #
# This software is a computer program whose purpose                                         #
# is to analyze high-throughput sequencing data.                                            #
# You can use, modify and/ or redistribute the software under                               #
# the terms of license (see the LICENSE file for more details).                             #
# The software is distributed in the hope that it will be useful,                           #
# but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.                                             #
# Users are therefore encouraged to test the software's suitabilityas regards               #
# their requirements in conditions enabling the security of their systems and/or data.      #
# The fact that you are presently reading this means that                                   #
# you have had knowledge of the license and that you accept its terms.                      #
#############################################################################################

import argparse
import csv
import sys
import re
import os


def parse_args():
    """
    Parsing input & outputs CSV files. Also takes in a boolean to indicate if
    the raw reads are single-end or paired-end
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--design", dest="design", help="Design file (csv)", default=None)
    parser.add_argument("-s", "--sampleplan", dest="sampleplan", help="SamplePlan file (csv)")
    parser.add_argument("--singleEnd", help="Specify that input reads are single-end", action="store_true")
    args = parser.parse_args()
    return args.design, args.sampleplan, args.singleEnd


def load_sample_plan(input_file, is_single_end=False):
    """
    Load SamplePlan file with sampleId,sampelName,fastqR1,[fastqr2]
    """
    sample_plan = {'SAMPLEID': [], 'SAMPLENAME': [], 'FASTQR1': []}
    if not is_single_end:
        sample_plan['FASTQR2'] = []

    with open(input_file, 'r') as dataFile:
        lines = csv.reader(dataFile)
        for sample in lines:
            sample_plan['SAMPLEID'].append(sample[0])
            sample_plan['SAMPLENAME'].append(sample[1])
            sample_plan['FASTQR1'].append(sample[2])
            if not is_single_end:
                sample_plan['FASTQR2'].append(sample[3])
    return sample_plan


def load_design(input_file, headers):
    """
    Load Design file using the defined headers
    """
    dict_design = dict.fromkeys(headers, '')
    with open(input_file, 'r') as designFile:
        lines = csv.reader(designFile)
        for sample in lines:
            for i in range(len(headers)):
                if dict_design[headers[i]] == '':
                    dict_design[headers[i]] = []
                else:
                    dict_design[headers[i]].append(sample[i])
    return dict_design


def check_headers(input_design, header_dict):
    """
    Check headers on the design file
    """
    # Checks for design file
    with open(input_design, 'r') as design_file:
        lines = csv.reader(design_file)
        header = next(lines)
        for i in range(0, len(header)):
            try:
                if not header[i] == [*header_dict][i]:
                    raise ()
            except:
                print('\nError: Headers are not valid, should be : {}'
                      .format([*header_dict]))
                sys.exit(1)


def check_column_content(column, values):
    """
    Check the content of a column
    """
    for val in column:
        if not val in values:
            print('\nError: The value \'{}\' is invalid, should be : {}'
                  .format(val, [*values]))
            sys.exit(1)


def check_columns_match(col1, col2, exclusive=False):
    """
    Check that values in col1 are (not) in col2
    """
    # Remove empty values from col1/col2
    col1 = [i for i in col1 if i]
    col2 = [i for i in col2 if i]

    match = []
    for ID in col1:
        if ID in col2:
            if exclusive:
                print(f"\nError: The value {ID} cannot be set in two columns")
                sys.exit(1)
            else:
                if not ID in match:
                    match.append(ID)

    if not exclusive and len(set(col1).difference(match)) != 0:
        print('\nError: Values {} are not found in {}'
              .format(set(col1).difference(match), set(col2)))
        sys.exit(1)


if __name__ == '__main__':
    # define your design header
    designHeader = ['GERMLINE_ID', 'TUMOR_ID', 'PAIR_ID', 'SEX']

    # Get args
    input_design, input_sample_plan, is_single_end = parse_args()

    # Load SamplePlan
    print("[SAMPLEPLAN] Load data ", end='...')
    sample_plan = load_sample_plan(input_sample_plan, is_single_end)
    print("ok")

    # Check Design headers
    print("[DESIGN] Check headers ", end='...')
    check_headers(input_design, designHeader)
    print("ok")

    # Load Design
    print("[DESIGN] Load data ", end="...")
    design = load_design(input_design, designHeader)
    print("ok")

    # Checks for design file
    # print("[DESIGN] Check peak type content ", end='...')
    # checkColumnContent(dictDesign['PEAKTYPE'], ['sharp', 'broad', 'very-broad'])
    # print("ok")

    # Check that a sample is not a control
    # print("[DESIGN] Check samples/controls IDs ", end='...')
    # checkColumnsMatch(dictDesign['SAMPLEID'], dictDesign['CONTROLID'], exclusive=True)
    # print("ok")

    # Check that all samples from samplePlan are also in the design file (and the reverse)
    # print("[DESIGN] Check samples matches between samplePlan and design ", end='...')
    # checkColumnsMatch(dictSamplePlan['SAMPLEID'], dictDesign['SAMPLEID'] + dictDesign['CONTROLID'])
    # checkColumnsMatch(dictDesign['SAMPLEID'] + dictDesign['CONTROLID'], dictSamplePlan['SAMPLEID'])
    # print("ok")
