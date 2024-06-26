#!/usr/bin/env python

import argparse
import csv
import sys
import re
import os

def argsParse():
    """
    Parsing input & outputs CSV files. Also takes in a boolean to indicate if
    the raw reads are single-end or paired-end
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--design", dest="design", help="Design file (csv)", default=None)
    parser.add_argument("-s", "--sampleplan", dest="sampleplan", help="SamplePlan file (csv)")
    args = parser.parse_args()
    inputDesign = args.design
    inputData = args.sampleplan
    return inputDesign, inputData

def loadSamplePlan(inputFile):
    """
    Load SamplePlan file with sampleId,sampelName,fastqR1,[fastqr2]
    """
    dictSamplePlan={'SAMPLEID':[], 'SAMPLENAME':[], 'INPUT1':[]}

    with open(inputFile, 'r') as dataFile:
        lines = csv.reader(dataFile)
        for sample in lines:
            if len(sample)>0:
                dictSamplePlan['SAMPLEID'].append(sample[0])
                dictSamplePlan['SAMPLENAME'].append(sample[1])
                dictSamplePlan['INPUT1'].append(sample[2])
                if len(sample) > 3:
                    if "INPUT2" not in dictSamplePlan:
                        dictSamplePlan['INPUT2']=[] 
                    dictSamplePlan['INPUT2'].append(sample[3])
    return(dictSamplePlan)


def loadDesign(inputFile, headers):
    """
    Load Design file using the defined headers
    """
    dictDesign = dict.fromkeys(headers, '')
    with open(inputFile, 'r') as designFile:
        lines = csv.reader(designFile)
        for sample in lines:
            if len(sample)>0:
                for i in range(len(headers)):
                    if dictDesign[headers[i]]=='':
                        dictDesign[headers[i]]=[]
                    else:
                        dictDesign[headers[i]].append(sample[i])
    return(dictDesign)


def checkHeaders(inputDesign, headerDict):
    """
    Check headers on the design file
    """
    ### Checks for design file
    with open(inputDesign, 'r') as designFile:
        lines = csv.reader(designFile)
        headers = next(lines)
        for i in range(0, len(headers)):
            try:
                if headers[i] not in headerDict: ##== [*headerDict][i]:
                    raise()
            except:
                print('\nError: Headers are not valid, should be : {}'
                      .format([*headerDict]))
                sys.exit(1)
    return headers


def checkColumnContent(column, values):
    """
    Check the content of a column
    """
    for val in column:
        if not val in values:
            print('\nError: The value \'{}\' is invalid, should be : {}'
                  .format(val, [*values]))
            sys.exit(1)

def checkColumnsMatch(col1, col2, exclusive=False):
    """
    Check that values in col1 are (not) in col2
    """
    ## Remove empty values from col1/col2
    col1 = [i for i in col1 if i]
    col2 = [i for i in col2 if i]

    match=[]
    for ID in col1:
        if ID in col2:
            if exclusive:
                print('\nError: The value {} cannot be set in two columns'
                      .format(ID))
                sys.exit(1)
            else:
                if not ID in match:
                    match.append(ID)

    if not exclusive and len(set(col1).difference(match)) != 0:
        print('\nError: Values {} are not found in {}'
              .format(set(col1).difference(match), set(col2)))
        sys.exit(1)


if __name__ == '__main__':

    designHeader=['GERMLINE_ID','TUMOR_ID','PAIR_ID','SEX']

    ## Get args
    inputDesign, inputSamplePlan = argsParse()
    
    ## Load SamplePlan
    print("[SAMPLEPLAN] Load data ", end='...')
    dictSamplePlan=loadSamplePlan(inputSamplePlan)
    print("ok") 

    ## Check Design headers
    print("[DESIGN] Check headers ", end='...')
    designHeader = checkHeaders(inputDesign, designHeader)
    print("ok")

    ## Load Design
    print("[DESIGN] Load data ", end="...")
    dictDesign=loadDesign(inputDesign, designHeader)
    print("ok")

    print(dictDesign)
    ## Checks for design file
    print("[DESIGN] Check peak type content ", end='...')
    checkColumnContent(dictDesign['SEX'], ['', 'XX', 'XY'])
    print("ok")

    ## Check that a sample is not a control
    print("[DESIGN] Check tumors/controls IDs ", end='...') 
    checkColumnsMatch(dictDesign['GERMLINE_ID'], dictDesign['TUMOR_ID'], exclusive=True)
    print("ok")

    ## Check that all samples from samplePlan are also in the design file (and the reverse)
    print("[DESIGN] Check samples matches between samplePlan and design ", end='...')
    checkColumnsMatch(dictSamplePlan['SAMPLEID'], dictDesign['GERMLINE_ID'] + dictDesign['TUMOR_ID'])
    checkColumnsMatch(dictDesign['TUMOR_ID'] + dictDesign['GERMLINE_ID'], dictSamplePlan['SAMPLEID'])
    print("ok")


