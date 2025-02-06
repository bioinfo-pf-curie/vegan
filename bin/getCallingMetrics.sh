#!/bin/bash

# getCallingMetrics.sh -i ${unfiltered} \
#                       -f ${sampleIdTN}_${variantCaller}_filtered.vcf.gz \
#                       -c ${sampleIdTN}_contamination.table
#                       -n ${sampleIdTN} > ${sampleIdTN}_Mutect2_callingMetrics.mqc

while getopts "i:f:c:n:" OPT
do
    case $OPT in
        i) vcf=$OPTARG;;
        f) fvcf=$OPTARG;;
        n) sname=$OPTARG;;
        h) help ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            exit 1
            ;;
    esac
done

if [ -z ${vcf} ];
then
    echo "Error - VCF file is missing"
    exit 1
fi

if [ -z ${sname} ];
then
    echo "Error - Sample Name is missing"
    exit 1
fi

## All variant callers
nbVar=$(zcat ${vcf} | grep -v "^#" | wc -l)
header="sample_ID,nb_var"
output="$sname,$nbVar"

## Mutect2
if [ ! -z ${fvcf} ]; then
    nbFilt=$(zcat ${fvcf} | grep -v "^#" | wc -l)
    if [[ $nbFilt -gt 0 ]]; then
	pFilt=$(echo "${nbVar} ${nbFilt}" | awk ' { printf "%.*f",2,($1-$2)/$1*100 } ')
    else
	pFilt=0
    fi
     nvar01=$(bcftools view -Oz -i 'FORMAT/AF[1:0]>=0.01' --threads 1 ${fvcf} | gunzip -c | grep -v "^#" | wc -l)
     nvar05=$(bcftools view -Oz -i 'FORMAT/AF[1:0]>=0.05' --threads 1 ${fvcf} | gunzip -c | grep -v "^#" | wc -l)
     nvar1=$(bcftools view -Oz -i 'FORMAT/AF[1:0]>=0.1' --threads 1 ${fvcf} | gunzip -c | grep -v "^#" | wc -l)
     
     header+=",nb_var_filt,perc_var_filt,nb_var_01,nb_var_05,nb_var_1"
     output+=",$nbFilt,$pFilt,$nvar01,$nvar05,$nvar1"
fi

echo -e $header
echo -e $output

