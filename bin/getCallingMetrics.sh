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
	c) conta=$OPTARG;;
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

nbVar=$(zcat ${vcf} | grep -v "^#" | wc -l)
nbFilt=$(zcat ${fvcf} | grep -v "^#" | grep PASS | wc -l)
pFilt=$(echo "${nbVar} ${nbFilt}" | awk ' { printf "%.*f",2,($1-$2)/$1*100 } ')
conta=$(awk 'NR==2{printf "%.*f",2,$2*100}' $conta)

echo "sample_ID,conta,nb_var,nb_var_filt,perc_var_filt"
echo $sname,$conta,$nbVar,$nbFilt,$pFilt
