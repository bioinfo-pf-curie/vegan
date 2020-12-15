#!/bin/bash

function usage {
    echo -e "usage : stats2multiqc.sh -s SAMPLE_PLAN -d DESIGN [-p][-h]"
    echo -e "Use option -h|--help for more information"
}

function help {
    usage;
    echo 
    echo "stat2multiqc.sh"
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -s SAMPLE_PLAN"
    echo "   -d DESIGN"
    echo "   [-p]: paired-end mode"
    echo "   [-h]: help"
    exit;
}

is_pe=0
while getopts "s:d:ph" OPT
do
    case $OPT in
        s) splan=$OPTARG;;
	d) design=$OPTARG;;
	p) is_pe=1;;
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

if  [[ -z $splan ]]; then
    usage
    exit
fi

all_samples=$(awk -F, '{print $1}' $splan)

echo -e "Sample_ID,Sample_name,Number_of_reads,Fragment_length,Number_of_aligned_reads,Percent_of_aligned_reads,Number_of_hq_mapped_reads,Percent_of_hq_mapped_reads,Number_of_lq_mapped_reads,Percent_of_lq_mapped_reads,Number_of_duplicates,Percent_of_duplicates,Number_reads_on_target,Percent_reads_on_target,Mean_depth" > mqc.stats

for sample in $all_samples
do
    #SAMPLE NAME
    sname=$(grep "$sample," $splan | awk -F, '{print $2}')

    #ALIGNMENT
    nb_reads=$(grep 'Total' Mapping/${sample}_bwa.log | awk -F "\t" '{print $2}')
    if [[ $is_pe == 1 ]]; then
	nb_frag=$(( $nb_reads / 2 ))
    else
	nb_frag=$nb_reads
    fi
    tail -n +3 Mapping/${sample}_bwa.log > Mapping/${sample}_bwa.mqc

    #Mapping stats (always in reads - so must be converted for PE)
    #These statistics are calculated after spike cleaning but before filtering
    nb_mapped=$(awk -F, '$1=="Mapped"{print $2}' Mapping/${sample}_mappingstats.mqc)
    nb_mapped_hq=$(awk -F, '$1=="HighQual"{print $2}' Mapping/${sample}_mappingstats.mqc)
    nb_mapped_lq=$(awk -F, '$1=="LowQual"{print $2}' Mapping/${sample}_mappingstats.mqc)
    perc_mapped=$(echo "${nb_mapped} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    perc_mapped_hq=$(echo "${nb_mapped_hq} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    perc_mapped_lq=$(echo "${nb_mapped_lq} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')

    #PICARD
    if [[ -e MarkDuplicates/${sample}.MarkDuplicates.metrics.txt ]]; then
	nb_dups_pair=$(grep -a2 "## METRICS" mapping/${sample}.MarkDuplicates.metrics.txt | tail -1 | awk -F"\t" '{print $7}')
	nb_dups_single=$(grep -a2 "## METRICS" mapping/${sample}.MarkDuplicates.metrics.txt | tail -1 | awk -F"\t" '{print $6}')
	nb_dups_optical=$(grep -a2 "## METRICS" mapping/${sample}.MarkDuplicates.metrics.txt | tail -1 | awk -F"\t" '{print $8}')
	nb_dups=$(( $nb_dups_pair * 2 + $nb_dups_single + $nb_dups_optical ))
	perc_dups=$(echo "${nb_dups} ${nb_mapped}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    else
	nb_dups='NA'
	perc_dups='NA'
    fi

    nb_ontarget='NA'
    perc_ontarget='NA'
    mean_depth='NA'
    frag_length='NA'

    echo -e ${sample},${sname},${nb_frag},${frag_length},${nb_mapped},${perc_mapped},${nb_mapped_hq},${perc_mapped_hq},${nb_mapped_lq},${perc_mapped_lq},${nb_dups},${perc_dups},${nb_ontarget},${perc_ontarget},${mean_depth} >> mqc.stats
done

