#!/bin/bash

function usage() {
  echo -e "usage : stats2multiqc.sh -s SAMPLE_PLAN -d DESIGN [-p][-t][-h]"
  echo -e "Use option -h|--help for more information"
}

function help() {
  usage
  echo
  echo "stat2multiqc.sh"
  echo "---------------"
  echo "OPTIONS"
  echo
  echo "   -s SAMPLE_PLAN"
  echo "   -d DESIGN"
  echo "   [-p] paired-end mode"
  echo "   [-h]: help"
  exit
}

is_pe=0
while getopts "s:d:ph" OPT; do
  case $OPT in
  s) splan=$OPTARG ;;
  d) design=$OPTARG ;;
  p) is_pe=1 ;;
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

if [[ -z $splan ]]; then
  usage
  exit
fi

all_samples=$(awk -F, '{print $1}' $splan)

echo -e "Sample_ID,Sample_name,Number_of_reads,Fragment_length,Number_of_aligned_reads,Percent_of_aligned_reads,Percent_of_overlap,Number_of_duplicates,Percent_of_duplicates,Number_reads_ontarget,Percent_reads_ontarget,Number_reads_after_filt,Percent_reads_after_filt,Mean_depth,30X_cov,50X_cov,100X_cov" > mqc.stats

for sample in $all_samples; do
  #SAMPLE NAME
  sname=$(grep "$sample," $splan | awk -F, '{print $2}')

  #ALIGNMENT
  if [[ -e mapping/${sample}_bwa.log && -e mapping/${sample}.flagstats ]]; then
    nb_reads=$(grep 'Total' mapping/${sample}_bwa.log | awk -F "\t" '{print $2}')
    if [[ $is_pe == 1 ]]; then
      nb_frag=$(($nb_reads / 2))
    else
      nb_frag=$nb_reads
    fi
    tail -n +3 mapping/${sample}_bwa.log > mapping/${sample}_bwa.mqc

    #Mapping stats (always in reads - so must be converted for PE)
    if [[ $is_pe == 1 ]]; then
	nb_paired_mapped=$(grep "with itself and mate mapped" mapping/${sample}.flagstats | awk '{print $1}')
	nb_single_mapped=$(grep "singletons" mapping/${sample}.flagstats | awk '{print $1}')
	nb_mapped=$(( $nb_paired_mapped + $nb_single_mapped ))
	perc_mapped=$(echo "${nb_mapped} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    else
	nb_mapped=$(grep "mapped (" mapping/${sample}.flagstats | awk '{print $1}')
	perc_mapped=$(echo "${nb_mapped} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    fi
  fi


  #duplicates
  if [[ -e MarkDuplicates/${sample}.md.flagstats ]]; then
    nb_dups=$(grep duplicates MarkDuplicates/${sample}.md.flagstats | awk '{print $1}')
    perc_dups=$(echo "${nb_dups} ${nb_mapped}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
  else
    nb_dups='NA'
    perc_dups='NA'
  fi

  #On target
  if [[ -e preprocessing/${sample}.onTarget.flagstats ]]; then
    nb_ontarget=$(grep "mapped (" preprocessing/${sample}.onTarget.flagstats | awk '{print $1}')
    perc_ontarget=$(echo "${nb_ontarget} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
  else
    nb_ontarget='NA'
    perc_ontarget='NA'
  fi

  #On target + filtering
  if ls preprocessing/${sample}*filtered.flagstats 1> /dev/null; then
    nb_filt=$(grep "mapped (" preprocessing/${sample}*filtered.flagstats | awk '{print $1}')
    perc_filt=$(echo "${nb_filt} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
  else
    nb_filt='NA'
    perc_filt='NA'
  fi


  if [[ -e coverage/${sample}.mosdepth.summary.txt ]]; then
    mean_depth=$(tail -n 1 coverage/${sample}.mosdepth.summary.txt | awk '{print $4}')
  else
    mean_depth='NA'
  fi

  if [[ -e coverage/${sample}.mosdepth.region.dist.txt ]]; then
    cov30=$(awk 'BEGIN{cov=0}$1=="total" && $2=="30"{cov=$3*100}END{printf "%.*f",2,cov}' coverage/${sample}.mosdepth.region.dist.txt)
    cov50=$(awk 'BEGIN{cov=0}$1=="total" && $2=="50"{cov=$3*100}END{printf "%.*f",2,cov}' coverage/${sample}.mosdepth.region.dist.txt)
    cov100=$(awk 'BEGIN{cov=0}$1=="total" && $2=="100"{cov=$3*100}END{printf "%.*f",2,cov}' coverage/${sample}.mosdepth.region.dist.txt)
  elif [[ -e coverage/${sample}.mosdepth.global.dist.txt ]]; then
    cov30=$(awk 'BEGIN{cov=0}$1=="total" && $2=="30"{cov=$3*100}END{printf "%.*f",2,cov}' coverage/${sample}.mosdepth.global.dist.txt)
    cov50=$(awk 'BEGIN{cov=0}$1=="total" && $2=="50"{cov=$3*100}END{printf "%.*f",2,cov}' coverage/${sample}.mosdepth.global.dist.txt)
    cov100=$(awk 'BEGIN{cov=0}$1=="total" && $2=="100"{cov=$3*100}END{printf "%.*f",2,cov}' coverage/${sample}.mosdepth.global.dist.txt)
  else
    cov30='NA'
    cov50='NA'
    cov100='NA'
  fi

  if [[ -e metrics/${sample}_insert_size_metrics.txt ]]; then
    frag_length=$(grep -A2 "## METRIC" metrics/${sample}_insert_size_metrics.txt | tail -n 1 | awk '{print $1}')
  else
    frag_length='NA'
  fi

  if ls metrics/${sample}*_collect_wgs_metrics.txt 1> /dev/null 2>&1; then
    ## PERC_EXC_OVERLAP is between 0 and 0.5 as this is the fraction of aligned bases that would be filtered out because they were the second observation from an insert with overlapping reads.
    ## So we multiply by 2 to have the % of base overlap
    perc_over=$(grep -A2 "## METRIC" metrics/${sample}*_collect_wgs_metrics.txt | tail -n 1 | awk '{print $11}')
    perc_over=$(echo "${perc_over}" | awk ' { printf "%.*f",2,$1*100 } ')
  else
    perc_over='NA'
  fi

  echo -e ${sample},${sname},${nb_frag},${frag_length},${nb_mapped},${perc_mapped},${perc_over},${nb_dups},${perc_dups},${nb_ontarget},${perc_ontarget},${nb_filt},${perc_filt},${mean_depth},${cov30},${cov50},${cov100} >>mqc.stats
done
