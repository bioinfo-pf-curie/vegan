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

n_header=0
#echo -e "Number_of_reads,Fragment_length,Number_of_aligned_reads,Percent_of_aligned_reads, Percent_of_overlap,Number_of_duplicates,Percent_of_duplicates,Number_reads_ontarget,Percent_reads_ontarget,Number_reads_after_filt,Percent_reads_after_filt,Mean_depth,30X_cov,50X_cov,100X_cov" > mqc.stats

for sample in $all_samples; do
  #SAMPLE NAME
  sname=$(grep "$sample," $splan | awk -F, 'NR==1{print $2}')
  #sname=$(echo $sname_raw | sed -e 's/\-/_/g' -e 's/\./_/g' -e 's/\ /_/g')
  #sample=$(echo $sample_raw | sed -e 's/\-/_/g' -e 's/\./_/g' -e 's/\ /_/g')

  header="Saple_ID,Sample_name"
  output="${sample},${sname}"

  #ALIGNMENT
  if [[ -e mapping/${sample}.flagstat ]]; then
    nb_reads=$(grep 'in total' mapping/${sample}.flagstat | awk '{print $1}')
    if [[ $is_pe == 1 ]]; then
      nb_frag=$(($nb_reads / 2))
    else
      nb_frag=$nb_reads
    fi

    #Mapping stats (always in reads - so must be converted for PE)
    if [[ $is_pe == 1 ]]; then
	nb_paired_mapped=$(grep "with itself and mate mapped" mapping/${sample}.flagstat | awk '{print $1}')
	nb_single_mapped=$(grep "singletons" mapping/${sample}.flagstat | awk '{print $1}')
	nb_mapped=$(( $nb_paired_mapped + $nb_single_mapped ))
	perc_mapped=$(echo "${nb_mapped} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    else
	nb_mapped=$(grep "primary mapped (" mapping/${sample}.flagstat | awk '{print $1}')
	perc_mapped=$(echo "${nb_mapped} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    fi
    header+=",Number_of_fragments,Number_of_aligned_reads,Percent_of_aligned_reads"
    output+=",${nb_frag},${nb_mapped},${perc_mapped}"
  fi


  ## Fragment length
  if [[ -e metrics/${sample}_insert_size_metrics.txt ]]; then
    frag_length=$(grep -A2 "## METRIC" metrics/${sample}_insert_size_metrics.txt | tail -n 1 | awk '{print $1}')
    header+=",Fragment_length"
    output+=",${frag_length}"
  fi

  ## Overlap of read pairs
  if ls metrics/${sample}*_collect_wgs_metrics.txt 1> /dev/null 2>&1; then
    ## PERC_EXC_OVERLAP is between 0 and 0.5 as this is the fraction of aligned bases that would be filtered out because they were the second observation from an insert with overlapping reads.
    ## So we multiply by 2 to have the % of base overlap
    perc_over=$(grep -A2 "## METRIC" metrics/${sample}*_collect_wgs_metrics.txt | tail -n 1 | awk '{print $11}')
    perc_over=$(echo "${perc_over}" | awk ' { printf "%.*f",2,$1*100 } ')
    header+=",Percent_of_overlap"
    output+=",${perc_over}"
  fi

  #duplicates
  if [[ -e preprocessing/${sample}.md.flagstat ]]; then
    nb_dups=$(grep 'primary duplicates' preprocessing/${sample}.md.flagstat | awk '{print $1}')
    perc_dups=$(echo "${nb_dups} ${nb_mapped}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    header+=",Number_of_duplicates,Percent_of_duplicates"
    output+=",${nb_dups},${perc_dups}"
  fi

  #On target
  if [[ -e preprocessing/${sample}.onTarget.stats ]]; then
    nb_ontarget=$(grep "reads mapped:" preprocessing/${sample}.onTarget.stats | awk '{print $4}')
    perc_ontarget=$(echo "${nb_ontarget} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    header+=",Number_reads_ontarget,Percent_reads_ontarget"
    output+=",${nb_ontarget},${perc_ontarget}"
  fi

  #Filtering
  if ls preprocessing/${sample}.filtered.flagstat 1> /dev/null; then
    nb_filt=$(grep "primary mapped (" preprocessing/${sample}*filtered.flagstat | awk '{print $1}')
    perc_filt=$(echo "${nb_filt} ${nb_reads}" | awk ' { printf "%.*f",2,$1*100/$2 } ')
    header+=",Number_reads_after_filt,Percent_reads_after_filt"
    output+=",${nb_filt},${perc_filt}"
  fi

  ## Coverage
  if [[ -e coverage/${sample}.regions.bed.gz ]]; then
    mean_depth=$(zcat coverage/${sample}.regions.bed.gz | awk '{ w=($3-$2); cov+=w*$5; tot+=w}END{print cov/tot}')
    cov30=$(zcat coverage/${sample}.regions.bed.gz | awk '{ tot+=($3-$2) } $5>=30{ z+=($3-$2)} END{printf "%.*f",2,z/tot*100 }')
    cov50=$(zcat coverage/${sample}.regions.bed.gz | awk '{ tot+=($3-$2) } $5>=50{ z+=($3-$2)} END{printf "%.*f",2,z/tot*100 }')
    cov100=$(zcat coverage/${sample}.regions.bed.gz | awk '{ tot+=($3-$2) } $5>=100{ z+=($3-$2)} END{printf "%.*f",2,z/tot*100 }')
    header+=",Mean_depth,30X_cov,50X_cov,100X_cov"
    output+=",${mean_depth},${cov30},${cov50},${cov100}"
  elif [[ -e coverage/${sample}.mosdepth.global.dist.txt && -e coverage/${sample}.mosdepth.summary.txt ]]; then
    mean_depth=$(tail -n 1 coverage/${sample}.mosdepth.summary.txt | awk '{print $4}')
    cov30=$(awk 'BEGIN{cov=0}$1=="total" && $2=="30"{cov=$3*100}END{printf "%.*f",2,cov}' coverage/${sample}.mosdepth.global.dist.txt)
    cov50=$(awk 'BEGIN{cov=0}$1=="total" && $2=="50"{cov=$3*100}END{printf "%.*f",2,cov}' coverage/${sample}.mosdepth.global.dist.txt)
    cov100=$(awk 'BEGIN{cov=0}$1=="total" && $2=="100"{cov=$3*100}END{printf "%.*f",2,cov}' coverage/${sample}.mosdepth.global.dist.txt)
    header+=",Mean_depth,30X_cov,50X_cov,100X_cov"
    output+=",${mean_depth},${cov30},${cov50},${cov100}"
  fi

  if [ $n_header == 0 ]; then
    echo -e $header > mqc.stats
    n_header=1
  fi
  echo -e $output >> mqc.stats
done
