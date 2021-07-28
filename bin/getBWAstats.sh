#!/bin/bash 

## BWA mapping statistics
## Note that the statistics are in number of reads (not pairs)

function usage {
    echo -e "usage : getBWAstats.sh -i INPUT [-p PROC] [-h]"
    echo -e "Use option -h|--help for more information"
}

function help {
    usage;
    echo 
    echo "getBWAstats.sh"
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "   -i INPUT"
    echo "   [-p] PROC"
    echo "   [-h]: help"
    exit;
}

proc=1
while getopts "i:p:h" OPT
do
    case $OPT in
        i) input=$OPTARG;;
	p) proc=$OPTARG;;
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

##########################################################
nb_pairs=$(samtools view -@ $proc -f 0x1 -c $input)
tot=$(samtools view -@ $proc -F 0x100 -F 0x800 -c $input)
mapped=$(samtools view -@ $proc -F 0x4 -F 0x100 -F 0x800 -c $input)
unmapped=$(( $tot - $mapped))
uniq=$(samtools view -@ $proc -q 1 -F 0x4 -F 0x100 -F 0x800 $input | grep -v -e 'XA:Z:' -e 'SA:Z:' | wc -l )
multi=$(( $mapped - $uniq ))
# determine number of reads with MapQ20
mapq20=$(samtools view -q 20 -@ $proc -c $input)
# determine number of duplicates
markdup=$(samtools view -c -@ $proc -f 1024 $input)

if [[ ${nb_pairs} -gt 0 ]]; then
    paired=$(samtools view -@ $proc -F 0x100 -F 0x800 -F 0x004 -F 0x0008 -f 0x001 -c $input)
    single=$(( $mapped - $paired ))
    uniq_paired=$(samtools view -@ $proc -q 1 -F 0x4 -F 0x100 -F 0x800 -F 0x0008 -f 0x1 $input | grep -v -e 'XA:Z:' -e 'SA:Z:' | wc -l)
    uniq_single=$(( $uniq - $uniq_paired ))
    multi_paired=$(( $paired - $uniq_paired ))
    multi_single=$(( $single - $uniq_single ))
    unmapped=$(samtools view -@ $proc -f 12 -c $input)

    # Back to fragment for PE
    uniq_paired=$(( $uniq_paired / 2 ))
    multi_paired=$(( $multi_paired / 2 ))
    unmapped=$(( $unmapped / 2 ))

    echo -e "Total\t${tot}" 
    echo -e "Mapped\t${mapped}" 
    echo -e "PE mapped uniquely\t${uniq_paired}" 
    echo -e "PE one mate mapped uniquely\t${uniq_single}" 
    echo -e "PE multi mapped\t${multi_paired}" 
    echo -e "PE one mate multi\t${multi_single}"
    echo -e "PE neither mate aligned\t${unmapped}"
    echo -e "Number of reads at mapq20\t${mapq20}"
    echo -e "Number of duplicated reads\t${markdup}"
else
    echo -e "Total\t${tot}" 
    echo -e "Mapped\t${mapped}" 
    echo -e "Uniquely mapped reads\t${uniq}" 
    echo -e "Multi mapped reads\t${multi}" 
    echo -e "Unmapped\t${unmapped}"
    echo -e "Number of reads at mapq20\t${mapq20}"
    echo -e "Number of duplicated reads\t${markdup}"
fi
