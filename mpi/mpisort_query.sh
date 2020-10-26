#!/bin/bash


#SBATCH -o /mnt/beegfs/home/fjarlier/NEXTFLOW_MPI/output_query_mpisort.%j.o
#SBATCH -e /mnt/beegfs/home/fjarlier/NEXTFLOW_MPI/output_query_mpisort.%j.e

module load gnu8/8.3.0
module load openmpi3

MPISORT=$2
SAM=$1
OUTPUTDIR=$3

mpirun $MPISORT $SAM $OUTPUTDIR -p -b -n -u -q 0 
