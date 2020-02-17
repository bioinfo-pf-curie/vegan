echo "export PATH=/bioinfo/local/build/singularity/singularity-3.5.2/bin:$PATH; cd ${PWD}; nextflow run main.nf -resume -profile singularity,cluster --step mapping --input /data/tmp/plarosa/data/Sarek-data/HPC-bench/tsv/HPC-bench-test-WGS.13X.tsv --genomeAnnotationPath /bioinfo/local/curie/ngs-data-analysis/centos/nf-sarek/annotations/GRCh38 --genome GRCh38 --singularityImagePath /data/u900pf-bioinfo/containers/commun/dev/singularity/sarek-2.5.2/images --tools strelka --targetBED /bioinfo/local/curie/ngs-data-analysis/centos/nf-sarek/data/HPC-bench/TargetRegions.bed --split_fastq 20000 -c nextflow.config --queue dev" | qsub -N SAREK-2.5.2_WGS -q dev
