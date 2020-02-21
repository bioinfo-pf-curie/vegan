echo "export PATH=/bioinfo/local/build/singularity/singularity-3.5.2/bin:$PATH; cd ${PWD}; nextflow run main.nf -resume -profile singularity,cluster --maxMemory '32.GB' --intervals /data/tmp/plarosa/data/Sarek-data/HPC-bench/wgs_calling_small_regions.hg38.bed --step mapping --input /data/tmp/plarosa/data/Sarek-data/HPC-bench/tsv/HPC-bench-test-WGS.13X.tsv --genomeAnnotationPath /bioinfo/local/curie/ngs-data-analysis/centos/nf-sarek/annotations/GRCh38 --genome GRCh38 --singularityImagePath /data/u900pf-bioinfo/containers/commun/dev/singularity/nf-wgswes/images --tools manta --skipQC bamQC -c nextflow.config --queue dev" | qsub -N NF-WGSWES -q dev
