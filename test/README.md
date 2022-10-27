# Test data for Vegan pipeline

## Test datasets
- exoBench dataset: KDI project ExoBench n°1560. Tumor Sample D251E02, Normal Sample D251E01. For this dataset, two versions exists, one with only 22000 reads. This one is available as fastq, bams or vcfs. And one with 10 Million reads but only as fastq.
- hpc-bench-wgs-13x: historical test dataset from TCGA that was downsampled

## Test to perform to ensure operational pipeline

Here is the command line to run on the cluster to test the pipeline:

```
nextflow run main.nf -profile test,multiconda \
--condaCacheDir /bioinfo/local/curie/ngs-data-analysis/centos/tools/containers/conda/vegan-1.2.0/ \
--genomeAnnotationPath '/data/annotations/pipelines/' \
--skipMultiQC --skipPreseq \
-w /PATH/TO/WORKDIR/ \
--outDir /PATH/TO/RESDIR/ \
-resume
```

In order to make sure Vegan is operational, several tests must be performed. They are listed below:
- Test the different steps of the pipeline with the parameter `--step mapping` or `filtering` or `calling` or `annotate`
- test the different tools of the pipeline with the parameter `--tools mutect2,snpeff,haplotypecaller,manta,tmb,msisensor,facets,ascat`. Be aware that with small test datasets msisensor, facets and Ascat won't run. They need bigger datasets.
- test the different genome versions available with the parameter `--genome hg19` or `hg38`
- test the use of a Target file (.bed) with the parameter `--targetBed /PATH/TO/FILE`
- test the different containers types with `-profile multiconda` or `singularity`
- test Vegan with WES or WGS data
- test vegan after deployment with geniac
