# VEGAN  
**V**ariant calling pipeline for whole **E**xome and whole **G**enome sequencing c**AN**cer data

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)
[![Install with](https://anaconda.org/anaconda/conda-build/badges/installer/conda.svg)](https://conda.anaconda.org/anaconda)
[![Singularity Container available](https://img.shields.io/badge/singularity-available-7E4C74.svg)](https://singularity.lbl.gov/)
[![Docker Container available](https://img.shields.io/badge/docker-available-003399.svg)](https://www.docker.com/)

### Introduction

This pipeline was built for **Whole Exome Sequencing** and **Whole Genome Sequencing** analysis. It provides a detailed quality controls of both frozen and FFPE samples as well as a first downstream analysis including mutation calling, structural variants and copy number analysis.  

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner.
It comes with conda / singularity containers making installation easier and results highly reproducible. The current workflow was inspired from the [nf-core Sarek pipeline](https://github.com/nf-core/sarek).

### Pipeline summary

1. Run quality control of raw sequencing reads ([`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Align reads on reference genome ([`BWA`](http://bio-bwa.sourceforge.net/))
3. Report mapping metrics ([`picard`](https://gatk.broadinstitute.org/hc/en-us/articles/360037055772-CollectInsertSizeMetrics-Picard-))
4. Mark duplicates ([`sambamba`](https://lomereiter.github.io/sambamba/))
5. Library complexity analysis ([`Preseq`](http://smithlabresearch.org/software/preseq/))
6. Filtering aligned BAM files ([`SAMTools`](http://www.htslib.org/))
7. Calculate insert size distribution ([`picard`](https://gatk.broadinstitute.org/hc/en-us/articles/360037055772-CollectInsertSizeMetrics-Picard-))
8. Calculate genes and genome coverage ([`mosdepth`](https://github.com/brentp/mosdepth))
9. Identity monitoring and samples similarity ([`bcftools`](http://samtools.github.io/bcftools/bcftools.html) / [`R`](https://www.r-project.org/))
10. GATK preprocessing ([`GATK`](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-))
11. Germline Variants calling ([`haplotypecaller`](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller))
  - HaplotypeCaller
  - GenotypeGVCFs
12. Somatic Variants calling ([`mutect2`](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2))
  - Mutect2
  - MergeMutectStats
  - GetPileupSummaries
  - GatherPileupSummaries
  - CalculateContamination
  - FilterMutectCalls
13. Variants annotation ([`snpeff`](https://pcingola.github.io/SnpEff/))
14. Copy-number analysis ([`ASCAT`](https://www.crick.ac.uk/research/labs/peter-van-loo/software), [`FACETS`](https://github.com/mskcc/facets))
15. Structural variants analysis ([`MANTA`](https://github.com/Illumina/manta))
16. Present all QC results in a final report ([`MultiQC`](http://multiqc.info/))

### Quick help

```bash
>>nextflow ~/GitLab/vegan/main.nf --help
N E X T F L O W  ~  version 20.10.0
Launching `/data/users/nservant/GitLab/vegan/main.nf` [ridiculous_borg] - revision: a44a457623
------------------------------------------------------------------------

    _   _   _____          __     __  _____    ____      _      _   _
   | \ | | |  ___|         \ \   / / | ____|  / ___|    / \    | \ | |
   |  \| | | |_     _____   \ \ / /  |  _|   | |  _    / _ \   |  \| |
   | |\  | |  _|   |_____|   \ V /   | |___  | |_| |  / ___ \  | |\  |
   |_| \_| |_|                \_/    |_____|  \____| /_/   \_\ |_| \_|


                   VEGAN v0.1.1dev
------------------------------------------------------------------------
========================================================================
Usage:

The typical command for running the pipeline is as follows:

nextflow run main.nf --reads/--samplePlan PATH --design PATH --profile STRING --genome STRING

MANDATORY ARGUMENTS:
    --design     PATH                                                                              Path to design file specifying the metadata associated with the samples
    --genome     STRING [hg19, hg38]                                                               Name of the reference genome.
    --profile    STRING [conda, cluster, docker, multiconda, conda, path, multipath, singularity]  Configuration profile to use. Can use multiple (comma separated).
    --reads      PATH                                                                              Path to input data (must be surrounded with quotes)
    --samplePlan PATH                                                                              Path to sample plan (csv format) with raw reads (if `--reads` is not specified), or intermediate 
                                                                                                   files according to the `--step` parameter

MAIN OPTIONS:
    --noIntervals                                                                   Disable usage of intervals
    --step        STRING [mapping, recalibrate, variantcalling, annotate]           Specify starting step
    --tools       STRING [facets, ascat, haplotypecaller, manta, mutect2, snpeff]   Specify tools to use for variant calling
    --singleEnd                                                                     For single-end input data

ALIGNMENT:
    --bwaOptions               STRING    Define BWA-mem option for read mapping
    --saveAlignedIntermediates           Save intermediates alignment files
    --targetBED                PATH      Target BED file for targeted or whole exome sequencing

FILTERING:
    --SNVFilters STRING [mapq, duplicates, singleton, multihits]   Specify which filter(s) to use for SNV
    --SVFilters  STRING [mapq, duplicates, singleton, multihits]   Specify which filter(s) to use for SV
    --mapQual    INTEGER                                           Minimum mapping quality

VARIANT CALLING:
    --baseQual INTEGER   Minimum base quality used by Facets for CNV calling
    --noGVCF             No g.vcf output from HaplotypeCaller

SKIP OPTIONS:
    --skipBQSR               Disable BQSR
    --skipIdentito           Disable Identito
    --skipMultiqc            Disable MultiQC
    --skipPreseq             Disable Preseq
    --skipQC                 Specify which QC tools to skip

REFERENCES:
    --acLoci                PATH     acLoci file
    --acLociGC              PATH     acLociGC file
    --bwaIndex              PATH     BWA indexes for reads alignment
    --dbsnp                 PATH     DBSNP file
    --dbsnpIndex            PATH     DBSNP index file
    --dict                  PATH     Dict from the fasta reference
    --fasta                 PATH     Fasta reference
    --fastafai              PATH     Fasta reference index
    --germlineResource      PATH     Germline Resource File
    --germlineResourceIndex PATH     Germline Resource Index File
    --intervals             PATH     List of genomic intervals to split the processing
    --knownIndels           PATH     knownIndels File
    --knownIndelsIndex      PATH     knownIndels Index File
    --pon                   PATH     Panel-of-normals VCF (bgzipped, indexed)
    --ponIndex              PATH     Index of pon panel-of-normals VCF
    --snpEffCache           PATH     Specity the path to snpEff cache
    --snpeffDb              STRING   snpeffDb version. Initialized from genomes configuration by default

OTHER OPTIONS:
    --multiqcConfig        PATH      Specify a custom config file for MultiQC
    --name                 STRING    Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
    --nucleotidesPerSecond DECIMAL   To estimate interval size
    --outDir               PATH      The output directory where the results will be saved
    --sequencingCenter     STRING    Name of sequencing center to be displayed in BAM file
======================================================================

    Available Profiles
      -profile test                     Run the test dataset
      -profile conda                    Build a new conda environment before running the pipeline. Use `--condaCacheDir` to define the conda cache path
      -profile multiconda               Build a new conda environment per process before running the pipeline. Use `--condaCacheDir` to define the conda cache path
      -profile path                     Use the installation path defined for all tools. Use `--globalPath` to define the insallation path
      -profile multipath                Use the installation paths defined for each tool. Use `--globalPath` to define the insallation path
      -profile docker                   Use the Docker images for each process
      -profile singularity              Use the Singularity images for each process. Use `--singularityPath` to define the insallation path
      -profile cluster                  Run the workflow on the cluster, instead of locally

```

### Quick run

The pipeline can be run on any infrastructure from a list of input files or from a sample plan as follow

#### Run the pipeline on the test dataset
See the conf/test.conf to set your test dataset.

```
nextflow run main.nf -profile test,multiconda --genomeAnnotationPath ANNOTATION_PATH

```

#### Run the pipeline for WES analysis using singularity containers

```
nextflow run main.nf -profile singularity,cluster --input samples-WES.tsv --targetBED capture.bed --tools mutect2,snpeff --genome hg19 -resume
```

#### Run the pipeline from a sample plan with specified tools and genome on the cluster, using the Singularity containers

```
nextflow run main.nf -profile singularity,cluster --input samples-WGS.tsv --tools manta,mutect2,haplotypecaller,ascat,snpeff --genome hg19 -resume

```

#### Run the pipeline on the cluster, building a new conda environment

```
nextflow run main.nf -profile multiconda,cluster --input samples-WGS.tsv --tools manta,mutect2,haplotypecaller,ascat,snpeff --genome hg38 -resume

```

#### Defining the '-profile'
By default (whithout any profile), Nextflow will excute the pipeline locally, expecting that all tools are available from your `PATH` variable.
In addition, we set up a few profiles that should allow you i/ to use containers instead of local installation, ii/ to run the pipeline on a cluster instead of on a local architecture.
The description of each profile is available on the help message (see above).
Here are a few examples of how to set the profile option.

```
## Run the pipeline locally, using a global environment where all tools are installed (build by conda for instance)
-profile path --globalPath INSTALLATION_PATH

## Run the pipeline on the cluster, using the Singularity containers
-profile cluster,singularity --singularityPath SINGULARITY_PATH

## Run the pipeline on the cluster, building a new conda environment
-profile cluster,conda --condaCacheDir CONDA_CACHE

```
#### Sample Plan

A sample plan is a csv file (comma separated) that list all samples with their biological IDs. The sample plan is expected to be created as below :


Sample ID | Sample Name | Path R1 .fastq file | [Path R2 .fastq file]

#### Design

A design file is a csv file that list all experimental samples, their IDs, the associated germinal sample, the sex of the patient and the status (tumor / normal). The design control is expected to be created as below :

SAMPLE_ID | GERMLINE ID | SAMPLE_NAME | SEX | STATUS

Both files will be checked by the pipeline and have to be rigorously defined in order to make the pipeline work.
Note that the control is optional if not available but is highly recommanded.
If the design file is not specified, the pipeline will run until the alignment. The variant calling and the annotation will be skipped.


### Full Documentation

1. [Installation](docs/installation.md)
2. [Reference genomes](docs/reference_genomes.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

#### Credits

This pipeline has been written by the Institut Curie bioinformatics platform (F. Allain, T. Gutman, P. La Rosa, P. Hupe, N. Servant).

#### Contacts

For any question, bug or suggestion, please send an issue or contact the bioinformatics core facility.
