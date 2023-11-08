# VEGAN
**V**ariant calling pipeline for whole **E**xome and whole **G**enome sequencing c**AN**cer data

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)
[![Install with conda](https://img.shields.io/badge/install%20with-conda-brightgreen.svg)](https://conda.anaconda.org/anaconda)
[![Singularity Container available](https://img.shields.io/badge/singularity-available-7E4C74.svg)](https://singularity.lbl.gov/)
[![Docker Container available](https://img.shields.io/badge/docker-available-003399.svg)](https://www.docker.com/)

### Introduction

This pipeline was built for **Whole Exome Sequencing** and **Whole Genome Sequencing** analysis. It provides a detailed quality controls of both frozen and FFPE samples as well as a first downstream analysis including mutation calling, structural variants and copy number analysis. Most of the pipeline steps can work for tumor/normal paired samples and tumor-only samples.  

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner.
It comes with conda / singularity containers making installation easier and results highly reproducible. The current workflow was inspired from the [nf-core Sarek pipeline](https://github.com/nf-core/sarek) with several common processes, and further modifications and new analysis steps.

### Pipeline summary

1. Run quality control of raw sequencing reads ([`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Align reads on reference genome ([`bwa-mem`](http://bio-bwa.sourceforge.net/), [`bwa-mem2`](https://github.com/bwa-mem2/bwa-mem2), [`dragmap`](https://github.com/Illumina/DRAGMAP))
3. Filtering and quality controls of aligned reads
  - Report mapping metrics ([`picard`](https://gatk.broadinstitute.org/hc/en-us/articles/360037055772-CollectInsertSizeMetrics-Picard-))
  - Mark and remove duplicates ([`sambamba`](https://lomereiter.github.io/sambamba/))
  - Library complexity analysis ([`Preseq`](http://smithlabresearch.org/software/preseq/))
  - Filtering aligned BAM files ([`SAMTools`](http://www.htslib.org/))
  - Insert size distribution ([`picard`](https://gatk.broadinstitute.org/hc/en-us/articles/360037055772-CollectInsertSizeMetrics-Picard-))
  -  Identity monitoring  ([`bcftools`](http://samtools.github.io/bcftools/bcftools.html) / [`R`](https://www.r-project.org/))
4. GATK preprocessing ([`GATK`](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-))
5. Germline Variants calling ([`haplotypecaller`](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) / [`bcftools`](http://samtools.github.io/bcftools/bcftools.html))
  - HaplotypeCaller
6. Somatic Variants calling ([`mutect2`](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2) / [`bcftools`](http://samtools.github.io/bcftools/bcftools.html))
  - Mutect2 (including learnReadOrientationModel, GetPileupSummaries, CalculateContamination)
  - FilterMutectCalls
7. Technical filters for somatic variants (DP, VAF, MAF) ([`SnpSift`](https://github.com/pcingola/SnpSift), [`bcftools`](http://samtools.github.io/bcftools/bcftools.html))
8. Variants annotation ([`SnpEff`](https://pcingola.github.io/SnpEff/) / [`SnpSift`](https://github.com/pcingola/SnpSift))
9. Copy-number analysis ([`ASCAT`](https://www.crick.ac.uk/research/labs/peter-van-loo/software), [`FACETS`](https://github.com/mskcc/facets))
10. Structural variants analysis ([`MANTA`](https://github.com/Illumina/manta))
11. Biomarkers analysis
  - Microsatellite instability analysis ([`MSIsensor-pro`](https://github.com/xjtu-omics/msisensor-pro))
  - Tumor Mutational Burden ([`pyTMB`](https://github.com/bioinfo-pf-curie/TMB))
12. Gather all QC results in a final report ([`MultiQC`](http://multiqc.info/))

### Quick help

```bash
nextflow run main.nf --help
N E X T F L O W  ~  version 21.10.6
Launching `main.nf` [lethal_torricelli] - revision: 4d570988d2
------------------------------------------------------------------------

    _   _   _____          __     __  _____    ____      _      _   _
   | \ | | |  ___|         \ \   / / | ____|  / ___|    / \    | \ | |
   |  \| | | |_     _____   \ \ / /  |  _|   | |  _    / _ \   |  \| |
   | |\  | |  _|   |_____|   \ V /   | |___  | |_| |  / ___ \  | |\  |
   |_| \_| |_|                \_/    |_____|  \____| /_/   \_\ |_| \_|


                   VEGAN v2.2.0
------------------------------------------------------------------------

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf --profile STRING --samplePlan PATH --design PATH --step STRING --genome STRING --genomeAnnotationPath PATH

MANDATORY ARGUMENTS:
    --design                   PATH                                                                              Path to designf ile specifying the metadata ssociated with the samples
    --genome                   STRING [hg19, hg19_base, hg38, hg38_base, mm10, mm39,...]                         Name of the reference genome.
    --genomeAnnotationPath PATH                                                                                  PATH to the reference genome folder.
    --profile                  STRING [test, multiconda, singularity, cluster, docker, conda, path, multipath]   Configuration profile to use. Can use multiple (comma separated).
    --step                     STRING [mapping, filtering, calling, annotate]                                    Specify starting step
    --outDir                   PATH                                                                              The output directory where the results will be saved
    --tools                    STRING [haplotypecaller, mutect2, manta, snpeff, facets, ascat, tmb, msisensor]   Specify tools to use for variant calling

INPUTS:
    --reads                    PATH      Path to input data (must be surrounded with quotes)
    --samplePlan               PATH      Path to sample plan (csv format) raw reads (if `--reads` is not secified), or intermediate files according to the `--step` parameter
    --singleEnd                          For single-end input data
    --splitFastq                         Split fastq files in chunks
	--fastqChunksSize          INTEGER   Reads chunks size

ALIGNMENT:
    --aligner                  STRING [bwa-mem, bwa-mem2, dragmap]   Specify tools to use for mapping
	--mapQual                  INTEGER                               Minimum mapping quality to consider for an alignment
	--saveAlignedIntermediates                                       Save intermediates alignment files
	--splitFastq                                                     Split fastq files in chunks

FILTERING:
    --keepDups                Specify to keep duplicate reads when filtering the alignment
	--keepMultiHits           Specify to keep multi hit reads when filtering the alignment
	--keepSingleton           Specify to keep singleton reads when filtering the alignment
	--targetBed     PATH      Target Bed file for targeted or whole exome sequencing

VARIANT CALLING:
    --baseQual                   INTEGER   Minimum base quality used by Facets for CNV calling
    --saveVcfIntermediates                 Save intermediate vcf files
    --saveVcfMetrics                       Save complementary vcf metrics files
    --skipMutectContamination              Do not apply the Contamination step for Mutect2 calls filtering
    --skipMutectOrientationModel           Do not apply the LearnOrientationModel step for Mutect2 calls filtering
					
TUMOR ONLY:
    --msiBaselineConfig PATH   PATH to Msisensor-pro baseline config file for tumor-only mode
    --pon               PATH   PATH to panels of normals (.vcf.gz)
    --ponIndex          PATH   PATH to panels of normals index file (.tbi)

VCF FILTERS:
    --filterSomaticDP  INTEGER   Minimum sequencing depth to consider a somatic variant
    --filterSomaticMAF INTEGER   Maximum variant frequency in the general population to consider a somatic variant
    --filterSomaticVAF INTEGER   Minimum variant allele frequency to consider a somatic variant

ANNOTATION:
    --annotDb STRING [cosmic, icgc, cancerhotspots, gnomad, dbnsfp]   Annotation databases to use with SnpEff and SnpSift
    --ffpe                                                            Specify to use the ffpe parameters and filters for TMB computation

SKIP OPTIONS:
    --skipBQSR                 Disable BQSR
    --skipBamQC                Disable QCs on BAM files
    --skipFastqc               Disable Fastqc
    --skipIdentito             Disable Identito
    --skipMultiqc              Disable MultiQC
    --skipSaturation           Disable Preseq				

OTHER OPTIONS:
    --disableAutoClean           Disable cleaning of work directory
    --multiqcConfig    PATH      Specify a custom config file for MultiQC
    --name             STRING    Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
    --sequencingCenter STRING    Name of sequencing center to be displayed in BAM file

=======================================================
Available Profiles
   -profile test                        Run the test dataset
   -profile conda                       Build a new conda environment before running the pipeline. Use `--condaCacheDir` to define the conda cache path
   -profile multiconda                  Build a new conda environment per process before running the pipeline. Use `--condaCacheDir` to define the conda cache path
   -profile path                        Use the installation path defined for all tools. Use `--globalPath` to define the insallation path
   -profile multipath                   Use the installation paths defined for each tool. Use `--globalPath` to define the insallation path
   -profile docker                      Use the Docker images for each process
   -profile singularity                 Use the Singularity images for each process. Use `--singularityPath` to define the insallation path
   -profile cluster                     Run the workflow on the cluster, instead of locally
```
### Quick run

The pipeline can be run on any infrastructure from a list of input files or from a sample plan as follow

#### Run the pipeline on the test dataset

The test dataset is a downsampled Whole Exome Sequencing.
It can be launched with the following command.

```
nextflow run main.nf -profile test,multiconda \
   --step mapping \ # or filtering, calling, annotate
   --condaCacheDir /bioinfo/local/curie/ngs-data-analysis/centos/tools/containers/conda/vegan-2.0.0/ \
   --genomeAnnotationPath /data/annotations/pipelines/
```

#### Run the pipeline for WES analysis from a sample plan with specified tools and genome on the cluster, using singularity containers

```
nextflow run main.nf -profile singularity,cluster \
    --samplePlan samples-WES.csv \
    --design samples.design.csv \
    --step mapping \
    --singularityImagePath /bioinfo/local/curie/ngs-data-analysis/centos/tools/containers/singularity/vegan-2.0.0/images/ \
    --targetBed capture.bed \
    --tools manta,mutect2,snpeff,facets,tmb,haplotypecaller,msisensor \
    --genome hg19 --genomeAnnotationPath /data/annotations/pipelines/ \
    -resume
```


#### Run the pipeline on the cluster, using existing conda

```
nextflow run main.nf -profile multiconda,cluster \
    --samplePlan samples-WES.csv \
    --design samples.design.csv \
    --step mapping \
    --targetBed capture.bed \
    --tools manta,mutect2,snpeff,facets,tmb,haplotypecaller,msisensor,ascat \
    --genome hg19_base \
    --genomeAnnotationPath /data/annotations/pipelines/ \
    --condaCacheDir /bioinfo/local/curie/ngs-data-analysis/centos/tools/containers/conda/vegan-1.2.0/ \

```

To build new conda environments, point to an empty folder for `--condaCacheDir` parameter

#### Defining the '-profile'
By default (whithout any profile), Nextflow will excute the pipeline locally, expecting that all tools are available from your `PATH` variable.
In addition, we set up a few profiles that should allow you
- 1) to use containers instead of local installation,
- 2) to run the pipeline on a cluster instead of on a local architecture.
The description of each profile is available on the help message (see above).

Here are a few examples of how to set the profile option. See the [full documentation](docs/profiles) for details.

```
## Run the pipeline locally, using a global environment where all tools are installed (build by conda for instance)
-profile path --globalPath INSTALLATION_PATH

## Run the pipeline on the cluster, using the Singularity containers
-profile cluster,singularity --singularityImagePath SINGULARITY_PATH

## Run the pipeline on the cluster, building new conda environments
-profile cluster,multiconda --condaCacheDir CONDA_CACHE

```
#### Sample Plan

A sample plan is a **csv** file (comma separated) that list all samples with their biological IDs, **with no header**.  
The sample plan is expected to be created as below :

SAMPLE_ID | SAMPLE_NAME | PATH_TO_R1_FASTQ | [PATH_TO_R2_FASTQ]

#### Design

A design file is a **csv** file that list all experimental samples, their IDs, the associated germinal sample, the sex of the patient and the status (tumor / normal). The design control is expected to have **the following header** :

GERMLINE_ID | TUMOR_ID | PAIR_ID | SEX

Both files will be checked by the pipeline and have to be rigorously defined in order to make the pipeline work.
Note that the control is optional if not available but is highly recommanded.
If the design file is not specified, the pipeline will run until the alignment. The variant calling and the annotation will be skipped.

### Full Documentation

1. [Installation](docs/installation.md)
2. [Geniac](docs/geniac.md)
3. [Reference genomes](docs/referenceGenomes.md)
4. [Running the pipeline](docs/usage.md)
5. [Profiles](docs/profiles.md)
6. [Output and how to interpret the results](docs/output.md)
7. [Troubleshooting](docs/troubleshooting.md)

#### Fundings

This pipeline has been written by the Institut Curie bioinformatics platform (T. Gutman, F. Allain, , P. La Rosa, P. Hupe, N. Servant). The project was funded by the European Union’s Horizon 2020 research and innovation programme and the Canadian Institutes of Health Research under the grant agreement No 825835 in the framework of the [European-Canadian Cancer Network](https://eucancan.com/), as well as the Canceropole Ile de France (GENOPROFILE - RIC2021) project.

#### Contacts

For any question, bug or suggestion, please send an issue or contact the bioinformatics core facility.
