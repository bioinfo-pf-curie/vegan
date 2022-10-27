# Usage

## Table of contents

* [General information](#general-information)
* [Use case](#use-case)
  * [Whole genome sequencing analysis](#whole-genome-sequencing-analysis)
  * [Whole exome sequencing analysis](#whole-exome-sequencing-analysis)
  * [Starting from intermediates results](#starting-from-intermediates-results)
* [Nextflow profiles](#nextflow-profiles)
* [Running the pipeline](#running-the-pipeline)
  * [Mandatory arguments](#mandatory-arguments)
    * [`--samplePlan`](#--samplePlan)
    * [`--design`](#--design)
    * [`--profile`](docs/profiles.md)
    * [`--step`](#--step)
    * [`--tools`](#--tools)
    * [`--outDir`](#--outDir)
  * [Inputs](#inputs)
    * [`--reads`](#--reads)
    * [`--singleEnd`](#--singleend)
  * [Alignment](#alignment)
    * [`--aligner`](#--aligner)
    * [`--mapQual`](#--mapQual)
    * [`--bwaOpts`](#--bwaOpts)
    * [`--saveAlignedIntermediates`](#--saveAlignedIntermediates)
  * [Filtering](#filtering)
    * [`--targetBed`](#--targetBed)
    * [`--keepDups`](#--kepDups)
    * [`--keepMultiHits`](#--keepMultiHits)
    * [`--keepSingleton`](#--keepSingleton)
  * [Variant calling](#variant_calling)
    * [`--baseQual`](#--baseQual)
    * [`--orientationBiais`](#--orientationBiais)
    * [`--saveVcfIntermediates`](#--saveVcfIntermediates)
    * [`--saveVcfMetrics`](#--saveVcfMetrics)
  * [Annotation](#Annotation)
    * [`--annotDb`](#--annotDb)
    * [`--ffpe`](#--ffpe)
  * [Reference genomes](#reference-genomes)
    * [`--genome`](#--genome)
    * [`--genomeAnnotationPath`](#--genomeAnnotationPath)
  * [Other command line parameters](#other-command-line-parameters)
    * [`--skip*`](#--skip*)
    * [`-w`](#-w)
    * [`-name`](#-name-single-dash)
    * [`-resume`](#-resume-single-dash)
    * [`-c`](#-c-single-dash)
    * [`--maxMemory`](#--maxMemory)
    * [`--maxTime`](#--maxTime)
    * [`--maxCpus`](#--maxCpus)
    * [`--multiqcConfig`](#--multiqcConfig)
* [Job resources](#job-resources)
* [Automatic resubmission](#automatic-resubmission)


## General inforomation

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted to your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Use case

### Whole genome sequencing analysis

Here is a typical command line to analyse whole-genome analysis data.

```bash
nextflow run main.nf -profile cluster,multiconda \
  --samplePlan splan.csv \
  --design design.csv \
  --step mapping \
  --tools ascat,manta,mutect2,snpeff \
  --genome h38 \
  --genomeAnnotationPath /PATH/TO/ANNOTATIONS \
  --condaCacheDir /PATH/TO/CONDA \
  --outDir /PATH/TO/RESULTS \
  -w /PATH/TO/WORK \
  -resume

```

### Whole exome sequencing analysis

Here is a typical command line to analyse whole-exome analysis data, thus focusing the analysis on a targetBed file.

```bash
nextflow run main.nf -profile cluster,multiconda \
  --samplePlan splan.csv \
  --design design.csv \
  --targetBed target.bed \
  --step mapping \
  --tools facets,mutect2,snpeff \
  --genome h38 \
  --genomeAnnotationPath /PATH/TO/ANNOTATIONS \
  --condaCacheDir /PATH/TO/CONDA \
  --outDir /PATH/TO/RESULTS \
  -w /PATH/TO/WORK \
  -resume
```

### Starting from intermediates results

`VEGAN` offers the possibility to start the pipeline from intermediate files.
In this case, the sample plan must be adapted.

#### Starting at Bam filtering step

In order to start at the filtering step, so just after reads mapping, you can use the option `--step 'filtering'` as follow :

```bash
nextflow run main.nf --samplePlan [RESULTS/resume/samplePlan.filtered.csv] --design [DESIGN] \
                     --step 'filtering' --genome 'hg38' --tools 'mutect2,ascat' \
                     -profile cluster,singularity --outDir [RESULTS_2] -w [RESULTS_2/work]
```

The sample plan contains the following information :

SAMPLE_ID | SAMPLE_NAME | PATH_TO_BAM_FILE | PATH_TO_BAM_INDEX

#### Starting at variant calling step

In order to start at the variant calling step, so just after BAM recalibration, you can use the option `--step 'calling'` as follow :

```bash
nextflow run main.nf --samplePlan [RESULTS/resume/samplePlan.recal.csv] --design [DESIGN] \
                     --step 'calling' --genome 'hg38' --tools 'mutect2,ascat' \
                     -profile cluster,singularity --outDir [RESULTS_3] -w [RESULTS_3/work]
```

The sample plan contains the following information :

SAMPLE_ID | SAMPLE_NAME | PATH_TO_BAM_FILE | PATH_TO_BAM_INDEX

#### Starting at the annotation step

In order to start at the annotation step, so just after variant calling, you can use the option `--step 'annotate'` as follow :

```bash
nextflow run main.nf --samplePlan [RESULTS/resume/samplePlan.vcf.csv] --design [DESIGN] \
                     --step 'annotate' --genome 'hg38' --tools 'mutect2,ascat' \
                     -profile cluster,singularity --outDir [RESULTS_3] -w [RESULTS_3/work]
```

The sample plan contains the following information :

SAMPLE_ID | SAMPLE_NAME | PATH_TO_VCF_FILE | PATH_TO_VCF_TABIX_INDEX

## Nextflow profiles

Different Nextflow profiles can be used. See [Profiles](profiles.md) for details.

## Running the pipeline

The typical command for running the pipeline is as follows:
```bash
nextflow run main.nf -profile cluster,multiconda \
  --samplePlan splan.csv \
  --design design.csv \
  --targetBed target.bed \
  --step mapping \
  --tools facets,mutect2,snpeff \
  --genome h38 \
  --genomeAnnotationPath /PATH/TO/ANNOTATIONS \
  --condaCacheDir /PATH/TO/CONDA \
  --outDir /PATH/TO/RESULTS \
  -w /PATH/TO/WORK \
  -resume
```

This will launch the pipeline with the `multiconda` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

You can change the output director and the work directory using respectively `--outDir` and `-w` options.

### Mandatory arguments

#### `--samplePlan`

Use this to specify a sample plan file instead of a regular expression to find fastq files. This is the preferred solution for Vegan pipeline. For example :

```bash
--samplePlan 'path/to/data/sample_plan.csv'
```

The sample plan is a csv file with the following information (no header):

SAMPLE_ID | SAMPLE_NAME | PATH_TO_R1_FASTQ_FILE | PATH_TO_R2_FASTQ_FILE

#### `--design`

Use this to specify a design file to list all experimental samples, their IDs, the associated germinal sample, the sex of the patient and the status (tumor / normal). For example :

```bash
--design 'path/to/design.csv'
```

The design file is a csv file with the following header :

GERMLINE_ID | TUMOR_ID | PAIR_ID | SEX

This file is mandatory in order to run the variant calling and annotation parts of the pipeline

#### `--profile`

see documentation about [profiles](docs/profiles.md)

#### `--step`

This parameter specify the starting step of the pipeline. Several entry point are available. Each step requires a specific sample_plan.

**Available:** mapping, filtering, calling, annotate.

#### `--tools`

Specify the tools to use for variant calling and downstream steps.

**Available:** facets, ascat, haplotypecaller, mutect2, manta, snpeff, msisensor, tmb

#### `--outDir`

The output directory where the results will be saved.

### `Inputs`

#### `--reads`

Use this to specify the location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*_{1,2}.fastq'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

#### `--singleEnd`

By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--singleEnd` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example:

```bash
--singleEnd --reads '*.fastq.gz'
```

It is not possible to run a mixture of single-end and paired-end files in one run.

### `Alignment`

#### `--aligner`

Change the aligner tool used for mapping step.

**Available:** bwa-mem (default), bwa-mem2 and dragmap   

#### `--mapQual`

Minimum mapping quality to consider for an alignment

**Default:** 20

#### `--bwaOpts`

Change default BWA-mem options ("-k 19 -T 30 -M") for reads alignment.

```bash
--bwaOpts "[NEW OPTIONS]"
```

#### `--saveAlignedIntermediates`

By default, only alignment files in BAM format after filtering are saved in the `--outDir` directory.  
Using this option, all intermediate BAM files are saved, including BWA-mem outputs, sambambam outputs, intersection with targets (for WES).
Note that activating this option usually consumes a large amount of disk space.

```bash
--saveAlignedIntermediates
```

### `Filtering`

#### `--targetBed`

Specify a target BED file for targeted or whole exome sequencing

```bash
--targetBed
```

#### `--keepDups`

Specify to keep duplicate reads when filtering the alignment

```bash
--keepDups
```

#### `--keepMultiHits`

Specify to keep multi hit reads when filtering the alignment

```bash
--keepMultiHits
```

#### `--keepSingleton`

Specify to keep singleton reads when filtering the alignment

```bash
--keepSingleton
```

### Variants calling

#### `--baseQual`

Define the base quality for `Facets` copy number analyis.

**Default:** 13

```bash
--baseQual 13
```

#### `--orientationBiais`

Specify to use the orientation biais filter for Mutect2 calls in order to reduce the number of false positives calls emerging from artifactual alterations in FFPE samples.

**Default:** true

```bash
--orientationBiais
```

#### `--saveVcfIntermediates`

Save intermediate vcf files for Manta, Mutect2 and HaplotypeCaller.

#### `--saveVcfMetrics`

Save complementary vcf metrics files for HaplotypeCaller and Mutect2.

### Annotation

#### `--annotDb`

Specify the database to use for the annotation with SnpsSift.

**Available:** cosmic, icgc, cancerhotspots, gnomad, dbnsfp

**Default:** all

#### `ffpe`

Specify to use the ffpe parameters and filters for TMB computation. This means a Variant Allele Frequency (VAF) of 10%.

Here the TMB is computed for WES samples (min coverage of 20X) and not for panels.

For more information about TMB compution: [documentation](https://github.com/bioinfo-pf-curie/TMB)

### Reference genomes

The pipeline config files come bundled with paths to the genomes reference files.

#### `--genome`

There are different species supported in the genomes references file. To run the pipeline, you must specify which to use with the `--genome` flag.

You can find the keys to specify the genomes in the [genomes config file](../conf/genomes.config). Common genomes that are supported are:

* Human
  * `--genome hg19`
  * `--genome hg38`
  * `--genome hg19_base` and `--genome hg38_base` (contains mitochondrial chromosome but extra contigs are removed)

* Mouse
  * `--genome mm10`

Note that you can use the same configuration setup to save sets of reference files for your own use, even if they are not part of the genomes resource.
See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions on where to save such a file.

The syntax for this reference configuration is as follows:

```nextflow
params {
  genomes {
    'hg19' {
      bwaIndex              = "/path/to/bwaIndex"
      chrLength             = "/path/to/chrom_hg19.sizes"
      dict                  = "/path/to/hg19.dict"
      fasta                 = "/path/to/hg19.fa"
      fastaFai              = "/path/to/hg19.fa.fai"
      gtf                   = "/path/to/gencode.v19.annotation_proteinCoding.gtf"
      dbsnp                 = "/path/to/dbsnp_138.hg19.vcf.gz"
      dbsnpIndex            = "/path/to/dbsnp_138.hg19.vcf.gz.tbi"
      acLoci                = "/path/to/1000G_phase3_20130502_SNP_maf0.3.loci"
      acLociGC              = "/path/to/1000G_phase3_20130502_SNP_maf0.3.loci.gc"
      polyms                = "/path/to/44polyms.bed"
      germlineResource      = "/path/to/af-only-gnomad_modified.raw.sites.vcf.gz"
      germlineResourceIndex = "/path/to/af-only-gnomad_modified.raw.sites.vcf.gz.tbi"
      intervals             = "/path/to/wgs_calling_regions.grch37.list.txt"
      knownIndels           = "/path/to/{1000G_phase1,Mills_and_1000G_gold_standard}.indels.hg19.sites.vcf.gz"
      knownIndelsIndex      = "/path/to/{1000G_phase1,Mills_and_1000G_gold_standard}.indels.hg19.sites.vcf.gz.tbi"
      snpeffDb              = "hg19"
      snpeffCache           = "/path/to/snpEff_v4_3"
    }
}
```

Note that all these paths can be updated on the command line using for example the following parameters:
- `--bwaIndex` - Path to Bwa index
- `gtf` - Path to GTF file
- ...

#### `--genomeAnnotationPath`

This parameter points to the folder containing all the annotations.

```bash
--genomeAnnotationPath /PATH/TO/ANNOTATIONS
```

### Other command line parameters

#### `--skip*`

The pipeline is made with a few *skip* options that allow to skip optional steps in the workflow.
The following options can be used:

- `--skipBQSR`  Disable BQSR
- `--skipIdentito` Disable Identito vigilance
- `--skipMultiqc` Disable MultiQC
- `--skipPreseq` Disable Preseq
- `--skipQC` Disable all QC tools
- `--skipMutectContamination` Disable mutect2 `CalculateContamination` step

#### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

#### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.
You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

#### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

#### `--maxMemory`

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

#### `--maxTime`

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

#### `--maxCpus`

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

#### `--multiqcConfig`

Specify a path to a custom MultiQC configuration file.

## Job resources

## Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time (see the `conf/process.conf` file).
For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.
