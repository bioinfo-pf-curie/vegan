# Usage

## Table of contents

* [Introduction](#general-nextflow-info)
* [Running the pipeline](#running-the-pipeline)
* [Main arguments](#main-arguments)
    * [`-profile`](#-profile-single-dash)
        * [`conda`](#conda)
        * [`multiconda`](#multiconda)
        * [`singularity`](#singularity)
        * [`docker`](#docker)
        * ['cluster'](#cluster)
        * [`test`](#test)
    * [`--reads`](#--reads)
	* [`--samplePlan`](#--samplePlan)
    * [`--design`](#--design)
	* [`--singleEnd`](#--singleend)
    * [`--noIntervals`](#--noIntervals)
    * [`--step`](#--step)
    * [`--targetBED`](#--targetBED)
    * [`--tools`](#--tools)
    * [`--SNVFilters`](#--SNVFilters)
    * [`--SVFilters`](#--SVFilters)
    * [`--condaCacheDir`](#--condaCacheDir)
    * [`--genomeAnnotationPath`](#--genomeAnnotationPath)

* [Reference genomes](#reference-genomes)
    * [`--genome`](#--genome)
* [Job resources](#job-resources)
* [Automatic resubmission](#automatic-resubmission)
* [Other command line parameters](#other-command-line-parameters)
    * [`--skip*`](#--skip*)
	* [`--outDir`](#--outDir)
    * [`--email`](#--email)
    * [`-name`](#-name-single-dash)
    * [`-resume`](#-resume-single-dash)
    * [`-c`](#-c-single-dash)
    * [`--max_memory`](#--max_memory)
    * [`--max_time`](#--max_time)
    * [`--max_cpus`](#--max_cpus)
    * [`--multiqc_config`](#--multiqc_config)

## General Nextflow info

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline
The typical command for running the pipeline is as follows:
```bash
nextflow run main.nf --reads '*_R{1,2}.fastq.gz' -profile 'singularity'
```

This will launch the pipeline with the `singularity` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

You can change the output director using the `--outDir/-w` options.

## Main arguments

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments. Note that multiple profiles can be loaded, for example: `-profile singularity,cluster` - the order of arguments is important!

If `-profile` is not specified at all the pipeline will be run locally and expects all software to be installed and available on the `PATH`.

* `conda`
    * A generic configuration profile to be used with [conda](https://conda.io/docs/)
    * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `multiconda`
    * A generic configuration profile to be used with [conda](https://conda.io/docs/)
    * Uses one conda environment per general task / tool
* `singularity`
    * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
    * Use the singularity images available on the cluster
* `docker`
    * A generic configuration profile to be used with [Docker](https://www.docker.com/)
* `cluster`
    * Run the workflow on the computational cluster
* `test`
    * A profile with a complete configuration for automated testing
    * Includes links to test data so needs no other parameters

### `--reads`

Use this to specify the location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*_{1,2}.fastq'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

### `--samplePlan`

Use this to specify a sample plan file instead of a regular expression to find fastq files. For example :

```bash
--samplePlan 'path/to/data/sample_plan.csv'
```

The sample plan is a csv file with the following information :

Sample ID | Sample Name | Path to R1 fastq file | Path to R2 fastq file

### `--design`

Use this to specify a design file to list all experimental samples, their IDs, the associated germinal sample, the sex of the patient and the status (tumor / normal). For example :

```bash
--design 'path/to/design.csv'
```

The design file is a csv file with the following information :

SAMPLE_ID | GERMLINE ID | SAMPLE_NAME | SEX | STATUS

### `--singleEnd`

By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--singleEnd` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example:

```bash
--singleEnd --reads '*.fastq.gz'
```

It is not possible to run a mixture of single-end and paired-end files in one run.

### `--noIntervals`
Disable usage of intervals file, and disable automatic generation of intervals file when none are provided.

For WGS samples, splitting in interval is advised to reduce running time.

### `--step`

This parameter specify the starting step of the pipeline. Several entry point are available. Each step require a specific sample_plan.

**Available:** mapping, recalibrate, variantcalling, annotate.

### `--targetBED`

Specify a target BED file for targeted or whole exome sequencing

### `--tools`

Specify the tools to use for variant calling and downstream steps.

**Available:** facets, ascat, haplotypecaller, manta, mutect2, snpeff

### `--SNVFilters`

Specify which filter to use for SNV calling.

**Available:** mapq, duplicates, singleton, multihits

**Default:** mapq and duplicates

### `--SVFilters`

Specify which filter to use for SV calling.

**Available:** mapq, duplicates, singleton, multihits

**Default:** duplicates

### `--condaCacheDir`

Specify the path to store conda environment create by Nextflow

### `--genomeAnnotationPath`

Specify the path to the annotations files required by the pipeline

## Reference genomes

The pipeline config files come bundled with paths to the genomes reference files.

### `--genome`

There are different species supported in the genomes references file. To run the pipeline, you must specify which to use with the `--genome` flag.

You can find the keys to specify the genomes in the [genomes config file](../conf/genomes.config). Common genomes that are supported are:

* Human
  * `--genome hg19`
  * `--genome hg38`
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
- `--gtf` - Path to GTF file

## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time (see the `conf/process.conf` file).
For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

## Other command line parameters

### `--skip*`

The pipeline is made with a few *skip* options that allow to skip optional steps in the workflow.
The following options can be used:

- `--skipBQSR`  Disable BQSR
- `--skipIdentito` Disable Identito vigilance
- `--skipMultiqc` Disable MultiQC
- `--skipPreseq` Disable Preseq
- `--skipQC` Specify which QC tools to skip from bamqc, fastqc, multiqc, samtoolsstats, versions
- `--skipMutectContamination` Disable mutect2 `CalculateContamination` step

### `--outDir`

The output directory where the results will be saved.

### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.
You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--max_memory`

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--multiqc_config`

Specify a path to a custom MultiQC configuration file.
