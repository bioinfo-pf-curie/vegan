#!/usr/bin/env nextflow

/*
Copyright Institut Curie 2019
This software is a computer program whose purpose is to analyze high-throughput sequencing data.
You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.
Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data.
The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.

This script is based on the nf-core guidelines. See https://nf-co.re/ for more information
*/

/*
========================================================================================
                         data-analysis_demo
========================================================================================
 data-analysis_demo analysis Pipeline.
 #### Homepage / Documentation
 ssh://git@gitlab.curie.fr:2222/plarosa/data-analysis_demo.git
----------------------------------------------------------------------------------------
*/
/*
================================================================================
                                 EUCANCAN/nf-wgswes 
================================================================================
Started Febery 2020.
--------------------------------------------------------------------------------
nf-wgswes:
  An open-source analysis pipeline to detect germline or somatic variants
  from whole genome or targeted sequencing
--------------------------------------------------------------------------------
 @Homepage
 https://gitlab.curie.fr/data-analysis/nf-wgswes 
--------------------------------------------------------------------------------
 @Documentation
 https://gitlab.curie.fr/data-analysis/nf-wgswes/README.md
--------------------------------------------------------------------------------
*/

def helpMessage() {
    //log.info nfcoreHeader()
    if ("${workflow.manifest.version}" =~ /dev/ ){
       log.info devMessageFile.text
    }

    log.info"""
    data-analysis_demo v${workflow.manifest.version}
    ======================================================================

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run main.nf --input samples-WGS.tsv -profile multiconda

    Mandatory arguments:
        --input                     Path to input TSV file on mapping, recalibrate and variantcalling steps
                                    Multiple TSV files can be specified with quotes
                                    Works also with the path to a directory on mapping step with a single germline sample only
                                    Alternatively, path to VCF input file on annotate step
                                    Multiple VCF files can be specified with quotes
        -profile                    Configuration profile to use
                                    Can use multiple (comma separated)
                                    Available: conda, docker, singularity, test and more

    Options:
        --genome                    Name of iGenomes reference
        --noGVCF                    No g.vcf output from HaplotypeCaller
        --noStrelkaBP               Will not use Manta candidateSmallIndels for Strelka as Best Practice
        --no_intervals              Disable usage of intervals
        --nucleotidesPerSecond      To estimate interval size
                                    Default: 1000.0
        --targetBED                 Target BED file for targeted or whole exome sequencing
        --step                      Specify starting step
                                    Available: Mapping, Recalibrate, VariantCalling, Annotate
                                    Default: Mapping
        --tools                     Specify tools to use for variant calling:
                                    Available: ASCAT, ControlFREEC, FreeBayes, HaplotypeCaller
                                    Manta, mpileup, Mutect2, Strelka, TIDDIT
                                    and/or for annotation:
                                    snpEff, VEP, merge
                                    Default: None
        --skipQC                    Specify which QC tools to skip when running Sarek
                                    Available: all, bamQC, BCFtools, FastQC, MultiQC, samtools, vcftools, versions
                                    Default: None
        --annotateTools             Specify from which tools Sarek will look for VCF files to annotate, only for step annotate
                                    Available: HaplotypeCaller, Manta, Mutect2, Strelka, TIDDIT
                                    Default: None
        --annotation_cache          Enable the use of cache for annotation, to be used with --snpEff_cache and/or --vep_cache
        --snpEff_cache              Specity the path to snpEff cache, to be used with --annotation_cache
        --vep_cache                 Specity the path to VEP cache, to be used with --annotation_cache
        --pon                       panel-of-normals VCF (bgzipped, indexed). See: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_mutect_CreateSomaticPanelOfNormals.php
        --pon_index                 index of pon panel-of-normals VCF

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
        --acLoci                    acLoci file
        --acLociGC                  acLoci GC file
        --bwaIndex                  bwa indexes
                                    If none provided, will be generated automatically from the fasta reference
        --dbsnp                     dbsnp file
        --dbsnpIndex                dbsnp index
                                    If none provided, will be generated automatically if a dbsnp file is provided
        --dict                      dict from the fasta reference
                                    If none provided, will be generated automatically from the fasta reference
        --fasta                     fasta reference
        --fastafai                  reference index
                                    If none provided, will be generated automatically from the fasta reference
        --germlineResource          Germline Resource File
        --germlineResourceIndex     Germline Resource Index
                                    If none provided, will be generated automatically if a germlineResource file is provided
        --intervals                 intervals
                                    If none provided, will be generated automatically from the fasta reference
                                    Use --no_intervals to disable automatic generation
        --knownIndels               knownIndels file
        --knownIndelsIndex          knownIndels index
                                    If none provided, will be generated automatically if a knownIndels file is provided
        --species                   species for VEP
        --snpeffDb                  snpeffDb version
        --vepCacheVersion           VEP Cache version

    Other options:
        --outputDir                 The output directory where the results will be saved
        --sequencing_center         Name of sequencing center to be displayed in BAM file
        --multiqc_config            Specify a custom config file for MultiQC
        --monochrome_logs           Logs will be without colors
        --email                     Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
        --maxMultiqcEmailFileSize   Theshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
        -name                       Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    AWSBatch options:
        --awsqueue                  The AWSBatch JobQueue that needs to be set when running on AWSBatch
        --awsregion                 The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help) exit 0, helpMessage()

// Print warning message
if (params.noReports) log.warn "The params `--noReports` is deprecated -- it will be removed in a future release.\n\tPlease check: https://github.com/nf-core/sarek/blob/master/docs/usage.md#--skipQC"
if (params.annotateVCF) log.warn "The params `--annotateVCF` is deprecated -- it will be removed in a future release.\n\tPlease check: https://github.com/nf-core/sarek/blob/master/docs/usage.md#--input"
if (params.genomeDict) log.warn "The params `--genomeDict` is deprecated -- it will be removed in a future release.\n\tPlease check: https://github.com/nf-core/sarek/blob/master/docs/usage.md#--dict"
if (params.genomeFile) log.warn "The params `--genomeFile` is deprecated -- it will be removed in a future release.\n\tPlease check: https://github.com/nf-core/sarek/blob/master/docs/usage.md#--fasta"
if (params.genomeIndex) log.warn "The params `--genomeIndex` is deprecated -- it will be removed in a future release.\n\tPlease check: https://github.com/nf-core/sarek/blob/master/docs/usage.md#--fastaFai"
if (params.sample) log.warn "The params `--sample` is deprecated -- it will be removed in a future release.\n\tPlease check: https://github.com/nf-core/sarek/blob/master/docs/usage.md#--input"
if (params.sampleDir) log.warn "The params `--sampleDir` is deprecated -- it will be removed in a future release.\n\tPlease check: https://github.com/nf-core/sarek/blob/master/docs/usage.md#--input"

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

stepList = defineStepList()
step = params.step ? params.step.toLowerCase() : ''

// Handle deprecation
if (step == 'preprocessing') step = 'mapping'

if (step.contains(',')) exit 1, 'You can choose only one step, see --help for more information'
if (!checkParameterExistence(step, stepList)) exit 1, "Unknown step ${step}, see --help for more information"

toolList = defineToolList()
tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase()} : []
if (!checkParameterList(tools, toolList)) exit 1, 'Unknown tool(s), see --help for more information'

skipQClist = defineSkipQClist()
skipQC = params.skipQC ? params.skipQC == 'all' ? skipQClist : params.skipQC.split(',').collect{it.trim().toLowerCase()} : []
if (!checkParameterList(skipQC, skipQClist)) exit 1, 'Unknown QC tool(s), see --help for more information'

annoList = defineAnnoList()
annotateTools = params.annotateTools ? params.annotateTools.split(',').collect{it.trim().toLowerCase()} : []
if (!checkParameterList(annotateTools,annoList)) exit 1, 'Unknown tool(s) to annotate, see --help for more information'

// Handle deprecation
if (params.noReports) skipQC = skipQClist

// Has the run name been specified by the user?
// This has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) custom_runName = workflow.runName

if (workflow.profile == 'awsbatch') {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outputDir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outputDir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (workflow.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
ch_output_docs = Channel.fromPath("${baseDir}/docs/output.md")

tsvPath = null
if (params.input && (hasExtension(params.input, "tsv") || hasExtension(params.input, "vcf") || hasExtension(params.input, "vcf.gz"))) tsvPath = params.input
if (params.input && (hasExtension(params.input, "vcf") || hasExtension(params.input, "vcf.gz"))) step = "annotate"

// Handle deprecation
if (params.annotateVCF) tsvPath = params.annotateVCF
if (params.sample) tsvPath = params.sample
if (params.sampleDir) tsvPath = params.sampleDir

// If no input file specified, trying to get TSV files corresponding to step in the TSV directory
// only for steps recalibrate and variantCalling
if (!params.input && step != 'mapping' && step != 'annotate') {
	tsvPath = step == 'recalibrate' ? "${params.outputDir}/Preprocessing/TSV/duplicateMarked.tsv" : "${params.outputDir}/Preprocessing/TSV/recalibrated.tsv"
}

inputSample = Channel.empty()
if (tsvPath) {
    tsvFile = file(tsvPath)
    switch (step) {
        case 'mapping': inputSample = extractFastq(tsvFile); break
        case 'recalibrate': inputSample = extractRecal(tsvFile); break
        case 'variantcalling': inputSample = extractBam(tsvFile); break
        case 'annotate': break
        default: exit 1, "Unknown step ${step}"
    }
} else if (params.input && !hasExtension(params.input, "tsv")) {
    log.info "No TSV file"
    if (step != 'mapping') exit 1, 'No other step than "mapping" support a dir as an input'
    log.info "Reading ${params.input} directory"
    inputSample = extractFastqFromDir(params.input)
    (inputSample, fastqTMP) = inputSample.into(2)
    fastqTMP.toList().subscribe onNext: {
        if (it.size() == 0) exit 1, "No FASTQ files found in --input directory '${params.input}'"
    }
    tsvFile = params.input  // used in the reports
} else if (tsvPath && step == 'annotate') {
    log.info "Annotating ${tsvPath}"
} else if (step == 'annotate') {
    log.info "Trying automatic annotation on file in the VariantCalling directory"
} else exit 1, 'No sample were defined, see --help'

(genderMap, statusMap, inputSample) = extractInfos(inputSample)

/*
================================================================================
                               CHECKING REFERENCES
================================================================================
*/

// Initialize each params in params.genomes, catch the command line first if it was defined
// params.fasta has to be the first one
params.fasta = params.genome && !('annotate' in step) ? params.genomes[params.genome].fasta ?: null : null
// The rest can be sorted
params.acLoci = params.genome && 'ascat' in tools ? params.genomes[params.genome].acLoci ?: null : null
params.acLociGC = params.genome && 'ascat' in tools ? params.genomes[params.genome].acLociGC ?: null : null
params.bwaIndex = params.genome && params.fasta && 'mapping' in step ? params.genomes[params.genome].bwaIndex ?: null : null
params.chrDir = params.genome && 'controlfreec' in tools ? params.genomes[params.genome].chrDir ?: null : null
params.chrLength = params.genome && 'controlfreec' in tools ? params.genomes[params.genome].chrLength ?: null : null
params.dbsnp = params.genome && ('mapping' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools) ? params.genomes[params.genome].dbsnp ?: null : null
params.dbsnpIndex = params.genome && params.dbsnp ? params.genomes[params.genome].dbsnpIndex ?: null : null
params.dict = params.genome && params.fasta ? params.genomes[params.genome].dict ?: null : null
params.fastaFai = params.genome && params.fasta ? params.genomes[params.genome].fastaFai ?: null : null
params.germlineResource = params.genome && 'mutect2' in tools ? params.genomes[params.genome].germlineResource ?: null : null
params.germlineResourceIndex = params.genome && params.germlineResource ? params.genomes[params.genome].germlineResourceIndex ?: null : null
params.intervals = params.genome && !('annotate' in step) ? params.genomes[params.genome].intervals ?: null : null
params.knownIndels = params.genome && 'mapping' in step ? params.genomes[params.genome].knownIndels ?: null : null
params.knownIndelsIndex = params.genome && params.knownIndels ? params.genomes[params.genome].knownIndelsIndex ?: null : null
params.snpeffDb = params.genome && 'snpeff' in tools ? params.genomes[params.genome].snpeffDb ?: null : null
params.species = params.genome && 'vep' in tools ? params.genomes[params.genome].species ?: null : null
params.vepCacheVersion = params.genome && 'vep' in tools ? params.genomes[params.genome].vepCacheVersion ?: null : null

// Initialize channels based on params
ch_acLoci = params.acLoci && 'ascat' in tools ? Channel.value(file(params.acLoci)) : "null"
ch_acLociGC = params.acLociGC && 'ascat' in tools ? Channel.value(file(params.acLociGC)) : "null"
ch_chrDir = params.chrDir && 'controlfreec' in tools ? Channel.value(file(params.chrDir)) : "null"
ch_chrLength = params.chrLength && 'controlfreec' in tools ? Channel.value(file(params.chrLength)) : "null"
ch_dbsnp = params.dbsnp && ('mapping' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools) ? Channel.value(file(params.dbsnp)) : "null"
ch_fasta = params.fasta && !('annotate' in step) ? Channel.value(file(params.fasta)) : "null"
ch_fastaFai = params.fastaFai && !('annotate' in step) ? Channel.value(file(params.fastaFai)) : "null"
ch_germlineResource = params.germlineResource && 'mutect2' in tools ? Channel.value(file(params.germlineResource)) : "null"
ch_intervals = params.intervals && !params.no_intervals && !('annotate' in step) ? Channel.value(file(params.intervals)) : "null"

// knownIndels is currently a list of file for smallGRCh37, so transform it in a channel
li_knownIndels = []
if (params.knownIndels && ('mapping' in step)) params.knownIndels.each { li_knownIndels.add(file(it)) }
ch_knownIndels = params.knownIndels && params.genome == 'smallGRCh37' ? Channel.value(li_knownIndels.collect()) : params.knownIndels ? Channel.value(file(params.knownIndels)) : "null"

ch_snpEff_cache = params.snpEff_cache ? Channel.value(file(params.snpEff_cache)) : "null"
ch_snpeffDb = params.snpeffDb ? Channel.value(params.snpeffDb) : "null"
ch_vepCacheVersion = params.vepCacheVersion ? Channel.value(params.vepCacheVersion) : "null"
ch_vep_cache = params.vep_cache ? Channel.value(file(params.vep_cache)) : "null"

// Optional files, not defined within the params.genomes[params.genome] scope
ch_cadd_InDels = params.cadd_InDels ? Channel.value(file(params.cadd_InDels)) : "null"
ch_cadd_InDels_tbi = params.cadd_InDels_tbi ? Channel.value(file(params.cadd_InDels_tbi)) : "null"
ch_cadd_WG_SNVs = params.cadd_WG_SNVs ? Channel.value(file(params.cadd_WG_SNVs)) : "null"
ch_cadd_WG_SNVs_tbi = params.cadd_WG_SNVs_tbi ? Channel.value(file(params.cadd_WG_SNVs_tbi)) : "null"
ch_pon = params.pon ? Channel.value(file(params.pon)) : "null"
ch_targetBED = params.targetBED ? Channel.value(file(params.targetBED)) : "null"

/*
================================================================================
                                PRINTING SUMMARY
================================================================================
*/

// Header log info
// log.info nfcoreHeader()
def summary = [:]
if (workflow.revision)          summary['Pipeline Release']    = workflow.revision
summary['Run Name']          = custom_runName ?: workflow.runName
summary['Max Resources']     = "${params.maxMemory} memory, ${params.maxCpus} cpus, ${params.maxTime} time per job"
if (workflow.containerEngine)   summary['Container']         = "${workflow.containerEngine} - ${workflow.container}"
if (params.input)               summary['Input']             = params.input
if (params.targetBED)           summary['Target BED']        = params.targetBED
if (step)                       summary['Step']              = step
if (params.tools)               summary['Tools']             = tools.join(', ')
if (params.skipQC)              summary['QC tools skip']     = skipQC.join(', ')

if (params.no_intervals && step != 'annotate') summary['Intervals']         = 'Do not use'
if ('haplotypecaller' in tools)                summary['GVCF']              = params.noGVCF ? 'No' : 'Yes'
if ('strelka' in tools && 'manta' in tools )   summary['Strelka BP']        = params.noStrelkaBP ? 'No' : 'Yes'
if (params.sequencing_center)                  summary['Sequenced by']      = params.sequencing_center
if (params.pon && 'mutect2' in tools)          summary['Panel of normals']  = params.pon

summary['Save Genome Index'] = params.saveGenomeIndex ? 'Yes' : 'No'
summary['Nucleotides/s']     = params.nucleotidesPerSecond
summary['Output dir']        = params.outputDir
summary['Launch dir']        = workflow.launchDir
summary['Working dir']       = workflow.workDir
summary['Script dir']        = workflow.projectDir
summary['User']              = workflow.userName
summary['genome']            = params.genome

if (params.fasta)                 summary['fasta']                 = params.fasta
if (params.fastaFai)              summary['fastaFai']              = params.fastaFai
if (params.dict)                  summary['dict']                  = params.dict
if (params.bwaIndex)              summary['bwaIndex']              = params.bwaIndex
if (params.germlineResource)      summary['germlineResource']      = params.germlineResource
if (params.germlineResourceIndex) summary['germlineResourceIndex'] = params.germlineResourceIndex
if (params.intervals)             summary['intervals']             = params.intervals
if (params.acLoci)                summary['acLoci']                = params.acLoci
if (params.acLociGC)              summary['acLociGC']              = params.acLociGC
if (params.chrDir)                summary['chrDir']                = params.chrDir
if (params.chrLength)             summary['chrLength']             = params.chrLength
if (params.dbsnp)                 summary['dbsnp']                 = params.dbsnp
if (params.dbsnpIndex)            summary['dbsnpIndex']            = params.dbsnpIndex
if (params.knownIndels)           summary['knownIndels']           = params.knownIndels
if (params.knownIndelsIndex)      summary['knownIndelsIndex']      = params.knownIndelsIndex
if (params.snpeffDb)              summary['snpeffDb']              = params.snpeffDb
if (params.species)               summary['species']               = params.species
if (params.vepCacheVersion)       summary['vepCacheVersion']       = params.vepCacheVersion
if (params.species)               summary['species']               = params.species
if (params.snpEff_cache)          summary['snpEff_cache']          = params.snpEff_cache
if (params.vep_cache)             summary['vep_cache']             = params.vep_cache

if (workflow.profile == 'awsbatch') {
    summary['AWS Region']        = params.awsregion
    summary['AWS Queue']         = params.awsqueue
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description)  summary['Config Description']  = params.config_profile_description
if (params.config_profile_contact)      summary['Config Contact']      = params.config_profile_contact
if (params.config_profile_url)          summary['Config URL']          = params.config_profile_url
if (params.email) {
    summary['E-mail Address']        = params.email
    summary['MultiQC maxsize']       = params.maxMultiqcEmailFileSize
}
log.info summary.collect { k, v -> "${k.padRight(18)}: $v" }.join("\n")
if (params.monochrome_logs) log.info "----------------------------------------------------"
else log.info "\033[2m----------------------------------------------------\033[0m"

// Check the hostnames against configured profiles
checkHostname()

/*
 * Parse software version numbers
 */
process GetSoftwareVersions {
    label 'onlyLinux'

    publishDir path:"${params.outputDir}/pipeline_info", mode: params.publishDirMode

    output:
        file 'software_versions_mqc.yaml' into yamlSoftwareVersion

    when: !('versions' in skipQC)

    script:
    """
    cp ${baseDir}/assets/software_versions/*.txt .
    echo "${workflow.manifest.version}" &> v_pipeline.txt 2>&1 || true
    echo "${workflow.nextflow.version}" &> v_nextflow.txt 2>&1 || true

    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

yamlSoftwareVersion = yamlSoftwareVersion.dump(tag:'SOFTWARE VERSIONS')

/*
================================================================================
                                BUILDING INDEXES
================================================================================
*/

// And then initialize channels based on params or indexes that were just built

process BuildBWAindexes {
    label 'bwa'
    tag {fasta}

    publishDir params.outputDir, mode: params.publishDirMode,
        saveAs: {params.saveGenomeIndex ? "reference_genome/BWAIndex/${it}" : null }

    input:
        file(fasta) from ch_fasta

    output:
        file("${fasta}.*") into bwaIndexes

    when: !(params.bwaIndex) && params.fasta && 'mapping' in step

    script:
    """
    bwa index ${fasta}
    """
}

ch_bwaIndex = params.bwaIndex ? Channel.value(file(params.bwaIndex)) : bwaIndexes

process BuildDict {
    label 'gatk'
    tag {fasta}

    publishDir params.outputDir, mode: params.publishDirMode,
        saveAs: {params.saveGenomeIndex ? "reference_genome/${it}" : null }

    input:
        file(fasta) from ch_fasta

    output:
        file("${fasta.baseName}.dict") into dictBuilt

    when: !(params.dict) && params.fasta && !('annotate' in step)

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        CreateSequenceDictionary \
        --REFERENCE ${fasta} \
        --OUTPUT ${fasta.baseName}.dict
    """
}

ch_dict = params.dict ? Channel.value(file(params.dict)) : dictBuilt

process BuildFastaFai {
    label 'samtools'
    tag {fasta}

    publishDir params.outputDir, mode: params.publishDirMode,
        saveAs: {params.saveGenomeIndex ? "reference_genome/${it}" : null }

    input:
        file(fasta) from ch_fasta

    output:
        file("${fasta}.fai") into fastaFaiBuilt

    when: !(params.fastaFai) && params.fasta && !('annotate' in step)

    script:
    """
    samtools faidx ${fasta}
    """
}

ch_fastaFai = params.fastaFai ? Channel.value(file(params.fastaFai)) : fastaFaiBuilt

process BuildDbsnpIndex {
    label 'tabix'
    tag {dbsnp}

    publishDir params.outputDir, mode: params.publishDirMode,
        saveAs: {params.saveGenomeIndex ? "reference_genome/${it}" : null }

    input:
        file(dbsnp) from ch_dbsnp

    output:
        file("${dbsnp}.tbi") into dbsnpIndexBuilt

    when: !(params.dbsnpIndex) && params.dbsnp && ('mapping' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools)

    script:
    """
    tabix -p vcf ${dbsnp}
    """
}

ch_dbsnpIndex = params.dbsnp ? params.dbsnpIndex ? Channel.value(file(params.dbsnpIndex)) : dbsnpIndexBuilt : "null"

process BuildGermlineResourceIndex {
    label 'tabix'
    tag {germlineResource}

    publishDir params.outputDir, mode: params.publishDirMode,
        saveAs: {params.saveGenomeIndex ? "reference_genome/${it}" : null }

    input:
        file(germlineResource) from ch_germlineResource

    output:
        file("${germlineResource}.tbi") into germlineResourceIndexBuilt

    when: !(params.germlineResourceIndex) && params.germlineResource && 'mutect2' in tools

    script:
    """
    tabix -p vcf ${germlineResource}
    """
}

ch_germlineResourceIndex = params.germlineResource ? params.germlineResourceIndex ? Channel.value(file(params.germlineResourceIndex)) : germlineResourceIndexBuilt : "null"

process BuildKnownIndelsIndex {
    label 'tabix'
    tag {knownIndels}

    publishDir params.outputDir, mode: params.publishDirMode,
        saveAs: {params.saveGenomeIndex ? "reference_genome/${it}" : null }

    input:
        each file(knownIndels) from ch_knownIndels

    output:
        file("${knownIndels}.tbi") into knownIndelsIndexBuilt

    when: !(params.knownIndelsIndex) && params.knownIndels && 'mapping' in step

    script:
    """
    tabix -p vcf ${knownIndels}
    """
}

ch_knownIndelsIndex = params.knownIndels ? params.knownIndelsIndex ? Channel.value(file(params.knownIndelsIndex)) : knownIndelsIndexBuilt.collect() : "null"

process BuildPonIndex {
    label 'tabix'
    tag {pon}

    publishDir params.outputDir, mode: params.publishDirMode,
        saveAs: {params.saveGenomeIndex ? "reference_genome/${it}" : null }

    input:
        file(pon) from ch_pon

    output:
        file("${pon}.tbi") into ponIndexBuilt

    when: !(params.pon_index) && params.pon && ('tnscope' in tools || 'mutect2' in tools)

    script:
    """
    tabix -p vcf ${pon}
    """
}

ch_ponIndex = params.pon_index ? Channel.value(file(params.pon_index)) : ponIndexBuilt

process BuildIntervals {
  label 'onlyLinux'
  tag {fastaFai}

  publishDir params.outputDir, mode: params.publishDirMode,
    saveAs: {params.saveGenomeIndex ? "reference_genome/${it}" : null }

  input:
    file(fastaFai) from ch_fastaFai

  output:
    file("${fastaFai.baseName}.bed") into intervalBuilt

  when: !(params.intervals) && !('annotate' in step) && !(params.no_intervals)

  script:
  """
  awk -v FS='\t' -v OFS='\t' '{ print \$1, \"0\", \$2 }' ${fastaFai} > ${fastaFai.baseName}.bed
  """
}

ch_intervals = params.no_intervals ? "null" : params.intervals && !('annotate' in step) ? Channel.value(file(params.intervals)) : intervalBuilt

/*
================================================================================
                                  PREPROCESSING
================================================================================
*/

// STEP 0: CREATING INTERVALS FOR PARALLELIZATION (PREPROCESSING AND VARIANT CALLING)

process CreateIntervalBeds {
    label 'onlyLinux' 
    tag {intervals.fileName}

    input:
        file(intervals) from ch_intervals

    output:
        file '*.bed' into bedIntervals mode flatten

    when: (!params.no_intervals) && step != 'annotate'

    script:
    // If the interval file is BED format, the fifth column is interpreted to
    // contain runtime estimates, which is then used to combine short-running jobs
    if (hasExtension(intervals, "bed"))
        """
        awk -vFS="\t" '{
          t = \$5  # runtime estimate
          if (t == "") {
            # no runtime estimate in this row, assume default value
            t = (\$3 - \$2) / ${params.nucleotidesPerSecond}
          }
          if (name == "" || (chunk > 600 && (chunk + t) > longest * 1.05)) {
            # start a new chunk
            name = sprintf("%s_%d-%d.bed", \$1, \$2+1, \$3)
            chunk = 0
            longest = 0
          }
          if (t > longest)
            longest = t
          chunk += t
          print \$0 > name
        }' ${intervals}
        """
    else if (hasExtension(intervals, "interval_list"))
        """
        grep -v '^@' ${intervals} | awk -vFS="\t" '{
          name = sprintf("%s_%d-%d", \$1, \$2, \$3);
          printf("%s\\t%d\\t%d\\n", \$1, \$2-1, \$3) > name ".bed"
        }'
        """
    else
        """
        awk -vFS="[:-]" '{
          name = sprintf("%s_%d-%d", \$1, \$2, \$3);
          printf("%s\\t%d\\t%d\\n", \$1, \$2-1, \$3) > name ".bed"
        }' ${intervals}
        """
}

bedIntervals = bedIntervals
    .map { intervalFile ->
        def duration = 0.0
        for (line in intervalFile.readLines()) {
            final fields = line.split('\t')
            if (fields.size() >= 5) duration += fields[4].toFloat()
            else {
                start = fields[1].toInteger()
                end = fields[2].toInteger()
                duration += (end - start) / params.nucleotidesPerSecond
            }
        }
        [duration, intervalFile]
        }.toSortedList({ a, b -> b[0] <=> a[0] })
    .flatten().collate(2)
    .map{duration, intervalFile -> intervalFile}

bedIntervals = bedIntervals.dump(tag:'bedintervals')

if (params.no_intervals && step != 'annotate') bedIntervals = Channel.from(file("no_intervals.bed"))

(intBaseRecalibrator, intApplyBQSR, intHaplotypeCaller, intMpileup, bedIntervals) = bedIntervals.into(5)

// PREPARING CHANNELS FOR PREPROCESSING AND QC

inputBam = Channel.create()
inputPairReads = Channel.create()

if (step in ['recalibrate', 'variantcalling', 'annotate']) {
    inputBam.close()
    inputPairReads.close()
} else inputSample.choice(inputPairReads, inputBam) {hasExtension(it[3], "bam") ? 1 : 0}

(inputBam, inputBamFastQC) = inputBam.into(2)

// Removing inputFile2 wich is null in case of uBAM
inputBamFastQC = inputBamFastQC.map {
    idPatient, idSample, idRun, inputFile1, inputFile2 ->
    [idPatient, idSample, idRun, inputFile1]
}

if (params.split_fastq){
    inputPairReads = inputPairReads
        // newly splitfastq are named based on split, so the name is easier to catch
        .splitFastq(by: params.split_fastq, compress:true, file:"split", pe:true)
        .map {idPatient, idSample, idRun, reads1, reads2 ->
            // The split fastq read1 is the 4th element (indexed 3) its name is split_3
            // The split fastq read2's name is split_4
            // It's followed by which split it's acutally based on the mother fastq file
            // Index start at 1
            // Extracting the index to get a new IdRun
            splitIndex = reads1.fileName.toString().minus("split_3.").minus(".gz")
            newIdRun = idRun + "_" + splitIndex
            // Giving the files a new nice name
            newReads1 = file("${idSample}_${newIdRun}_R1.fastq.gz")
            newReads2 = file("${idSample}_${newIdRun}_R2.fastq.gz")
            [idPatient, idSample, newIdRun, reads1, reads2]}
}

inputPairReads = inputPairReads.dump(tag:'INPUT')

(inputPairReads, inputPairReadsFastQC) = inputPairReads.into(2)

// STEP 0.5: QC ON READS

// TODO: Use only one process for FastQC for FASTQ files and uBAM files
// FASTQ and uBAM files are renamed based on the sample name

process FastQCFQ { 
    label 'fastqc'
    label 'FastQC'
    label 'cpus_2'

    tag {idPatient + "-" + idRun}

    publishDir "${params.outputDir}/Reports/${idSample}/FastQC/${idSample}_${idRun}", mode: params.publishDirMode

    input:
        set idPatient, idSample, idRun, file("${idSample}_${idRun}_R1.fastq.gz"), file("${idSample}_${idRun}_R2.fastq.gz") from inputPairReadsFastQC

    output:
        file("*.{html,zip}") into fastQCFQReport

    when: !('fastqc' in skipQC)
    
    script:
    """
    fastqc -t 2 -q ${idSample}_${idRun}_R1.fastq.gz ${idSample}_${idRun}_R2.fastq.gz
    """
}

process FastQCBAM { 
    label 'fastqc'
    label 'FastQC'
    label 'cpus_2'

    tag {idPatient + "-" + idRun}

    publishDir "${params.outputDir}/Reports/${idSample}/FastQC/${idSample}_${idRun}", mode: params.publishDirMode

    input:
        set idPatient, idSample, idRun, file("${idSample}_${idRun}.bam") from inputBamFastQC

    output:
        file("*.{html,zip}") into fastQCBAMReport

    when: !('fastqc' in skipQC)

    script:
    """
    fastqc -t 2 -q ${idSample}_${idRun}.bam
    """
}

fastQCReport = fastQCFQReport.mix(fastQCBAMReport)

fastQCReport = fastQCReport.dump(tag:'FastQC')

// STEP 1: MAPPING READS TO REFERENCE GENOME WITH BWA MEM

inputPairReads = inputPairReads.dump(tag:'INPUT')

inputPairReads = inputPairReads.mix(inputBam)

process MapReads {
    label 'gatk_bwa_samtools' 
    label 'cpus_max'

    tag {idPatient + "-" + idRun}

    input:
        set idPatient, idSample, idRun, file(inputFile1), file(inputFile2) from inputPairReads
        file(bwaIndex) from ch_bwaIndex
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fastaFai

    output:
        set idPatient, idSample, idRun, file("${idSample}_${idRun}.bam") into bamMapped
        set idPatient, val("${idSample}_${idRun}"), file("${idSample}_${idRun}.bam") into bamMappedBamQC

    script:
    // -K is an hidden option, used to fix the number of reads processed by bwa mem
    // Chunk size can affect bwa results, if not specified,
    // the number of threads can change which can give not deterministic result.
    // cf https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md
    // and https://github.com/gatk-workflows/gatk4-data-processing/blob/8ffa26ff4580df4ac3a5aa9e272a4ff6bab44ba2/processing-for-variant-discovery-gatk4.b37.wgs.inputs.json#L29
    CN = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ""
    readGroup = "@RG\\tID:${idRun}\\t${CN}PU:${idRun}\\tSM:${idSample}\\tLB:${idSample}\\tPL:illumina"
    // adjust mismatch penalty for tumor samples
    status = statusMap[idPatient, idSample]
    extra = status == 1 ? "-B 3" : ""
    convertToFastq = hasExtension(inputFile1, "bam") ? "gatk --java-options -Xmx${task.memory.toGiga()}g SamToFastq --INPUT=${inputFile1} --FASTQ=/dev/stdout --INTERLEAVE=true --NON_PF=true | \\" : ""
    input = hasExtension(inputFile1, "bam") ? "-p /dev/stdin - 2> >(tee ${inputFile1}.bwa.stderr.log >&2)" : "${inputFile1} ${inputFile2}"
    """
        ${convertToFastq}
        bwa mem ${params.bwaOptions} -R \"${readGroup}\" ${extra} -t ${task.cpus} -M ${fasta} \
        ${input} | \
        samtools sort --threads ${task.cpus} -m 2G - > ${idSample}_${idRun}.bam
    """
}

bamMapped = bamMapped.dump(tag:'Mapped BAM')
// Sort BAM whether they are standalone or should be merged

singleBam = Channel.create()
multipleBam = Channel.create()
bamMapped.groupTuple(by:[0, 1])
    .choice(singleBam, multipleBam) {it[2].size() > 1 ? 1 : 0}
singleBam = singleBam.map {
    idPatient, idSample, idRun, bam ->
    [idPatient, idSample, bam]
}
singleBam = singleBam.dump(tag:'Single BAM')


// STEP 1': MERGING BAM FROM MULTIPLE LANES

process MergeBamMapped {
    label 'samtools' 
    label 'cpus_8'

    tag {idPatient + "-" + idSample}

    input:
        set idPatient, idSample, idRun, file(bam) from multipleBam

    output:
        set idPatient, idSample, file("${idSample}.bam") into mergedBam, mergedBamU

    script:
    """
    samtools merge --threads ${task.cpus} ${idSample}.bam ${bam}
    """
}

mergedBam = mergedBam.dump(tag:'Merged BAM')

mergedBam = mergedBam.mix(singleBam)
//mergedBamU = mergedBamU.mix(singleBam)

mergedBam = mergedBam.dump(tag:'BAMs for MD')

(mergedBam, mergedBamToIndex, mergedBamU) = mergedBam.into(3)

process IndexBamFile {
    label 'samtools'  
    label 'cpus_8'

    tag {idPatient + "-" + idSample}

    input:
        set idPatient, idSample, file(bam) from mergedBamToIndex

    output:
        set idPatient, idSample, file(bam), file("*.bai") into indexedBam

    when: !params.knownIndels

    script:
    """
    samtools index ${bam}
    mv ${bam}.bai ${bam.baseName}.bai
    """
}

mapMbam = mergedBam 

// STEP 2: FILTRES : BWAMEM UNIQ
// Mapping Quality Filter
process BwaMemUniq {
    label 'samtools'
    label 'cpus_2'

    tag {idPatient + "-" + idSample}

    publishDir "${params.outputDir}/Reports/${idSample}/Uniq", mode: params.publishDirMode

    input:
        set idPatient, idSample, file(bam) from mapMbam 

    output:
        set idPatient, idSample, file("${idSample}.bam") into memUbam
        file("${idSample}.mapping.stats") into mapUReport

    when: !('uniq' in skipQC)

    script:

    """
    #removed unmapped also with -F 4
    samtools view  -@ ${task.cpus} -h -F 4 ${bam} | grep -v \"XA:Z\" | samtools view  -@ ${task.cpus} -bS > ${idSample}.temp.bam 2> ${idSample}.temp.txt 
    samtools sort -@ ${task.cpus} -o ${idSample}.bam ${idSample}.temp.bam
    samtools index ${idSample}.bam
    samtools index ${bam}

    UniqueHits=\$(samtools idxstats ${idSample}.bam |  awk '{ UNIQ_HIT+=\$3 } END { print UNIQ_HIT }')
    samtools idxstats ${bam} |  awk -v Unique_hits="\$UniqueHits" '{
    Total_reads+=\$3+\$4; Mapped_reads+=\$3; Unmapped+=\$4 } END {
          printf("Total_reads\\t%d\\nMapped_reads\\t%d\\nUnique_hits\\t%d\\nMulti_hits\\t%d\\nUnmapped\\t%d\\n.uniq(%%)\\t%.2f \\n", \
          Total_reads, Mapped_reads, Unique_hits, (Mapped_reads - Unique_hits), Unmapped, (Unique_hits*100/Total_reads))
    }' > ${idSample}.mapping.stats 
    # clean
    rm ./${idSample}.temp.* ./*.bam.bai 

    """
}

if ('uniq' in skipQC) {
 	memUbam = mergedBamU  
}

// STEP 2: FILTRES : MARKING DUPLICATES 

process MarkDuplicates { 
    label 'sambamba' 
    label 'cpus_16'

    tag {idPatient + "-" + idSample}

    publishDir params.outputDir, mode: params.publishDirMode,
        saveAs: {
            if (it == "${idSample}.bam.metrics" && 'markduplicates' in skipQC) null
            else if (it == "${idSample}.bam.metrics") "Reports/${idSample}/MarkDuplicates/${it}"
            else "Preprocessing/${idSample}/DuplicateMarked/${it}"
        }

    input:
        set idPatient, idSample, file("${idSample}.bam") from memUbam 

    output:
        set idPatient, idSample, file("${idSample}.md.bam"), file("${idSample}.md.bam.bai") into duplicateMarkedBams, duplicateMarkedBamsMQ
        file ("${idSample}.bam.metrics") into markDuplicatesReport

    when: params.knownIndels

    script:
    """

    sambamba markdup --remove-duplicates --nthreads ${task.cpus} --tmpdir . ${idSample}.bam ${idSample}.md.bam 
    echo "markDuplicatesReport"  > ${idSample}.bam.metrics

    """
}

if ('markduplicates' in skipQC) markDuplicatesReport.close()

// Mapping Quality Filter
process MapQ {
    label 'samtools'
    label 'cpus_2'

    tag {idPatient + "-" + idSample}

    publishDir "${params.outputDir}/Reports/${idSample}/MapQ", mode: params.publishDirMode
   // publishDir "${params.outputDir}/Reports/${idSample}/MapQ", pattern: '*.{bam,bam.bai}', mode: 'copy', overwrite: true


    input:
        set idPatient, idSample, file(bam), file(bai) from duplicateMarkedBamsMQ 

    output:
        set idPatient, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bam.bai") into mapQbam 
        file("${bam.baseName}.${params.mapQual}.mapping.stats") into mapQReport

    when: !('mapq' in skipQC)

    script:

    """
    samtools view -@ ${task.cpus} -q ${params.mapQual} -b ${bam} > ${idSample}.recal.bam
    samtools index ${idSample}.recal.bam 
    samtools idxstats ${idSample}.recal.bam |  awk -v id_sample="${idSample}" -v map_qual="${params.mapQual}" '{
    mapped+=\$3; unmapped+=\$4 } END {
          printf("SAMPLE\\t%s\\nNB\\t%d\\nNB_MAPPED\\t%d\\n.q%d(%%)\\t%.2f \\n", id_sample, mapped+unmapped, mapped, map_qual, (mapped*100/(mapped+unmapped)))
    }' > ${bam.baseName}.${params.mapQual}.mapping.stats

    """
}

duplicateMarkedBams = duplicateMarkedBams.dump(tag:'MD BAM')
markDuplicatesReport = markDuplicatesReport.dump(tag:'MD Report')

if ('mapq' in skipQC) {
 	mapQbam = duplicateMarkedBams 
}

(bamMD, bamMDToJoin) = mapQbam.into(2) // duplicateMarked + MapQ 

bamBaseRecalibrator = bamMD.combine(intBaseRecalibrator)

bamBaseRecalibrator = bamBaseRecalibrator.dump(tag:'BAM FOR BASERECALIBRATOR')

// STEP 3: CREATING RECALIBRATION TABLES

process BaseRecalibrator {
    label 'gatk'
    label 'cpus_1'

    tag {idPatient + "-" + idSample + "-" + intervalBed.baseName}

    input:
        set idPatient, idSample, file(bam), file(bai), file(intervalBed) from bamBaseRecalibrator
        file(dbsnp) from ch_dbsnp
        file(dbsnpIndex) from ch_dbsnpIndex
        file(fasta) from ch_fasta
        file(dict) from ch_dict
        file(fastaFai) from ch_fastaFai
        file(knownIndels) from ch_knownIndels
        file(knownIndelsIndex) from ch_knownIndelsIndex

    output:
        set idPatient, idSample, file("${prefix}${idSample}.recal.table") into tableGatherBQSRReports
        set idPatient, idSample into recalTableTSVnoInt

    when: params.knownIndels

    script:
    dbsnpOptions = params.dbsnp ? "--known-sites ${dbsnp}" : ""
    knownOptions = params.knownIndels ? knownIndels.collect{"--known-sites ${it}"}.join(' ') : ""
    prefix = params.no_intervals ? "" : "${intervalBed.baseName}_"
    intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
    // TODO: --use-original-qualities ???
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        BaseRecalibrator \
        -I ${bam} \
        -O ${prefix}${idSample}.recal.table \
        --tmp-dir /tmp \
        -R ${fasta} \
        ${intervalsOptions} \
        ${dbsnpOptions} \
        ${knownOptions} \
        --verbosity INFO
    """
}

if (!params.no_intervals) tableGatherBQSRReports = tableGatherBQSRReports.groupTuple(by:[0, 1])

tableGatherBQSRReports = tableGatherBQSRReports.dump(tag:'BQSR REPORTS')

if (params.no_intervals) {
    (tableGatherBQSRReports, tableGatherBQSRReportsNoInt) = tableGatherBQSRReports.into(2)
    recalTable = tableGatherBQSRReportsNoInt
} else recalTableTSVnoInt.close()

// STEP 3.5: MERGING RECALIBRATION TABLES

process GatherBQSRReports {
    label 'gatk'
    label 'memory_singleCPU_2_task'
    label 'cpus_2'

    tag {idPatient + "-" + idSample}

    publishDir "${params.outputDir}/Preprocessing/${idSample}/DuplicateMarked", mode: params.publishDirMode, overwrite: false

    input:
        set idPatient, idSample, file(recal) from tableGatherBQSRReports

    output:
        set idPatient, idSample, file("${idSample}.recal.table") into recalTable
        set idPatient, idSample into recalTableTSV

    when: !(params.no_intervals)

    script:
    input = recal.collect{"-I ${it}"}.join(' ')
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        GatherBQSRReports \
        ${input} \
        -O ${idSample}.recal.table \
    """
}

recalTable = recalTable.dump(tag:'RECAL TABLE')

(recalTableTSV, recalTableSampleTSV) = recalTableTSV.mix(recalTableTSVnoInt).into(2)

// Create TSV files to restart from this step
recalTableTSV.map { idPatient, idSample ->
    status = statusMap[idPatient, idSample]
    gender = genderMap[idPatient]
    bam = "${params.outputDir}/Preprocessing/${idSample}/DuplicateMarked/${idSample}.md.bam"
    bai = "${params.outputDir}/Preprocessing/${idSample}/DuplicateMarked/${idSample}.md.bai"
    recalTable = "${params.outputDir}/Preprocessing/${idSample}/DuplicateMarked/${idSample}.recal.table"
    "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\t${recalTable}\n"
}.collectFile(
    name: 'duplicateMarked.tsv', sort: true, storeDir: "${params.outputDir}/Preprocessing/TSV"
)

recalTableSampleTSV
    .collectFile(storeDir: "${params.outputDir}/Preprocessing/TSV/") {
        idPatient, idSample ->
        status = statusMap[idPatient, idSample]
        gender = genderMap[idPatient]
        bam = "${params.outputDir}/Preprocessing/${idSample}/DuplicateMarked/${idSample}.md.bam"
        bai = "${params.outputDir}/Preprocessing/${idSample}/DuplicateMarked/${idSample}.md.bai"
        recalTable = "${params.outputDir}/Preprocessing/${idSample}/DuplicateMarked/${idSample}.recal.table"
        ["duplicateMarked_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\t${recalTable}\n"]
}

bamApplyBQSR = bamMDToJoin.join(recalTable, by:[0,1])

if (step == 'recalibrate') bamApplyBQSR = inputSample

bamApplyBQSR = bamApplyBQSR.dump(tag:'BAM + BAI + RECAL TABLE')
// [DUMP: recal.table] ['normal', 'normal', normal.md.bam, normal.md.bai, normal.recal.table]

bamApplyBQSR = bamApplyBQSR.combine(intApplyBQSR)

bamApplyBQSR = bamApplyBQSR.dump(tag:'BAM + BAI + RECAL TABLE + INT')
// [DUMP: BAM + BAI + RECAL TABLE + INT] ['normal', 'normal', normal.md.bam, normal.md.bai, normal.recal.table, 1_1-200000.bed]

// STEP 4: RECALIBRATING

process ApplyBQSR {
    label 'gatk'
    label 'memory_singleCPU_2_task'
    label 'cpus_2'

    tag {idPatient + "-" + idSample + "-" + intervalBed.baseName}

    input:
        set idPatient, idSample, file(bam), file(bai), file(recalibrationReport), file(intervalBed) from bamApplyBQSR
        file(dict) from ch_dict
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fastaFai

    output:
        set idPatient, idSample, file("${prefix}${idSample}.recal.bam") into bamMergeBamRecal

    script:
    prefix = params.no_intervals ? "" : "${intervalBed.baseName}_"
    intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        ApplyBQSR \
        -R ${fasta} \
        --input ${bam} \
        --output ${prefix}${idSample}.recal.bam \
        ${intervalsOptions} \
        --bqsr-recal-file ${recalibrationReport}
    """
}

bamMergeBamRecal = bamMergeBamRecal.groupTuple(by:[0, 1])
(bamMergeBamRecal, bamMergeBamRecalNoInt) = bamMergeBamRecal.into(2)

// EP 4.5: MERGING THE RECALIBRATED BAM FILES

process MergeBamRecal { 
    label 'samtools'
    label 'cpus_8'

    tag {idPatient + "-" + idSample}

    publishDir "${params.outputDir}/Preprocessing/${idSample}/Recalibrated", mode: params.publishDirMode

    input:
        set idPatient, idSample, file(bam) from bamMergeBamRecal

    output:
        set idPatient, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bam.bai") into bamRecal
        set idPatient, idSample, file("${idSample}.recal.bam") into bamRecalQC
        set idPatient, idSample into bamRecalTSV

    when: !(params.no_intervals)

    script:
    """
    samtools merge --threads ${task.cpus} ${idSample}.recal.bam ${bam}
    samtools index ${idSample}.recal.bam
    """
}

// STEP 4.5': INDEXING THE RECALIBRATED BAM FILES

process IndexBamRecal {
    label 'samtools'
    label 'cpus_8'

    tag {idPatient + "-" + idSample}

    publishDir "${params.outputDir}/Preprocessing/${idSample}/Recalibrated", mode: params.publishDirMode

    input:
        set idPatient, idSample, file("${idSample}.recal.bam") from bamMergeBamRecalNoInt

    output:
        set idPatient, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bam.bai") into bamRecalNoInt
        set idPatient, idSample, file("${idSample}.recal.bam") into bamRecalQCnoInt
        set idPatient, idSample into bamRecalTSVnoInt

    when: params.no_intervals

    script:
    """
    samtools index ${idSample}.recal.bam
    """
}

bamRecal = bamRecal.mix(bamRecalNoInt)
bamRecalQC = bamRecalQC.mix(bamRecalQCnoInt)
bamRecalTSV = bamRecalTSV.mix(bamRecalTSVnoInt)

(bamRecalBamQC, bamRecalSamToolsStats) = bamRecalQC.into(2)
(bamRecalTSV, bamRecalSampleTSV) = bamRecalTSV.into(2)

// Creating a TSV file to restart from this step
bamRecalTSV.map { idPatient, idSample ->
    gender = genderMap[idPatient]
    status = statusMap[idPatient, idSample]
    bam = "${params.outputDir}/Preprocessing/${idSample}/Recalibrated/${idSample}.recal.bam"
    bai = "${params.outputDir}/Preprocessing/${idSample}/Recalibrated/${idSample}.recal.bam.bai"
    "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"
}.collectFile(
    name: 'recalibrated.tsv', sort: true, storeDir: "${params.outputDir}/Preprocessing/TSV"
)

bamRecalSampleTSV
    .collectFile(storeDir: "${params.outputDir}/Preprocessing/TSV") {
        idPatient, idSample ->
        status = statusMap[idPatient, idSample]
        gender = genderMap[idPatient]
        bam = "${params.outputDir}/Preprocessing/${idSample}/Recalibrated/${idSample}.recal.bam"
        bai = "${params.outputDir}/Preprocessing/${idSample}/Recalibrated/${idSample}.recal.bam.bai"
        ["recalibrated_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"]
}

// STEP 5: QC

process SamtoolsStats {
    label 'samtools'
    label 'cpus_2'

    tag {idPatient + "-" + idSample}

    publishDir "${params.outputDir}/Reports/${idSample}/SamToolsStats", mode: params.publishDirMode

    input:
        set idPatient, idSample, file(bam) from bamRecalSamToolsStats

    output:
        file ("${bam}.samtools.stats.out") into samtoolsStatsReport

    when: !('samtoolsStats' in skipQC)

    script:
    """
    samtools stats ${bam} > ${bam}.samtools.stats.out
    """
}

samtoolsStatsReport = samtoolsStatsReport.dump(tag:'SAMTools')

// bamBamQC = bamMappedBamQC // Mapreads only
bamBamQC = bamMappedBamQC.mix(bamRecalBamQC) // Mapreads + MapQ + MarkDuplicates + ApplyBQSR 

process BamQC {
    label 'qualimap'
    label 'memory_max'
    label 'cpus_16'

    tag {idPatient + "-" + idSample}

    publishDir "${params.outputDir}/Reports/${idSample}/bamQC", mode: params.publishDirMode

    input:
        set idPatient, idSample, file(bam) from bamBamQC
        file(targetBED) from ch_targetBED

    output:
        file("${bam.baseName}") into bamQCReport

    when: !('bamqc' in skipQC)

    script:
    use_bed = params.targetBED ? "-gff ${targetBED}" : ''
    """
    qualimap --java-mem-size=${task.memory.toGiga()}G \
        bamqc \
        -bam ${bam} \
        --paint-chromosome-limits \
        --genome-gc-distr HUMAN \
        $use_bed \
        -nt ${task.cpus} \
        -skip-duplicated \
        --skip-dup-mode 0 \
        -outdir ${bam.baseName} \
        -outformat HTML
    """
}

bamQCReport = bamQCReport.dump(tag:'BamQC')


/*
================================================================================
                            GERMLINE VARIANT CALLING
================================================================================
*/

// When no knownIndels for mapping, Channel bamRecal is indexedBam
bamRecal = (params.knownIndels && step == 'mapping') ? bamRecal : indexedBam

// When starting with variant calling, Channel bamRecal is inputSample
if (step == 'variantcalling') bamRecal = inputSample

bamRecal = bamRecal.dump(tag:'BAM')

// Here we have a recalibrated bam set
// The TSV file is formatted like: "idPatient status idSample bamFile baiFile"
// Manta will be run in Germline mode, or in Tumor mode depending on status
// HaplotypeCaller, TIDDIT and Strelka will be run for Normal and Tumor samples

(bamMantaSingle, bamStrelkaSingle, bamTIDDIT, bamRecalAll, bamRecalAllTemp) = bamRecal.into(5)

// To speed Variant Callers up we are chopping the reference into smaller pieces
// Do variant calling by this intervals, and re-merge the VCFs

bamHaplotypeCaller = bamRecalAllTemp.combine(intHaplotypeCaller)

// STEP GATK HAPLOTYPECALLER.1

process HaplotypeCaller {
    label 'gatk'
    label 'memory_singleCPU_task_sq'
    label 'cpus_2'

    tag {idSample + "-" + intervalBed.baseName}

    input:
        set idPatient, idSample, file(bam), file(bai), file(intervalBed) from bamHaplotypeCaller
        file(dbsnp) from ch_dbsnp
        file(dbsnpIndex) from ch_dbsnpIndex
        file(dict) from ch_dict
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fastaFai

    output:
        set val("HaplotypeCallerGVCF"), idPatient, idSample, file("${intervalBed.baseName}_${idSample}.g.vcf") into gvcfHaplotypeCaller
        set idPatient, idSample, file(intervalBed), file("${intervalBed.baseName}_${idSample}.g.vcf") into gvcfGenotypeGVCFs

    when: 'haplotypecaller' in tools

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g -Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        HaplotypeCaller \
        -R ${fasta} \
        -I ${bam} \
        -L ${intervalBed} \
        -D ${dbsnp} \
        -O ${intervalBed.baseName}_${idSample}.g.vcf \
        -ERC GVCF
    """
}

gvcfHaplotypeCaller = gvcfHaplotypeCaller.groupTuple(by:[0, 1, 2])

if (params.noGVCF) gvcfHaplotypeCaller.close()
else gvcfHaplotypeCaller = gvcfHaplotypeCaller.dump(tag:'GVCF HaplotypeCaller')

// STEP GATK HAPLOTYPECALLER.2

process GenotypeGVCFs {
    label 'gatk'
    tag {idSample + "-" + intervalBed.baseName}

    input:
        set idPatient, idSample, file(intervalBed), file(gvcf) from gvcfGenotypeGVCFs
        file(dbsnp) from ch_dbsnp
        file(dbsnpIndex) from ch_dbsnpIndex
        file(dict) from ch_dict
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fastaFai

    output:
    set val("HaplotypeCaller"), idPatient, idSample, file("${intervalBed.baseName}_${idSample}.vcf") into vcfGenotypeGVCFs

    when: 'haplotypecaller' in tools

    script:
    // Using -L is important for speed and we have to index the interval files also
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        IndexFeatureFile -F ${gvcf}

    gatk --java-options -Xmx${task.memory.toGiga()}g \
        GenotypeGVCFs \
        -R ${fasta} \
        -L ${intervalBed} \
        -D ${dbsnp} \
        -V ${gvcf} \
        -O ${intervalBed.baseName}_${idSample}.vcf
    """
}
vcfGenotypeGVCFs = vcfGenotypeGVCFs.groupTuple(by:[0, 1, 2])

// STEP MANTA.1 - SINGLE MODE

process MantaSingle {
    label 'manta'
    label 'cpus_max'
    label 'memory_max'

    tag {idSample}

    publishDir "${params.outputDir}/VariantCalling/${idSample}/Manta", mode: params.publishDirMode

    input:
        set idPatient, idSample, file(bam), file(bai) from bamMantaSingle
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fastaFai
        file(targetBED) from ch_targetBED

    output:
        set val("Manta"), idPatient, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") into vcfMantaSingle

    when: 'manta' in tools

    script:
    beforeScript = params.targetBED ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
    options = params.targetBED ? "--exome --callRegions call_targets.bed.gz" : ""
    status = statusMap[idPatient, idSample]
    inputbam = status == 0 ? "--bam" : "--tumorBam"
    vcftype = status == 0 ? "diploid" : "tumor"
    """
    ${beforeScript}
    configManta.py \
        ${inputbam} ${bam} \
        --reference ${fasta} \
        ${options} \
        --runDir Manta

    python Manta/runWorkflow.py -m local -j ${task.cpus}

    mv Manta/results/variants/candidateSmallIndels.vcf.gz \
        Manta_${idSample}.candidateSmallIndels.vcf.gz
    mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
        Manta_${idSample}.candidateSmallIndels.vcf.gz.tbi
    mv Manta/results/variants/candidateSV.vcf.gz \
        Manta_${idSample}.candidateSV.vcf.gz
    mv Manta/results/variants/candidateSV.vcf.gz.tbi \
        Manta_${idSample}.candidateSV.vcf.gz.tbi
    mv Manta/results/variants/${vcftype}SV.vcf.gz \
        Manta_${idSample}.${vcftype}SV.vcf.gz
    mv Manta/results/variants/${vcftype}SV.vcf.gz.tbi \
        Manta_${idSample}.${vcftype}SV.vcf.gz.tbi
    """
}

vcfMantaSingle = vcfMantaSingle.dump(tag:'Single Manta')

/*
================================================================================
                             SOMATIC VARIANT CALLING
================================================================================
*/
// Ascat
(bamAscat, bamMpileup, bamMpileupNoInt, bamRecalAll) = bamRecalAll.into(4)

// separate BAM by status
bamNormal = Channel.create()
bamTumor = Channel.create()

bamRecalAll
    .choice(bamTumor, bamNormal) {statusMap[it[0], it[1]] == 0 ? 1 : 0}

// Crossing Normal and Tumor to get a T/N pair for Somatic Variant Calling
// Remapping channel to remove common key idPatient
pairBam = bamNormal.cross(bamTumor).map {
    normal, tumor ->
    [normal[0], normal[1], normal[2], normal[3], tumor[1], tumor[2], tumor[3]]
}

pairBam = pairBam.dump(tag:'BAM Somatic Pair')

// Manta,  Mutect2
(pairBamManta, pairBamStrelka, pairBamStrelkaBP, pairBamCalculateContamination, pairBamFilterMutect2, pairBamTNscope, pairBam) = pairBam.into(7)

intervalPairBam = pairBam.spread(bedIntervals)

bamMpileup = bamMpileup.spread(intMpileup)


// intervals for Mutect2 calls, FreeBayes and pileups for Mutect2 filtering
(pairBamMutect2, pairBamFreeBayes, pairBamPileupSummaries) = intervalPairBam.into(3)


// STEP MANTA.2 - SOMATIC PAIR

process Manta {
    label 'manta'
    label 'cpus_max'
    label 'memory_max'

    tag {idSampleTumor + "_vs_" + idSampleNormal}

    publishDir "${params.outputDir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/Manta", mode: params.publishDirMode

    input:
        set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from pairBamManta
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fastaFai
        file(targetBED) from ch_targetBED

    output:
        set val("Manta"), idPatient, val("${idSampleTumor}_vs_${idSampleNormal}"), file("*.vcf.gz"), file("*.vcf.gz.tbi") into vcfManta
        set idPatient, idSampleNormal, idSampleTumor, file("*.candidateSmallIndels.vcf.gz"), file("*.candidateSmallIndels.vcf.gz.tbi") into mantaToStrelka

    when: 'manta' in tools

    script:
    beforeScript = params.targetBED ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
    options = params.targetBED ? "--exome --callRegions call_targets.bed.gz" : ""
    """
    ${beforeScript}
    configManta.py \
        --normalBam ${bamNormal} \
        --tumorBam ${bamTumor} \
        --reference ${fasta} \
        ${options} \
        --runDir Manta

    python Manta/runWorkflow.py -m local -j ${task.cpus}

    mv Manta/results/variants/candidateSmallIndels.vcf.gz \
        Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSmallIndels.vcf.gz
    mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
        Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSmallIndels.vcf.gz.tbi
    mv Manta/results/variants/candidateSV.vcf.gz \
        Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSV.vcf.gz
    mv Manta/results/variants/candidateSV.vcf.gz.tbi \
        Manta_${idSampleTumor}_vs_${idSampleNormal}.candidateSV.vcf.gz.tbi
    mv Manta/results/variants/diploidSV.vcf.gz \
        Manta_${idSampleTumor}_vs_${idSampleNormal}.diploidSV.vcf.gz
    mv Manta/results/variants/diploidSV.vcf.gz.tbi \
        Manta_${idSampleTumor}_vs_${idSampleNormal}.diploidSV.vcf.gz.tbi
    mv Manta/results/variants/somaticSV.vcf.gz \
        Manta_${idSampleTumor}_vs_${idSampleNormal}.somaticSV.vcf.gz
    mv Manta/results/variants/somaticSV.vcf.gz.tbi \
        Manta_${idSampleTumor}_vs_${idSampleNormal}.somaticSV.vcf.gz.tbi
    """
}

vcfManta = vcfManta.dump(tag:'Manta')

// STEP GATK MUTECT2.1 - RAW CALLS

process Mutect2 {
    tag {idSampleTumor + "_vs_" + idSampleNormal + "-" + intervalBed.baseName}
    label 'gatk'
    label 'cpus_1'

    input:
        set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), file(intervalBed) from pairBamMutect2
        file(dict) from ch_dict
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fastaFai
        file(germlineResource) from ch_germlineResource
        file(germlineResourceIndex) from ch_germlineResourceIndex
        file(intervals) from ch_intervals

    output:
        set val("Mutect2"), 
            idPatient,
            val("${idSampleTumor}_vs_${idSampleNormal}"),
            file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf") into mutect2Output
        set idPatient,
            idSampleTumor,
            idSampleNormal,
            file("${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf.stats") optional true into mutect2Stats

    when: 'mutect2' in tools

    script:
    // please make a panel-of-normals, using at least 40 samples
    // https://gatkforums.broadinstitute.org/gatk/discussion/11136/how-to-call-somatic-mutations-using-gatk4-mutect2
    PON = params.pon ? "--panel-of-normals ${pon}" : ""
    """
    # Get raw calls
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
      Mutect2 \
      -R ${fasta}\
      -I ${bamTumor}  -tumor ${idSampleTumor} \
      -I ${bamNormal} -normal ${idSampleNormal} \
      -L ${intervalBed} \
      --germline-resource ${germlineResource} \
      ${PON} \
      -O ${intervalBed.baseName}_${idSampleTumor}_vs_${idSampleNormal}.vcf
    """
}

mutect2Output = mutect2Output.groupTuple(by:[0,1,2])
(mutect2Output, mutect2OutForStats) = mutect2Output.into(2)

(mutect2Stats, intervalStatsFiles) = mutect2Stats.into(2)
mutect2Stats = mutect2Stats.groupTuple(by:[0,1,2])

// STEP GATK MUTECT2.2 - MERGING STATS

process MergeMutect2Stats {
    tag {idSampleTumor + "_vs_" + idSampleNormal}
    label 'gatk'

    publishDir "${params.outputDir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/Mutect2", mode: params.publishDirMode

    input:
        set caller, idPatient, idSampleTumor_vs_idSampleNormal, file(vcfFiles) from mutect2OutForStats // corresponding small VCF chunks
        set idPatient, idSampleTumor, idSampleNormal, file(statsFiles) from mutect2Stats               // the actual stats files
        file(dict) from ch_dict
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fastaFai
        file(germlineResource) from ch_germlineResource
        file(germlineResourceIndex) from ch_germlineResourceIndex
        file(intervals) from ch_intervals

    output:
        file("${idSampleTumor_vs_idSampleNormal}.vcf.gz.stats") into mergedStatsFile

    when: 'mutect2' in tools

    script:     
      stats = statsFiles.collect{ "-stats ${it} " }.join(' ')
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        MergeMutectStats \
        ${stats} \
        -O ${idSampleTumor}_vs_${idSampleNormal}.vcf.gz.stats
    """
}

// we are merging the VCFs that are called separatelly for different intervals
// so we can have a single sorted VCF containing all the calls for a given caller

// STEP MERGING VCF - FREEBAYES, GATK HAPLOTYPECALLER & GATK MUTECT2 (UNFILTERED)

//vcfConcatenateVCFs = mutect2Output.mix(vcfFreeBayes, vcfGenotypeGVCFs, gvcfHaplotypeCaller)
vcfConcatenateVCFs = mutect2Output.mix(vcfGenotypeGVCFs, gvcfHaplotypeCaller)
vcfConcatenateVCFs = vcfConcatenateVCFs.dump(tag:'VCF to merge')

process ConcatVCF {
    label 'bcftools'
    label 'cpus_8'

    tag {variantCaller + "-" + idSample}

    publishDir "${params.outputDir}/VariantCalling/${idSample}/${"$variantCaller"}", mode: params.publishDirMode

    input:
        set variantCaller, idPatient, idSample, file(vcFiles) from vcfConcatenateVCFs
        file(fastaFai) from ch_fastaFai
        file(targetBED) from ch_targetBED

    output:
    // we have this funny *_* pattern to avoid copying the raw calls to publishdir
        set variantCaller, idPatient, idSample, file("*_*.vcf.gz"), file("*_*.vcf.gz.tbi") into vcfConcatenated

    when: ('haplotypecaller' in tools || 'mutect2' in tools || 'freebayes' in tools)

    script:
    if (variantCaller == 'HaplotypeCallerGVCF') 
      outputFile = "HaplotypeCaller_${idSample}.g.vcf"
    else if (variantCaller == "Mutect2") 
      outputFile = "unfiltered_${variantCaller}_${idSample}.vcf"
    else 
      outputFile = "${variantCaller}_${idSample}.vcf"
    options = params.targetBED ? "-t ${targetBED}" : ""
    """
    concatenateVCFs.sh -i ${fastaFai} -c ${task.cpus} -o ${outputFile} ${options}
    """
}

(vcfConcatenated, vcfConcatenatedForFilter) = vcfConcatenated.into(2)
vcfConcatenated = vcfConcatenated.dump(tag:'VCF')

// STEP GATK MUTECT2.3 - GENERATING PILEUP SUMMARIES

process PileupSummariesForMutect2 {
    tag {idSampleTumor + "_vs_" + idSampleNormal + "_" + intervalBed.baseName }
    label 'gatk'
    label 'cpus_1'

    input:
        set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor), file(intervalBed) from pairBamPileupSummaries 
        set idPatient, idSampleNormal, idSampleTumor, file(statsFile) from intervalStatsFiles
        file(germlineResource) from ch_germlineResource
        file(germlineResourceIndex) from ch_germlineResourceIndex

    output:
        set idPatient,
            idSampleTumor,
            file("${intervalBed.baseName}_${idSampleTumor}_pileupsummaries.table") into pileupSummaries

    when: 'mutect2' in tools

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        GetPileupSummaries \
        -I ${bamTumor} \
        -V ${germlineResource} \
        -L ${intervalBed} \
        -O ${intervalBed.baseName}_${idSampleTumor}_pileupsummaries.table
    """
}

pileupSummaries = pileupSummaries.groupTuple(by:[0,1])

// STEP GATK MUTECT2.4 - MERGING PILEUP SUMMARIES

process MergePileupSummaries {
    label 'gatk'
    label 'cpus_1'

    tag {idPatient + "_" + idSampleTumor}

    publishDir "${params.outputDir}/VariantCalling/${idSampleTumor}/Mutect2", mode: params.publishDirMode

    input:
        set idPatient, idSampleTumor, file(pileupSums) from pileupSummaries
        file(dict) from ch_dict

    output:
        file("${idSampleTumor}_pileupsummaries.table.tsv") into mergedPileupFile

    when: 'mutect2' in tools
    script:
        allPileups = pileupSums.collect{ "-I ${it} " }.join(' ')
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        GatherPileupSummaries \
        --sequence-dictionary ${dict} \
        ${allPileups} \
        -O ${idSampleTumor}_pileupsummaries.table.tsv
    """
}

// STEP GATK MUTECT2.5 - CALCULATING CONTAMINATION

process CalculateContamination {
    label 'gatk'
    label 'cpus_1'

    tag {idSampleTumor + "_vs_" + idSampleNormal}

    publishDir "${params.outputDir}/VariantCalling/${idSampleTumor}/Mutect2", mode: params.publishDirMode

    input:
        set idPatient, idSampleNormal, file(bamNormal), file(baiNormal), idSampleTumor, file(bamTumor), file(baiTumor) from pairBamCalculateContamination 
        file("${idSampleTumor}_pileupsummaries.table") from mergedPileupFile
  
    output:
        file("${idSampleTumor}_contamination.table") into contaminationTable

    when: 'mutect2' in tools

    script:     
    """
    # calculate contamination
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        CalculateContamination \
        -I ${idSampleTumor}_pileupsummaries.table \
        -O ${idSampleTumor}_contamination.table
    """
}

// STEP GATK MUTECT2.6 - FILTERING CALLS

process FilterMutect2Calls {
    label 'gatk'
    label 'cpus_1'

    tag {idSampleTN}

    publishDir "${params.outputDir}/VariantCalling/${idSampleTN}/${"$variantCaller"}", mode: params.publishDirMode

    input:
        set variantCaller, idPatient, idSampleTN, file(unfiltered), file(unfilteredIndex) from vcfConcatenatedForFilter
        file("${idSampleTN}.vcf.gz.stats") from mergedStatsFile
        file("${idSampleTN}_contamination.table") from contaminationTable
        file(dict) from ch_dict
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fastaFai
        file(germlineResource) from ch_germlineResource
        file(germlineResourceIndex) from ch_germlineResourceIndex
        file(intervals) from ch_intervals
        
    output:
        set val("Mutect2"), idPatient, idSampleTN,
            file("filtered_${variantCaller}_${idSampleTN}.vcf.gz"),
            file("filtered_${variantCaller}_${idSampleTN}.vcf.gz.tbi"),
            file("filtered_${variantCaller}_${idSampleTN}.vcf.gz.filteringStats.tsv") into filteredMutect2Output

    when: 'mutect2' in tools

    script:
    """
    # do the actual filtering
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        FilterMutectCalls \
        -V ${unfiltered} \
        --contamination-table ${idSampleTN}_contamination.table \
        --stats ${idSampleTN}.vcf.gz.stats \
        -R ${fasta} \
        -O filtered_${variantCaller}_${idSampleTN}.vcf.gz
    """
}

// STEP ASCAT.1 - ALLELECOUNTER

// Run commands and code from Malin Larsson
// Based on Jesper Eisfeldt's code
process AlleleCounter {
    label 'canceritAllelecount'
    label 'memory_singleCPU_2_task'

    tag {idSample}

    input:
        set idPatient, idSample, file(bam), file(bai) from bamAscat
        file(acLoci) from ch_acLoci
        file(dict) from ch_dict
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fastaFai

    output:
        set idPatient, idSample, file("${idSample}.alleleCount") into alleleCounterOut

    when: 'ascat' in tools

    script:
    """
    alleleCounter \
        -l ${acLoci} \
        -r ${fasta} \
        -b ${bam} \
        -o ${idSample}.alleleCount;
    """
}

alleleCountOutNormal = Channel.create()
alleleCountOutTumor = Channel.create()

alleleCounterOut
    .choice(alleleCountOutTumor, alleleCountOutNormal) {statusMap[it[0], it[1]] == 0 ? 1 : 0}

alleleCounterOut = alleleCountOutNormal.combine(alleleCountOutTumor)

alleleCounterOut = alleleCounterOut.map {
    idPatientNormal, idSampleNormal, alleleCountOutNormal,
    idPatientTumor, idSampleTumor, alleleCountOutTumor ->
    [idPatientNormal, idSampleNormal, idSampleTumor, alleleCountOutNormal, alleleCountOutTumor]
}
// STEP ASCAT.2 - CONVERTALLELECOUNTS

// R script from Malin Larssons bitbucket repo:
// https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
process ConvertAlleleCounts {
    label 'ascat'
    label 'memory_singleCPU_2_task'

    tag {idSampleTumor + "_vs_" + idSampleNormal}

    publishDir "${params.outputDir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/ASCAT", mode: params.publishDirMode

    input:
        set idPatient, idSampleNormal, idSampleTumor, file(alleleCountNormal), file(alleleCountTumor) from alleleCounterOut

    output:
        set idPatient, idSampleNormal, idSampleTumor, file("${idSampleNormal}.BAF"), file("${idSampleNormal}.LogR"), file("${idSampleTumor}.BAF"), file("${idSampleTumor}.LogR") into convertAlleleCountsOut

    when: 'ascat' in tools

    script:
    gender = genderMap[idPatient]
    """
    Rscript ${workflow.projectDir}/bin/convertAlleleCounts.r ${idSampleTumor} ${alleleCountTumor} ${idSampleNormal} ${alleleCountNormal} ${gender}
    """
}

// STEP ASCAT.3 - ASCAT

// R scripts from Malin Larssons bitbucket repo:
// https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
process Ascat {
    label 'ascat'
    label 'memory_singleCPU_2_task'

    tag {idSampleTumor + "_vs_" + idSampleNormal}

    publishDir "${params.outputDir}/VariantCalling/${idSampleTumor}_vs_${idSampleNormal}/ASCAT", mode: params.publishDirMode

    input:
        set idPatient, idSampleNormal, idSampleTumor, file(bafNormal), file(logrNormal), file(bafTumor), file(logrTumor) from convertAlleleCountsOut
        file(acLociGC) from ch_acLociGC

    output:
        set val("ASCAT"), idPatient, idSampleNormal, idSampleTumor, file("${idSampleTumor}.*.{png,txt}") into ascatOut

    when: 'ascat' in tools

    script:
    """
    # get rid of "chr" string if there is any
    for f in *BAF *LogR; do sed 's/chr//g' \$f > tmpFile; mv tmpFile \$f;done
    Rscript ${workflow.projectDir}/bin/run_ascat.r ${bafTumor} ${logrTumor} ${bafNormal} ${logrNormal} ${idSampleTumor} ${baseDir} ${acLociGC}
    """
}

ascatOut.dump(tag:'ASCAT')

// Remapping channels for QC and annotation

(vcfMantaSomaticSV, vcfMantaDiploidSV) = vcfManta.into(2)

vcfKeep = Channel.empty().mix(
    vcfMantaSingle.map {
        variantcaller, idPatient, idSample, vcf, tbi ->
        [variantcaller, idSample, vcf[2]]
    },
    vcfMantaDiploidSV.map {
        variantcaller, idPatient, idSample, vcf, tbi ->
        [variantcaller, idSample, vcf[2]]
    },
    vcfMantaSomaticSV.map {
        variantcaller, idPatient, idSample, vcf, tbi ->
        [variantcaller, idSample, vcf[3]]
    })


(vcfBCFtools, vcfVCFtools, vcfAnnotation) = vcfKeep.into(3)

/*
================================================================================
                                   ANNOTATION
================================================================================
*/

if (step == 'annotate') {
    vcfToAnnotate = Channel.create()
    vcfNoAnnotate = Channel.create()

    if (tsvPath == []) {
    // Sarek, by default, annotates all available vcfs that it can find in the VariantCalling directory
    // Excluding vcfs from FreeBayes, and g.vcf from HaplotypeCaller
    // Basically it's: results/VariantCalling/*/{HaplotypeCaller,Manta,Mutect2,SentieonDNAseq,SentieonDNAscope,SentieonTNscope,Strelka,TIDDIT}/*.vcf.gz
    // Without *SmallIndels.vcf.gz from Manta, and *.genome.vcf.gz from Strelka
    // The small snippet `vcf.minus(vcf.fileName)[-2]` catches idSample
    // This field is used to output final annotated VCFs in the correct directory
      Channel.empty().mix(
        Channel.fromPath("${params.outputDir}/VariantCalling/*/HaplotypeCaller/*.vcf.gz")
          .flatten().map{vcf -> ['HaplotypeCaller', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outputDir}/VariantCalling/*/Manta/*[!candidate]SV.vcf.gz")
          .flatten().map{vcf -> ['Manta', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outputDir}/VariantCalling/*/Mutect2/*.vcf.gz")
          .flatten().map{vcf -> ['Mutect2', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outputDir}/VariantCalling/*/SentieonDNAseq/*.vcf.gz")
          .flatten().map{vcf -> ['SentieonDNAseq', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outputDir}/VariantCalling/*/SentieonDNAscope/*.vcf.gz")
          .flatten().map{vcf -> ['SentieonDNAscope', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outputDir}/VariantCalling/*/SentieonTNscope/*.vcf.gz")
          .flatten().map{vcf -> ['SentieonTNscope', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outputDir}/VariantCalling/*/Strelka/*{somatic,variant}*.vcf.gz")
          .flatten().map{vcf -> ['Strelka', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outputDir}/VariantCalling/*/TIDDIT/*.vcf.gz")
          .flatten().map{vcf -> ['TIDDIT', vcf.minus(vcf.fileName)[-2].toString(), vcf]}
      ).choice(vcfToAnnotate, vcfNoAnnotate) {
        annotateTools == [] || (annotateTools != [] && it[0] in annotateTools) ? 0 : 1
      }
    } else if (annotateTools == []) {
    // Annotate user-submitted VCFs
    // If user-submitted, Sarek assume that the idSample should be assumed automatically
      vcfToAnnotate = Channel.fromPath(tsvPath)
        .map{vcf -> ['userspecified', vcf.minus(vcf.fileName)[-2].toString(), vcf]}
    } else exit 1, "specify only tools or files to annotate, not both"

    vcfNoAnnotate.close()
    vcfAnnotation = vcfAnnotation.mix(vcfToAnnotate)
}
// as now have the list of VCFs to annotate, the first step is to annotate with allele frequencies, if there are any

(vcfSnpeff, vcfVep) = vcfAnnotation.into(2)

vcfVep = vcfVep.map {
  variantCaller, idSample, vcf ->
  [variantCaller, idSample, vcf, null]
}

// STEP SNPEFF

process Snpeff {
    tag {"${idSample} - ${variantCaller} - ${vcf}"}

    publishDir params.outputDir, mode: params.publishDirMode, saveAs: {
        if (it == "${reducedVCF}_snpEff.ann.vcf") null
        else "Reports/${idSample}/snpEff/${it}"
    }

    input:
        set variantCaller, idSample, file(vcf) from vcfSnpeff
        file(dataDir) from ch_snpEff_cache
        val snpeffDb from ch_snpeffDb

    output:
        set file("${reducedVCF}_snpEff.txt"), file("${reducedVCF}_snpEff.html"), file("${reducedVCF}_snpEff.csv") into snpeffReport
        set variantCaller, idSample, file("${reducedVCF}_snpEff.ann.vcf") into snpeffVCF

    when: 'snpeff' in tools || 'merge' in tools

    script:
    reducedVCF = reduceVCF(vcf.fileName)
    cache = (params.snpEff_cache && params.annotation_cache) ? "-dataDir \${PWD}/${dataDir}" : ""
    """
    snpEff -Xmx${task.memory.toGiga()}g \
        ${snpeffDb} \
        -csvStats ${reducedVCF}_snpEff.csv \
        -nodownload \
        ${cache} \
        -canon \
        -v \
        ${vcf} \
        > ${reducedVCF}_snpEff.ann.vcf

    mv snpEff_summary.html ${reducedVCF}_snpEff.html
    mv ${reducedVCF}_snpEff.genes.txt ${reducedVCF}_snpEff.txt
    """
}

snpeffReport = snpeffReport.dump(tag:'snpEff report')

// STEP COMPRESS AND INDEX VCF.1 - SNPEFF

process CompressVCFsnpEff {
    tag {"${idSample} - ${vcf}"}

    publishDir "${params.outputDir}/Annotation/${idSample}/snpEff", mode: params.publishDirMode

    input:
        set variantCaller, idSample, file(vcf) from snpeffVCF

    output:
        set variantCaller, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") into (compressVCFsnpEffOut)

    script:
    """
    bgzip < ${vcf} > ${vcf}.gz
    tabix ${vcf}.gz
    """
}

compressVCFsnpEffOut = compressVCFsnpEffOut.dump(tag:'VCF')






/*
================================================================================
                                     MultiQC
================================================================================
*/

// STEP MULTIQC

process MultiQC { 
    label 'multiqc' 
    publishDir "${params.outputDir}/Reports/MultiQC", mode: params.publishDirMode

    input:
        file (multiqcConfig) from Channel.value(params.multiqc_config ? file(params.multiqc_config) : "")
        file (versions) from yamlSoftwareVersion
        file ('bamQC/*') from bamQCReport.collect().ifEmpty([])
        file ('FastQC/*') from fastQCReport.collect().ifEmpty([])
        file ('MarkDuplicates/*') from markDuplicatesReport.collect().ifEmpty([])
        file ('SamToolsStats/*') from samtoolsStatsReport.collect().ifEmpty([])
        file ('snpEff/*') from snpeffReport.collect().ifEmpty([])

    output:
        set file("*multiqc_report.html"), file("*multiqc_data") into multiQCOut

    when: !('multiqc' in skipQC)

    script:
    """
    multiqc -f -v .
    """
}

multiQCOut.dump(tag:'MultiQC')

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[EUCANCAN/nf-wgswes] Successful: ${workflow.runName}"
    if (!workflow.success) subject = "[EUCANCAN/nf-wgswes] FAILED: ${workflow.runName}"
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    if (workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[EUCANCAN/nf-wgswes] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[EUCANCAN/nf-wgswes] Could not attach MultiQC report to summary email"
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[EUCANCAN/nf-wgswes] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            [ 'mail', '-s', subject, params.email ].execute() << email_txt
            log.info "[EUCANCAN/nf-wgswes] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outputDir}/pipeline_info/")
    if (!output_d.exists()) output_d.mkdirs()
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_reset  = params.monochrome_logs ? '' : "\033[0m";
    c_red    = params.monochrome_logs ? '' : "\033[0;31m";
    c_green  = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "${c_purple}Warning, pipeline completed, but with errored process(es)${c_reset}"
        log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCountFmt}${c_reset}"
        log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCountFmt}${c_reset}"
    }

    if (workflow.success) log.info "${c_purple}[EUCANCAN/nf-wgswes]${c_green} Pipeline completed successfully${c_reset}"
    else {
        checkHostname()
        log.info "${c_purple}[EUCANCAN/nf-wgswes]${c_red} Pipeline completed with errors${c_reset}"
    }
}

/*
================================================================================
                                functions
================================================================================
*/

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-sarek-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'EUCANCAN/nf-wgswes Workflow Summary'
    section_href: 'https://github.com/EUCANCAN/nf-wgswes'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k, v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}

def nfcoreHeader() {
    // Log colors ANSI codes
    c_reset  = params.monochrome_logs ? '' : "\033[0m";
    c_dim    = params.monochrome_logs ? '' : "\033[2m";
    c_black  = params.monochrome_logs ? '' : "\033[0;30m";
    c_red    = params.monochrome_logs ? '' : "\033[0;31m";
    c_green  = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue   = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan   = params.monochrome_logs ? '' : "\033[0;36m";
    c_white  = params.monochrome_logs ? '' : "\033[0;37m";

    return """    ${c_dim}----------------------------------------------------${c_reset}
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
        ${c_white}____${c_reset}
      ${c_white}.´ _  `.${c_reset}
     ${c_white}/  ${c_green}|\\${c_reset}`-_ \\${c_reset}     ${c_blue} __        __   ___     ${c_reset}
    ${c_white}|   ${c_green}| \\${c_reset}  `-|${c_reset}    ${c_blue}|__`  /\\  |__) |__  |__/${c_reset}
     ${c_white}\\ ${c_green}|   \\${c_reset}  /${c_reset}     ${c_blue}.__| /¯¯\\ |  \\ |___ |  \\${c_reset}
      ${c_white}`${c_green}|${c_reset}____${c_green}\\${c_reset}´${c_reset}

    ${c_purple}  EUCANCAN/nf-wgswes v${workflow.manifest.version}${c_reset}
    ${c_dim}----------------------------------------------------${c_reset}
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}

// Check if a row has the expected number of item
def checkNumberOfItem(row, number) {
    if (row.size() != number) exit 1, "Malformed row in TSV file: ${row}, see --help for more information"
    return true
}

// Check parameter existence
def checkParameterExistence(it, list) {
    if (!list.contains(it)) {
        log.warn "Unknown parameter: ${it}"
        return false
    }
    return true
}

// Compare each parameter with a list of parameters
def checkParameterList(list, realList) {
    return list.every{ checkParameterExistence(it, realList) }
}

// Check if params.item exists and return params.genomes[params.genome].item otherwise
def checkParamReturnFile(item) {
    // Handle deprecation
    if (params.genomeDict && item == "dict") return file(params.genomeDict)
    if (params.genomeFile && item == "fasta") return file(params.genomeFile)
    if (params.genomeIndex && item == "fastaFai") return file(params.genomeIndex)

    params."${item}" = params.genomes[params.genome]."${item}"
    return file(params."${item}")
}

// Define list of available tools to annotate
def defineAnnoList() {
    return [
        'HaplotypeCaller',
        'Manta',
        'Mutect2',
        'Strelka',
        'TIDDIT'
    ]
}

// Define list of skipable QC tools
def defineSkipQClist() {
    return [
        'uniq',
        'bamqc',
        'mapq',
        'bcftools',
        'fastqc',
        'markduplicates',
        'multiqc',
        'samtools',
        'samtoolsstats',
        'vcftools',
        'versions'
    ]
}

// Define list of available step
def defineStepList() {
    return [
        'annotate',
        'mapping',
        'recalibrate',
        'variantcalling'
    ]
}

// Define list of available tools
def defineToolList() {
    return [
        'ascat',
        'freebayes',
        'haplotypecaller',
        'manta',
        'merge',
        'mpileup',
        'mutect2',
        'snpeff',
        'strelka',
        'tiddit',
        'vep'
    ]
}

// Channeling the TSV file containing BAM.
// Format is: "subject gender status sample bam bai"
def extractBam(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            checkNumberOfItem(row, 6)
            def idPatient = row[0]
            def gender    = row[1]
            def status    = returnStatus(row[2].toInteger())
            def idSample  = row[3]
            def bamFile   = returnFile(row[4])
            def baiFile   = returnFile(row[5])

            if (!hasExtension(bamFile, "bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
            if (!hasExtension(baiFile, "bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"

            return [idPatient, gender, status, idSample, bamFile, baiFile]
        }
}

// Create a channel of germline FASTQs from a directory pattern: "my_samples/*/"
// All FASTQ files in subdirectories are collected and emitted;
// they must have _R1_ and _R2_ in their names.
def extractFastqFromDir(pattern) {
    def fastq = Channel.create()
    // a temporary channel does all the work
    Channel
        .fromPath(pattern, type: 'dir')
        .ifEmpty { error "No directories found matching pattern '${pattern}'" }
        .subscribe onNext: { sampleDir ->
            // the last name of the sampleDir is assumed to be a unique sample id
            sampleId = sampleDir.getFileName().toString()

            for (path1 in file("${sampleDir}/**_R1_*.fastq.gz")) {
                assert path1.getName().contains('_R1_')
                path2 = file(path1.toString().replace('_R1_', '_R2_'))
                if (!path2.exists()) error "Path '${path2}' not found"
                (flowcell, lane) = flowcellLaneFromFastq(path1)
                patient = sampleId
                gender = 'ZZ'  // unused
                status = 0  // normal (not tumor)
                rgId = "${flowcell}.${sampleId}.${lane}"
                result = [patient, gender, status, sampleId, rgId, path1, path2]
                fastq.bind(result)
            }
    }, onComplete: { fastq.close() }
    fastq
}

// Extract gender and status from Channel
def extractInfos(channel) {
    def genderMap = [:]
    def statusMap = [:]
    channel = channel.map{ it ->
        def idPatient = it[0]
        def gender = it[1]
        def status = it[2]
        def idSample = it[3]
        genderMap[idPatient] = gender
        statusMap[idPatient, idSample] = status
        [idPatient] + it[3..-1]
    }
    [genderMap, statusMap, channel]
}

// Channeling the TSV file containing FASTQ or BAM
// Format is: "subject gender status sample lane fastq1 fastq2"
// or: "subject gender status sample lane bam"
def extractFastq(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            def idPatient  = row[0]
            def gender     = row[1]
            def status     = returnStatus(row[2].toInteger())
            def idSample   = row[3]
            def idRun      = row[4]
            def file1      = returnFile(row[5])
            def file2      = "null"
            if (hasExtension(file1, "fastq.gz") || hasExtension(file1, "fq.gz")) {
                checkNumberOfItem(row, 7)
                file2 = returnFile(row[6])
            if (!hasExtension(file2, "fastq.gz") && !hasExtension(file2, "fq.gz")) exit 1, "File: ${file2} has the wrong extension. See --help for more information"
        }
        else if (hasExtension(file1, "bam")) checkNumberOfItem(row, 6)
        else "No recognisable extention for input file: ${file1}"

        [idPatient, gender, status, idSample, idRun, file1, file2]
    }
}

// Channeling the TSV file containing Recalibration Tables.
// Format is: "subject gender status sample bam bai recalTables"
def extractRecal(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            checkNumberOfItem(row, 7)
            def idPatient  = row[0]
            def gender     = row[1]
            def status     = returnStatus(row[2].toInteger())
            def idSample   = row[3]
            def bamFile    = returnFile(row[4])
            def baiFile    = returnFile(row[5])
            def recalTable = returnFile(row[6])

            if (!hasExtension(bamFile, "bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
            if (!hasExtension(baiFile, "bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"
            if (!hasExtension(recalTable, "recal.table")) exit 1, "File: ${recalTable} has the wrong extension. See --help for more information"

            [idPatient, gender, status, idSample, bamFile, baiFile, recalTable]
    }
}

// Parse first line of a FASTQ file, return the flowcell id and lane number.
def flowcellLaneFromFastq(path) {
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    InputStream fileStream = new FileInputStream(path.toFile())
    InputStream gzipStream = new java.util.zip.GZIPInputStream(fileStream)
    Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
    BufferedReader buffered = new BufferedReader(decoder)
    def line = buffered.readLine()
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(' ')[0].split(':')
    String fcid
    int lane
    if (fields.size() == 7) {
        // CASAVA 1.8+ format
        fcid = fields[2]
        lane = fields[3].toInteger()
    } else if (fields.size() == 5) {
        fcid = fields[0]
        lane = fields[1].toInteger()
    }
    [fcid, lane]
}

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// Return file if it exists
def returnFile(it) {
    if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
    return file(it)
}

// Remove .ann .gz and .vcf extension from a VCF file
def reduceVCF(file) {
    return file.fileName.toString().minus(".ann").minus(".vcf").minus(".gz")
}

// Return status [0,1]
// 0 == Normal, 1 == Tumor
def returnStatus(it) {
    if (!(it in [0, 1])) exit 1, "Status is not recognized in TSV file: ${it}, see --help for more information"
    return it
}
