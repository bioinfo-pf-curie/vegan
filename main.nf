#!/usr/bin/env nextflow

/*
================================================================================
                                  nf-core/sarek
================================================================================
Started March 2016.
Ported to nf-core May 2019.
--------------------------------------------------------------------------------
nf-core/sarek:
  An open-source analysis pipeline to detect germline or somatic variants
  from whole genome or targeted sequencing
--------------------------------------------------------------------------------
 @Homepage
 https://sarek.scilifelab.se/
--------------------------------------------------------------------------------
 @Documentation
 https://github.com/nf-core/sarek/README.md
--------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/sarek --input sample.tsv -profile docker

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
        --sentieon                  If sentieon is available, will enable it for preprocessing, and variant calling
                                    Adds the following tools for --tools: DNAseq, DNAscope and TNscope
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
        --outdir                    The output directory where the results will be saved
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
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
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
    if (params.sentieon) {
        if (step == 'variantcalling') tsvPath =  "${params.outdir}/Preprocessing/TSV/recalibrated_sentieon.tsv"
        else exit 1, "Not possible to restart from that step"
    }
    else {
        tsvPath = step == 'recalibrate' ? "${params.outdir}/Preprocessing/TSV/duplicateMarked.tsv" : "${params.outdir}/Preprocessing/TSV/recalibrated.tsv"
    }
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
log.info nfcoreHeader()
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
summary['Output dir']        = params.outdir
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

    publishDir path:"${params.outdir}/pipeline_info", mode: params.publishDirMode

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

    publishDir params.outdir, mode: params.publishDirMode,
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

    publishDir params.outdir, mode: params.publishDirMode,
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

    publishDir params.outdir, mode: params.publishDirMode,
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

    publishDir params.outdir, mode: params.publishDirMode,
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

    publishDir params.outdir, mode: params.publishDirMode,
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

    publishDir params.outdir, mode: params.publishDirMode,
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

    publishDir params.outdir, mode: params.publishDirMode,
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

  publishDir params.outdir, mode: params.publishDirMode,
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

    publishDir "${params.outdir}/Reports/${idSample}/FastQC/${idSample}_${idRun}", mode: params.publishDirMode

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

    publishDir "${params.outdir}/Reports/${idSample}/FastQC/${idSample}_${idRun}", mode: params.publishDirMode

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

(inputPairReads, inputPairReadsSentieon) = inputPairReads.into(2)
if (params.sentieon) inputPairReads.close()
else inputPairReadsSentieon.close()

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
        bwa mem -K 100000000 -R \"${readGroup}\" ${extra} -t ${task.cpus} -M ${fasta} \
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

// STEP 1': MAPPING READS TO REFERENCE GENOME WITH SENTIEON BWA MEM

process SentieonMapReads {
    label 'sentieon'
    label 'cpus_max'
    label 'memory_max'

    tag {idPatient + "-" + idRun}

    input:
        set idPatient, idSample, idRun, file(inputFile1), file(inputFile2) from inputPairReadsSentieon
        file(bwaIndex) from ch_bwaIndex
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fastaFai

    output:
        set idPatient, idSample, idRun, file("${idSample}_${idRun}.bam") into bamMappedSentieon
        set idPatient, idSample, file("${idSample}_${idRun}.bam") into bamMappedSentieonBamQC

    when: params.sentieon

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
    """
    sentieon bwa mem -K 100000000 -R \"${readGroup}\" ${extra} -t ${task.cpus} -M ${fasta} \
    ${inputFile1} ${inputFile2} | \
    sentieon util sort -r ${fasta} -o ${idSample}_${idRun}.bam -t ${task.cpus} --sam2bam -i - 
    """
}

bamMappedSentieon = bamMappedSentieon.dump(tag:'Sentieon Mapped BAM')
// Sort BAM whether they are standalone or should be merged

singleBamSentieon = Channel.create()
multipleBamSentieon = Channel.create()
bamMappedSentieon.groupTuple(by:[0, 1])
    .choice(singleBamSentieon, multipleBamSentieon) {it[2].size() > 1 ? 1 : 0}
singleBamSentieon = singleBamSentieon.map {
    idPatient, idSample, idRun, bam ->
    [idPatient, idSample, bam]
}
singleBamSentieon = singleBamSentieon.dump(tag:'Single BAM')

// STEP 1.5: MERGING BAM FROM MULTIPLE LANES

multipleBam = multipleBam.mix(multipleBamSentieon)

process MergeBamMapped {
    label 'samtools' 
    label 'cpus_8'

    tag {idPatient + "-" + idSample}

    input:
        set idPatient, idSample, idRun, file(bam) from multipleBam

    output:
        set idPatient, idSample, file("${idSample}.bam") into mergedBam

    script:
    """
    samtools merge --threads ${task.cpus} ${idSample}.bam ${bam}
    """
}

mergedBam = mergedBam.dump(tag:'Merged BAM')

mergedBam = mergedBam.mix(singleBam,singleBamSentieon)

(mergedBam, mergedBamForSentieon) = mergedBam.into(2)

if (!params.sentieon) mergedBamForSentieon.close()
else mergedBam.close()

mergedBam = mergedBam.dump(tag:'BAMs for MD')
mergedBamForSentieon = mergedBamForSentieon.dump(tag:'Sentieon BAMs to Index')

process IndexBamMergedForSentieon { 
    label 'samtools' 
    label 'cpus_8'

    tag {idPatient + "-" + idSample}

    input:
        set idPatient, idSample, file(bam) from mergedBamForSentieon

    output:
        set idPatient, idSample, file(bam), file("${idSample}.bam.bai") into bamForSentieonDedup

    script:
    """
    samtools index ${bam}
    """
}

(mergedBam, mergedBamToIndex) = mergedBam.into(2)

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

// STEP 2: MARKING DUPLICATES

process MarkDuplicates { 
    label 'gatk' 
    label 'cpus_16'

    tag {idPatient + "-" + idSample}

    publishDir params.outdir, mode: params.publishDirMode,
        saveAs: {
            if (it == "${idSample}.bam.metrics" && 'markduplicates' in skipQC) null
            else if (it == "${idSample}.bam.metrics") "Reports/${idSample}/MarkDuplicates/${it}"
            else "Preprocessing/${idSample}/DuplicateMarked/${it}"
        }

    input:
        set idPatient, idSample, file("${idSample}.bam") from mergedBam

    output:
        set idPatient, idSample, file("${idSample}.md.bam"), file("${idSample}.md.bai") into duplicateMarkedBams
        file ("${idSample}.bam.metrics") into markDuplicatesReport

    when: params.knownIndels

    script:
    markdup_java_options = task.memory.toGiga() > 8 ? params.markdupJavaOptions : "\"-Xms" +  (task.memory.toGiga() / 2).trunc() + "g -Xmx" + (task.memory.toGiga() - 1) + "g\""
    """
    gatk --java-options ${markdup_java_options} \
        MarkDuplicates \
        --MAX_RECORDS_IN_RAM 50000 \
        --INPUT ${idSample}.bam \
        --METRICS_FILE ${idSample}.bam.metrics \
        --TMP_DIR . \
        --ASSUME_SORT_ORDER coordinate \
        --CREATE_INDEX true \
        --OUTPUT ${idSample}.md.bam
    """
}

if ('markduplicates' in skipQC) markDuplicatesReport.close()

duplicateMarkedBams = duplicateMarkedBams.dump(tag:'MD BAM')
markDuplicatesReport = markDuplicatesReport.dump(tag:'MD Report')

(bamMD, bamMDToJoin) = duplicateMarkedBams.into(2)

bamBaseRecalibrator = bamMD.combine(intBaseRecalibrator)

bamBaseRecalibrator = bamBaseRecalibrator.dump(tag:'BAM FOR BASERECALIBRATOR')

// STEP 2': SENTIEON DEDUP

process SentieonDedup { 
    label 'sentieon'
    label 'cpus_max'
    label 'memory_max'

    tag {idPatient + "-" + idSample}

    publishDir params.outdir, mode: params.publishDirMode,
        saveAs: {
            if (it == "${idSample}_*.txt" && 'sentieon' in skipQC) null
            else if (it == "${idSample}_*.txt") "Reports/${idSample}/Sentieon/${it}"
            else null
        }

    input:
        set idPatient, idSample, file(bam), file(bai) from bamForSentieonDedup
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fastaFai

    output:
        set idPatient, idSample, file("${idSample}.deduped.bam"), file("${idSample}.deduped.bam.bai") into bamDedupedSentieon
        file("${idSample}_*.txt") into bamDedupedSentieonQC

    when: params.sentieon

    script:
    """
    sentieon driver \
        -t ${task.cpus} \
        -i ${bam} \
        -r ${fasta} \
        --algo GCBias --summary ${idSample}_gc_summary.txt ${idSample}_gc_metric.txt \
        --algo MeanQualityByCycle ${idSample}_mq_metric.txt \
        --algo QualDistribution ${idSample}_qd_metric.txt \
        --algo InsertSizeMetricAlgo ${idSample}_is_metric.txt  \
        --algo AlignmentStat ${idSample}_aln_metric.txt

    sentieon driver \
        -t ${task.cpus} \
        -i ${bam} \
        --algo LocusCollector \
        --fun score_info ${idSample}_score.gz

    sentieon driver \
        -t ${task.cpus} \
        -i ${bam} \
        --algo Dedup \
        --rmdup \
        --score_info ${idSample}_score.gz  \
        --metrics ${idSample}_dedup_metric.txt ${idSample}.deduped.bam
    """
}




 
/*
================================================================================
                                     MultiQC
================================================================================
*/

// STEP MULTIQC

process MultiQC { 
    label 'multiqc' 
    publishDir "${params.outdir}/Reports/MultiQC", mode: params.publishDirMode

    input:
        file (multiqcConfig) from Channel.value(params.multiqc_config ? file(params.multiqc_config) : "")
        file (versions) from yamlSoftwareVersion
        file ('FastQC/*') from fastQCReport.collect().ifEmpty([])
        file ('MarkDuplicates/*') from markDuplicatesReport.collect().ifEmpty([])

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
    def subject = "[nf-core/sarek] Successful: ${workflow.runName}"
    if (!workflow.success) subject = "[nf-core/sarek] FAILED: ${workflow.runName}"
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
                log.warn "[nf-core/sarek] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/sarek] Could not attach MultiQC report to summary email"
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
            log.info "[nf-core/sarek] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            [ 'mail', '-s', subject, params.email ].execute() << email_txt
            log.info "[nf-core/sarek] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
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

    if (workflow.success) log.info "${c_purple}[nf-core/sarek]${c_green} Pipeline completed successfully${c_reset}"
    else {
        checkHostname()
        log.info "${c_purple}[nf-core/sarek]${c_red} Pipeline completed with errors${c_reset}"
    }
}

/*
================================================================================
                                nf-core functions
================================================================================
*/

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-sarek-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/sarek Workflow Summary'
    section_href: 'https://github.com/nf-core/sarek'
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

    ${c_purple}  nf-core/sarek v${workflow.manifest.version}${c_reset}
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

/*
================================================================================
                                 sarek functions
================================================================================
*/

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
        'bamqc',
        'bcftools',
        'fastqc',
        'markduplicates',
        'multiqc',
        'samtools',
        'sentieon',
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
        'controlfreec',
        'dnascope',
        'dnaseq',
        'freebayes',
        'haplotypecaller',
        'manta',
        'merge',
        'mpileup',
        'mutect2',
        'snpeff',
        'strelka',
        'tiddit',
        'tnscope',
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
