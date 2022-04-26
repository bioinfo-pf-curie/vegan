#!/usr/bin/env nextflow

/*
Copyright Institut Curie 2019-2022
This software is a computer program whose purpose is to analyze high-throughput sequencing data.
You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.
Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data.
The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.
*/

/*
========================================================================================
                         project : EUCANCAN/vegan
========================================================================================
VEGAN: Variant calling pipeline for whole Exome and whole Genome sequencing cANcer data Pipeline.
  An open-source analysis pipeline to detect germline or somatic variants
  from whole genome or targeted sequencing

@Homepage
https://gitlab.curie.fr/data-analysis/vegan
@Documentation
https://gitlab.curie.fr/data-analysis/vegan/README.md
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

// Initialize lintedParams and paramsWithUsage
NFTools.welcome(workflow, params)

/*
===================================
  SET UP CONFIGURATION VARIABLES
===================================
*/

// Use lintedParams as default params object
paramsWithUsage = NFTools.readParamsFromJsonSettings("${projectDir}/parameters.settings.json")
params.putAll(NFTools.lint(params, paramsWithUsage))

// Run name
customRunName = NFTools.checkRunName(workflow.runName, params.name)

// Custom functions/variables
mqcReport = []
include {checkAlignmentPercent} from './lib/functions'
include {loadDesign} from './lib/functions'


tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase()} : []


/*
================================================================================
                               CHECKING REFERENCES
================================================================================
*/

// ???
// CHECKING VARIABLES ??
// Init Channels ??
// Intervals ?

// Genome-based variables
if (!params.genome){
  exit 1, "No genome provided. The --genome option is mandatory"
}

if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
  exit 1, "The provided genome '${params.genome}' is not available in the genomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Initialize variable from the genome.conf file
params.bwaIndex = NFTools.getGenomeAttribute(params, 'bwaIndex')
params.chrLength = NFTools.getGenomeAttribute(params, 'chrLength')
params.dict = NFTools.getGenomeAttribute(params, 'dict')
params.fasta = NFTools.getGenomeAttribute(params, 'fasta')
params.fastaFai = NFTools.getGenomeAttribute(params, 'fastaFai')
params.gtf = NFTools.getGenomeAttribute(params, 'gtf')
params.dbsnp = NFTools.getGenomeAttribute(params, 'dbsnp')
params.dbsnpIndex = NFTools.getGenomeAttribute(params, 'dbsnpIndex')
params.acLoci = NFTools.getGenomeAttribute(params, 'acLoci')
params.acLociGC = NFTools.getGenomeAttribute(params, 'acLociGC')
params.polyms = NFTools.getGenomeAttribute(params, 'polyms')
params.germlineResource = NFTools.getGenomeAttribute(params, 'germlineResource')
params.germlineResourceIndex = NFTools.getGenomeAttribute(params, 'germlineResourceIndex')
params.intervals = NFTools.getGenomeAttribute(params, 'intervals')
params.knownIndels = NFTools.getGenomeAttribute(params, 'knownIndels')
params.knownIndelsIndex = NFTools.getGenomeAttribute(params, 'knownIndelsIndex')
params.snpeffDb = NFTools.getGenomeAttribute(params, 'snpeffDb')
params.snpeffCache = NFTools.getGenomeAttribute(params, 'snpeffCache')


// Stage config files
chMultiqcConfig = Channel.fromPath(params.multiqcConfig)
chOutputDocs = Channel.fromPath("$baseDir/docs/output.md")
chOutputDocsImages = file("$baseDir/docs/images/", checkIfExists: true)

/*
==========================
 VALIDATE INPUTS
==========================
*/

if ((params.reads && params.samplePlan) || (params.readPaths && params.samplePlan)){
  exit 1, "Input reads must be defined using either '--reads' or '--samplePlan' parameter. Please choose one way"
}

/*
==========================
 BUILD CHANNELS
==========================
*/

// Channels from genome.config
if ( params.bwaIndex ){
  Channel
    .fromPath("${params.bwaIndex}")
    .ifEmpty { exit 1, "Bwa index not found: ${params.bwaIndex}" }
    .set{chBwaIndex}
}else{
  exit 1, "No genome index specified!"
}

if (params.chrLength) {
  Channel
    .fromPath(params.chrLength, checkIfExists: true)
    .set{chrLength}
}else {
  chrLength = Channel.empty()
}

// Dict file
if ( params.dict ){
  Channel
    .fromPath(params.dict, checkIfExists: true)
    .set{chDict}
}
else{
  exit 1, "Fasta file not found: ${params.fasta}"
}

// Genome Fasta file
if ( params.fasta ){
  Channel
    .fromPath(params.fasta, checkIfExists: true)
    .set{chFasta}
}
else{
  exit 1, "Fasta file not found: ${params.fasta}"
}

if ( params.fastaFai ){
  Channel
    .fromPath(params.fastaFai, checkIfExists: true)
    .set{chFastaFai}
}
else{
  exit 1, "Fasta file not found: ${params.fasta}"
}

if (params.gtf) {
  Channel
    .fromPath(params.gtf, checkIfExists: true)
    .set{chGtf}
}else {
  exit 1, "GTF annotation file not specified!"
}

if (params.dbsnp) {
  Channel
    .fromPath(params.dbsnp, checkIfExists: true)
    .set{chDbsnp}
}else {
  exit 1, "Dbsnp annotation file not specified!"
}

if (params.dbsnpIndex) {
  Channel
    .fromPath(params.dbsnpIndex, checkIfExists: true)
    .set{chDbsnpIndex}
}else {
  exit 1, "Dbsnp annotation file not specified!"
}

if (params.acLoci) {
  Channel
    .fromPath(params.acLoci, checkIfExists: true)
    .set{chAcLoci}
}else {
  exit 1, "Ascat annotation file not specified!"
}

if (params.acLociGC) {
  Channel
    .fromPath(params.acLociGC, checkIfExists: true)
    .set{chAcLociGC}
}else {
  exit 1, "Ascat annotation file not specified!"
}

if (params.polyms) {
  Channel
    .fromPath(params.polyms, checkIfExists: true)
    .set{chPolyms}
}else {
  exit 1, "Identito vigilance annotation file not specified!"
}

if (params.germlineResource) {
  Channel
    .fromPath(params.germlineResource, checkIfExists: true)
    .set{chGermlineResource}
}else {
  exit 1, "Mutect2 annotation file not specified!"
}

if (params.germlineResourceIndex) {
  Channel
    .fromPath(params.germlineResourceIndex, checkIfExists: true)
    .set{chGermlineResourceIndex}
}else {
  exit 1, "Mutect2 annotation file not specified!"
}

if (params.intervals) {
  Channel
    .fromPath(params.intervals, checkIfExists: true)
    .set{chrIntervals}
}else {
  chrIntervals = Channel.empty()
}

if (params.knownIndels) {
  Channel
    .fromPath(params.knownIndels, checkIfExists: true)
    .set{chKnownIndels}
}else {
  exit 1, "BQSR annotation file not specified!"
}

if (params.knownIndelsIndex) {
  Channel
    .fromPath(params.knownIndelsIndex, checkIfExists: true)
    .set{chKnownIndelsIndex}
}else {
  exit 1, "BQSR annotation file not specified!"
}

// if (params.snpeffDb) {
//   Channel
//     .fromPath(params.snpeffDb, checkIfExists: true)
//     .set{chsnpeffDb}
// }else {
//   exit 1, "SnpEff annotation file not specified!"
// }

if (params.snpeffCache) {
  Channel
    .fromPath(params.snpeffCache, checkIfExists: true)
    .set{chSnpeffCache}
}else {
  exit 1, "SnpEff annotation file not specified!"
}

// Other channels

if ( params.metadata ){
  Channel
    .fromPath( params.metadata )
    .ifEmpty { exit 1, "Metadata file not found: ${params.metadata}" }
    .set { chMetadata }
}

if (params.targetBed) {
  Channel
    .fromPath(params.targetBed, checkIfExists: true)
    .set{chBed}
}else {
  chBed = Channel.empty()
}


/*
===========================
   SUMMARY
===========================
*/

summary = [
  'Pipeline Release': workflow.revision ?: null,
  'Run Name': customRunName,
  'Step' : params.step ?: null,
  'Profile' : workflow.profile,
  'Inputs' : params.samplePlan ?: params.reads ?: null,
  'Genome' : params.genome,
  'Tools' : params.tools ?: null,
  'Target Bed' : params.targetBed ?: null,
  'Identito' : params.polym ?: null,
  // 'SNV filters': params.SNVFilters ?: null,
  // 'SV filters': params.SVFilters ?: null,
  'QC tools skip' : params.skipQC ? 'Yes' : 'No',
  'Script dir': workflow.projectDir,
  'Launch Dir' : workflow.launchDir,
  'Output Dir' : params.outDir,
  'Working Dir': workflow.workDir,
  'Max Resources': "${params.maxMemory} memory, ${params.maxCpus} cpus, ${params.maxTime} time per job",
  'Container': workflow.containerEngine && workflow.container ? "${workflow.containerEngine} - ${workflow.container}" : null,
  'User': workflow.userName,
  'Config Description': params.configProfileDescription ?: null,
  'Config Contact': params.configProfileContact ?: null,
  'Config URL': params.configProfileUrl ?: null,
].findAll{ it.value != null }

workflowSummaryCh = NFTools.summarize(summary, workflow, params)

/*
==============================
  LOAD INPUT DATA
==============================
*/

// Load raw reads
chRawReads = NFTools.getInputData(params.samplePlan, params.reads, params.readPaths, params.singleEnd, params)
//chRawReads.view()

// Make samplePlan if not available
chSplan = NFTools.getSamplePlan(params.samplePlan, params.reads, params.readPaths, params.singleEnd)

// Load design file
if (params.design){
  Channel
    .fromPath(params.design)
    .ifEmpty { exit 1, "Design file not found: ${params.design}" }
    .set { chDesignFile }
  chDesign = loadDesign(params.design)
}else{
  chDesignFile = Channel.empty()
  chDesign = Channel.empty()
}

/*
==================================
           INCLUDE
==================================
*/

// Workflows
include { bwaMapping } from './nf-modules/local/subworkflow/bwaFlow'
include { bamFilters } from './nf-modules/local/subworkflow/bamFilteringFlow'
include { bamQcFlow } from './nf-modules/local/subworkflow/bamQcFlow'
include { identitoFlow } from './nf-modules/common/subworkflow/identito'
include { bqsrFlow } from './nf-modules/local/subworkflow/bqsrFlow'

// Processes
include { getSoftwareVersions } from './nf-modules/common/process/getSoftwareVersions'
include { outputDocumentation } from './nf-modules/common/process/outputDocumentation'
include { fastqc } from './nf-modules/common/process/fastqc'
include { multiqc } from './nf-modules/local/process/multiqc'
include { preseq } from './nf-modules/common/process/preseq'


/*
=====================================
            WORKFLOW
=====================================
*/

workflow {
  chVersions = Channel.empty()

  main:
    // Init MultiQC Channels
    //chRawReads.view()
    //chDesign.view()
    chFastqcMqc = Channel.empty()
    chPreseqMqc = Channel.empty()
    chIdentitoMqc = Channel.empty()

    // subroutines
    outputDocumentation(
      chOutputDocs,
      chOutputDocsImages
    )

    // PROCESS: fastqc
    if (! params.skipFastqc){
      fastqc(
        chRawReads
      )
      chFastqcMqc = fastqc.out.results.collect()
      chVersions = chVersions.mix(fastqc.out.versions)
    }

    //*******************************************
    // SUB-WORFKLOW : MAPPING BWA
    bwaMapping(
      chRawReads,
      chBwaIndex
    )

    chAlignedBam = bwaMapping.out.bam
    //chAlignedBam.view()
    chVersions = chVersions.mix(bwaMapping.out.versions)

    //*******************************************
    // Process : Preseq

    if (!params.skipPreseq){
    preseq(
      chAlignedBam
    )
    chPreseqMqc = preseq.out.results.collect()
    chVersions = chVersions.mix(preseq.out.versions)
  }

  //*******************************************
  // SUB-WORKFLOW : bamFiltering

    bamFilters(
      chAlignedBam,
      chBed
    )

  chVersions = chVersions.mix(bamFilters.out.versions)

  //*******************************************
  // SUB-WORKFLOW : GATK
  //gatk(
  //  bwaMapping.out.bam
  //)

  //*******************************************
  //SUB-WORKFLOW : bamQcFlow

  chFilteredBam = bamFilters.out.bam

  if (!params.skipQC){
    bamQcFlow(
      chFilteredBam,
      chBed,
      chGtf,
      chFasta,
      chDict
    )
    chVersions = chVersions.mix(bamQcFlow.out.versions)
  }

  // SUBWORKFLOW: Identito - polym and Monitoring
  if (!params.skipIdentito){
    identitoFlow(
      chFilteredBam,
      chFasta.collect(),
      chFastaFai.collect(),
      chPolyms.collect()
    )

    chIdentitoMqc = identitoFlow.out.results.collect()
    chVersions = chVersions.mix(identitoFlow.out.versions)
  }

  if('haplotypecaller' in tools || 'mutect2' in tools){
    bqsrFlow(
      chFilteredBam,
      chBed,
      chDbsnp,
      chDbsnpIndex,
      chFasta,
      chFastaFai,
      chKnownIndels,
      chKnownIndelsIndex,
      chDict
      )

    chVersions = chVersions.mix(bqsrFlow.out.versions)
  }
    //*******************************************
    // MULTIQC

    // Warnings that will be printed in the mqc report
    chWarn = Channel.empty()

    if (!params.skipMultiQC){

      getSoftwareVersions(
        chVersions.unique().collectFile()
      )

      multiqc(
        customRunName,
        chSplan.collect(),
        chMetadata.ifEmpty([]),
        chMultiqcConfig.ifEmpty([]),
        chFastqcMqc.ifEmpty([]),
        getSoftwareVersions.out.versionsYaml.collect().ifEmpty([]),
        workflowSummaryCh.collectFile(name: "workflow_summary_mqc.yaml"),
        chWarn.collect().ifEmpty([])
      )
      mqcReport = multiqc.out.report.toList()
    }
}

workflow.onComplete {
  NFTools.makeReports(workflow, params, summary, customRunName, mqcReport)
}
