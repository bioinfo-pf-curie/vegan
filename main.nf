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
========================================
 REFERENCES PARAMS FROM GENOME CONFIG
========================================
*/

// Genome-based variables
if (!params.genome){
  exit 1, "No genome provided. The --genome option is mandatory"
}

if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
  exit 1, "The provided genome '${params.genome}' is not available in the genomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Initialize variable from the genome.conf file
params.bwaIndex = NFTools.getGenomeAttribute(params, 'bwaIndex')
params.bwamem2Index = NFTools.getGenomeAttribute(params, 'bwamem2Index')
params.dragmapIndex = NFTools.getGenomeAttribute(params, 'dragmapIndex')
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
params.pileupSum = NFTools.getGenomeAttribute(params, 'pileupSum')
params.pileupSumIndex = NFTools.getGenomeAttribute(params, 'pileupSumIndex')
params.intervals = NFTools.getGenomeAttribute(params, 'intervals')
params.knownIndels = NFTools.getGenomeAttribute(params, 'knownIndels')
params.knownIndelsIndex = NFTools.getGenomeAttribute(params, 'knownIndelsIndex')
params.snpeffDb = NFTools.getGenomeAttribute(params, 'snpeffDb')
params.snpeffCache = NFTools.getGenomeAttribute(params, 'snpeffCache')
params.cosmicDb = NFTools.getGenomeAttribute(params, 'cosmicDb')
params.cosmicDbIndex = NFTools.getGenomeAttribute(params, 'cosmicDbIndex')
params.icgcDb = NFTools.getGenomeAttribute(params, 'icgcDb')
params.icgcDbIndex = NFTools.getGenomeAttribute(params, 'icgcDbIndex')
params.cancerhotspotsDb = NFTools.getGenomeAttribute(params, 'cancerhotspotsDb')
params.cancerhotspotsDbIndex = NFTools.getGenomeAttribute(params, 'cancerhotspotsDbIndex')
params.gnomadDb = NFTools.getGenomeAttribute(params, 'gnomadDb')
params.gnomadDbIndex = NFTools.getGenomeAttribute(params, 'gnomadDbIndex')
params.dbnsfp = NFTools.getGenomeAttribute(params, 'dbnsfp')
params.dbnsfpIndex = NFTools.getGenomeAttribute(params, 'dbnsfpIndex')
params.effGenomeSize = NFTools.getGenomeAttribute(params, 'effGenomeSize')

// Stage config files
chMultiqcConfig = Channel.fromPath(params.multiqcConfig)
chOutputDocs = Channel.fromPath("$baseDir/docs/output.md")
chOutputDocsImages = file("$baseDir/docs/images/", checkIfExists: true)

/*
==========================
 VALIDATE INPUTS
==========================
*/

if (params.noIntervals && params.targetBed){
  exit 1, "Options '--noIntervals' cannot be used together with '--targetBed' !"
}

if (params.pon && !params.ponIndex){
  exit 1, "Panel of normals (--pon) option requires an index file. Please use '--ponIndex' to specify it !"
}

if (!params.skipBQSR && ('haplotypecaller' in tools || 'mutect2' in tools) && (!params.dbsnp || !params.knownIndels || !params.dbsnpIndex || !params.knownIndelsIndex)){
  exit 1, "Missing annotation file(s) for GATK Base Recalibrator (dbSNP, knownIndels): Please use '--skipBQSR'"
}

if (!params.skipMutectContamination && 'mutect2' in tools && !params.germlineResource){
  exit 1, "Missing annotation file(s) for Mutect2 filtering (germlineResource): Please use '--skipMutectContamination'"
}

if (tools && ('ascat' in tools) && (!params.acLoci || !params.acLociGC)) {
  log.info """\
========================================================================
  ${colors.redBold}WARNING${colors.reset}: Missing annotation file(s) for ASCAT (acLoci, acLociGC).
========================================================================"""
  tools.removeAll(["ascat"])
}

if (tools && ('facets' in tools) && (!params.dbsnp)) {
  log.info """\
========================================================================
  ${colors.redBold}WARNING${colors.reset}: Missing annotation file(s) for Facets (dbsnp).
========================================================================"""
  tools.removeAll(["facets"])
}

/*
=====================================
 INIT FILE CHANNELS BASED ON PARAMS
=====================================
*/

chMetadata              = params.metadata              ? Channel.fromPath(params.metadata, checkIfExists: true).collect()               : Channel.empty()
chChrLength             = params.chrLength             ? Channel.fromPath(params.chrLength, checkIfExists: true).collect()              : Channel.empty()
chFasta                 = params.fasta                 ? Channel.fromPath(params.fasta, checkIfExists: true).collect()                  : Channel.empty()
chFastaFai              = params.fastaFai              ? Channel.fromPath(params.fastaFai, checkIfExists: true).collect()               : Channel.empty()
chDict                  = params.dict                  ? Channel.fromPath(params.dict, checkIfExists: true).collect()                   : Channel.empty()
chGtf                   = params.gtf                   ? Channel.fromPath(params.gtf, checkIfExists: true).collect()                    : Channel.value([]) //optional
chDbsnp                 = params.dbsnp                 ? Channel.fromPath(params.dbsnp, checkIfExists: true).collect()                  : Channel.value([]) //optional
chDbsnpIndex            = params.dbsnpIndex            ? Channel.fromPath(params.dbsnpIndex, checkIfExists: true).collect()             : Channel.value([]) //optional
chAcLoci                = params.acLoci                ? Channel.fromPath(params.acLoci, checkIfExists: true).collect()                 : Channel.value([]) //optional
chAcLociGC              = params.acLociGC              ? Channel.fromPath(params.acLociGC, checkIfExists: true).collect()               : Channel.value([]) //optional
chPolyms                = params.polyms                ? Channel.fromPath(params.polyms, checkIfExists: true).collect()                 : Channel.value([]) //optional
chGermlineResource      = params.germlineResource      ? Channel.fromPath(params.germlineResource, checkIfExists: true).collect()       : Channel.value([]) //optional
chGermlineResourceIndex = params.germlineResourceIndex ? Channel.fromPath(params.germlineResourceIndex, checkIfExists: true).collect()  : Channel.value([]) //optional
chPileupSum             = params.pileupSum             ? Channel.fromPath(params.pileupSum, checkIfExists: true).collect()              : Channel.value([]) //optional
chPileupSumIndex        = params.pileupSumIndex        ? Channel.fromPath(params.pileupSumIndex, checkIfExists: true).collect()         : Channel.value([]) //optional
chPon                   = params.pon                   ? Channel.fromPath(params.pon, checkIfExists: true).collect()                    : Channel.value([]) //optional
chPonIndex              = params.ponIndex              ? Channel.fromPath(params.ponIndex, checkIfExists: true).collect()               : Channel.value([]) //optional
chKnownIndels           = params.knownIndels           ? Channel.fromPath(params.knownIndels, checkIfExists: true).collect()            : Channel.value([]) //optional
chKnownIndelsIndex      = params.knownIndelsIndex      ? Channel.fromPath(params.knownIndelsIndex, checkIfExists: true).collect()       : Channel.value([]) //optional
chSnpeffDb              = params.snpeffDb              ? Channel.value(params.snpeffDb)                                                 : Channel.empty()
chSnpeffCache           = params.snpeffCache           ? Channel.fromPath(params.snpeffCache, checkIfExists: true).collect()            : Channel.value([]) //optional
chCosmicDb              = params.cosmicDb              ? Channel.fromPath(params.cosmicDb, checkIfExists: true).collect()               : Channel.empty()
chCosmicDbIndex         = params.cosmicDbIndex         ? Channel.fromPath(params.cosmicDbIndex, checkIfExists: true).collect()          : Channel.empty()
chIcgcDb                = params.icgcDb                ? Channel.fromPath(params.icgcDb, checkIfExists: true).collect()                 : Channel.empty()
chIcgcDbIndex           = params.icgcDbIndex           ? Channel.fromPath(params.icgcDbIndex, checkIfExists: true).collect()            : Channel.empty()
chCancerhotspotsDb      = params.cancerhotspotsDb      ? Channel.fromPath(params.cancerhotspotsDb, checkIfExists: true).collect()       : Channel.empty()
chCancerhotspotsDbIndex = params.cancerhotspotsDbIndex ? Channel.fromPath(params.cancerhotspotsDbIndex, checkIfExists: true).collect()  : Channel.empty()
chGnomadDb              = params.gnomadDb              ? Channel.fromPath(params.gnomadDb, checkIfExists: true).collect()               : Channel.empty()
chGnomadDbIndex         = params.gnomadDbIndex         ? Channel.fromPath(params.gnomadDbIndex, checkIfExists: true).collect()          : Channel.empty()
chDbnsfp                = params.dbnsfp                ? Channel.fromPath(params.dbnsfp, checkIfExists: true).collect()                 : Channel.empty()
chDbnsfpIndex           = params.dbnsfpIndex           ? Channel.fromPath(params.dbnsfpIndex, checkIfExists: true).collect()            : Channel.empty()
chTargetBed             = params.targetBed             ? Channel.fromPath(params.targetBed, checkIfExists: true).collect()              : Channel.value([]) //optional
chIntervals             = params.intervals             ? Channel.fromPath(params.intervals, checkIfExists: true).collect()              : Channel.value([]) //optional
chBwaIndex              = params.bwaIndex              ? Channel.fromPath(params.bwaIndex)                                              : Channel.empty()
chBwaMem2Index          = params.bwamem2Index          ? Channel.fromPath(params.bwamem2Index)                                          : Channel.empty()
chDragmapIndex          = params.dragmapIndex          ? Channel.fromPath(params.dragmapIndex)                                          : Channel.empty()
chEffGenomeSize         = params.effGenomeSize         ? Channel.value(params.effGenomeSize)                                            : Channel.value([]) //optionals

/*
===========================
   SUMMARY
===========================
*/

summary = [
  'Pipeline' : workflow.manifest.name ?: null,
  'Version': workflow.manifest.version ?: null,
  'DOI': workflow.manifest.doi ?: null,
  'Run Name': customRunName,
  'Step' : params.step ?: null,
  'Profile' : workflow.profile,
  'Inputs' : params.samplePlan ?: params.reads ?: null,
  'Chunks' : params.splitFastq ? params.fastqChunksSize : null,
  'Design' : params.design ?: null,
  'Genome' : params.genome,
  'Intervals' : params.noIntervals || params.targetBed ? 'no' : 'yes',
  'Target Bed' : params.targetBed ?: null,
  'Pon' : params.pon ?: null,
  'Aligner':params.aligner ?: null,
  'Tools' : params.tools ?: null,
  'Databases' : params.annotDb ?: null,
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
if (params.step == "mapping"){
  chRawReads = NFTools.getInputData(params.samplePlan, params.reads, params.readPaths, params.singleEnd, params)
}else if (params.step == "filtering"){
  chRawReads = Channel.empty()
  chAlignedBam = NFTools.getIntermediatesData(params.samplePlan, ['.bam','.bai'],  params).map{it ->[it[0], it[1][0], it[1][1]]}
}else if (params.step == "calling"){
  chRawReads = Channel.empty()
  chAlignedBam = Channel.empty()
  chFilteredBam = NFTools.getIntermediatesData(params.samplePlan, ['.bam','.bai'],  params).map{it ->[it[0], it[1][0], it[1][1]]}
}else if (params.step == "annotate"){
  chRawReads = Channel.empty()
  chAlignedBam = Channel.empty()
  chFilteredBam = Channel.empty()
  chAllVcf = NFTools.getIntermediatesData(params.samplePlan, ['.vcf.gz','.tbi'],  params).map{it ->[it[0], it[1][0], it[1][1]]}
}

// Make samplePlan if not available
chSplan = NFTools.getSamplePlan(params.samplePlan, params.reads, params.readPaths, params.singleEnd)

// Load design file
if (params.design){
  Channel
    .fromPath(params.design)
    .ifEmpty { exit 1, "Design file not found: ${params.design}" }
    .set { chDesignFile }

  chDesign = loadDesign(params.design)

  //Separate the design in germline only / tumor only / germline + tumor
  chDesign.branch{
    paired: it[0] != '' && it[1] != ''
    tumorOnly: it[0] != '' && it[1] == ''
    germlineOnly: it[0] == '' && it[1] != ''
  }.set{ chDesign }

//  if (chDesign.tumorOnly.count() > 0 && !params.pon){
//    exit 1, "Tumor only samples detected without panels of normal: Please use '--pon'"
//  }

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
include { prepareIntervalsFlow } from './nf-modules/local/subworkflow/prepareIntervals'
include { mappingFlow } from './nf-modules/local/subworkflow/mapping'
include { bamFiltersFlow } from './nf-modules/local/subworkflow/bamFiltering'
include { bamQcFlow } from './nf-modules/local/subworkflow/bamQc'
include { identitoFlow } from './nf-modules/common/subworkflow/identito'
include { bqsrFlow } from './nf-modules/local/subworkflow/bqsr'
include { haplotypeCallerFlow } from './nf-modules/local/subworkflow/haplotypeCaller'
include { mutect2PairsFlow } from './nf-modules/local/subworkflow/mutect2Pairs'
include { mutect2TumorOnlyFlow } from './nf-modules/local/subworkflow/mutect2TumorOnly'
include { vcfQcFlow } from './nf-modules/local/subworkflow/vcfQc'
include { annotateSomaticFlow as  annotateFlow} from './nf-modules/local/subworkflow/annotateSomatic'
include { tableReportFlow } from './nf-modules/local/subworkflow/tableReport'
include { mantaFlow } from './nf-modules/local/subworkflow/manta'
include { tmbFlow } from './nf-modules/local/subworkflow/tmb'
include { msiFlow } from './nf-modules/local/subworkflow/msi'
include { facetsFlow } from './nf-modules/local/subworkflow/facets'
include { ascatFlow } from './nf-modules/local/subworkflow/ascat'

// Processes
include { getSoftwareVersions } from './nf-modules/common/process/utils/getSoftwareVersions'
include { outputDocumentation } from './nf-modules/common/process/utils/outputDocumentation'
include { fastqc } from './nf-modules/common/process/fastqc/fastqc'
include { multiqc } from './nf-modules/local/process/multiqc'
include { preseq } from './nf-modules/common/process/preseq/preseq'

/*
=====================================
 WORKFLOW
=====================================
*/

workflow {
  chVersions = Channel.empty()

  main:
    // Init MultiQC Channels
    chFastqcMqc = Channel.empty()
    chMappingMqc = Channel.empty()
    chMappingStats = Channel.empty()
    chMappingFlagstat = Channel.empty()
    chOntargetStatsMqc = Channel.empty()
    chFilteringStatsMqc = Channel.empty()
    chPreseqMqc = Channel.empty()
    chIdentitoMqc = Channel.empty()
    chGeneCovMqc = Channel.empty()
    chMosdepthMqc = Channel.empty()
    chFragSizeMqc = Channel.empty()
    chWgsMetricsMqc = Channel.empty()
    chHaplotypecallerMetricsMqc = Channel.empty()
    chMutect2MetricsMqc = Channel.empty()
    chTsTvMqc = Channel.empty()
    chSnpEffMqc = Channel.empty()
    chMSIMqc = Channel.empty()
    chTMBMqc = Channel.empty()

    // subroutines
    outputDocumentation(
      chOutputDocs,
      chOutputDocsImages
    )

    //********************************************
    // SUB-WORKFLOW : Prepare intervals list

    if (!params.noIntervals && !params.targetBed){
      prepareIntervalsFlow(
        chFastaFai.collect(),
        chIntervals
      )
      chVersions = chVersions.mix(prepareIntervalsFlow.out.versions)
      chIntervalBeds = prepareIntervalsFlow.out.intervals
    }else if (params.targetBed){
      chIntervalBeds = chTargetBed
    }else{
      chIntervalBeds = Channel.value([])
    }

    // PROCESS: fastqc
    fastqc(
      chRawReads
    )
    chFastqcMqc = fastqc.out.results.collect()
    chVersions = chVersions.mix(fastqc.out.versions)

    //*******************************************
    // SUB-WORFKLOW : MAPPING WITH BWA-MEM/BWA-MEM2/DRAGMAP

    if (params.step == "mapping"){
      chAlignerIndex = params.aligner == 'bwa-mem' ? chBwaIndex : params.aligner == 'bwa-mem2' ? chBwaMem2Index : chDragmapIndex
    
      mappingFlow(
        chRawReads,
        chAlignerIndex
      )
      chAlignedBam = mappingFlow.out.bam
      chMappingMqc = mappingFlow.out.logs
      chMappingStats = mappingFlow.out.stats
      chMappingFlagstat = mappingFlow.out.flagstat
      chVersions = chVersions.mix(mappingFlow.out.versions)
    }

    //*******************************************
    // Process : Preseq

    if (!params.skipPreseq){
      preseq(
        chAlignedBam
      )
      chPreseqMqc = preseq.out.results.collect()
      chVersions = chVersions.mix(preseq.out.versions)
    }

    if (params.step == "mapping" || params.step == "filtering"){

      //*******************************************
      // SUB-WORKFLOW : bamFiltering

      bamFiltersFlow(
        chAlignedBam,
        chTargetBed
      )
      chVersions = chVersions.mix(bamFiltersFlow.out.versions)
      chMarkdupStatsMqc = bamFiltersFlow.out.markdupFlagstats
      chOnTargetStatsMqc = bamFiltersFlow.out.onTargetFlagstats
      chFilteringStatsMqc = bamFiltersFlow.out.filteringFlagstats
      chFilteredBam = bamFiltersFlow.out.bam

      //*******************************************
      //SUB-WORKFLOW : bamQcFlow

      if (!params.skipQC){
        bamQcFlow(
          chFilteredBam,
          chTargetBed,
          chGtf,
          chFasta,
          chDict
        )
        chGeneCovMqc = bamQcFlow.out.geneCovMqc
        chMosdepthMqc = bamQcFlow.out.depth
        chFragSizeMqc = bamQcFlow.out.fragSize
        chWgsMetricsMqc = bamQcFlow.out.wgsMetrics
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

      //*****************************
      // GATK4 - PRE-PROCESSING

      bqsrFlow(
        chFilteredBam,
        chIntervalBeds,
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

    if (params.step != "annotate"){
      chProcBam = (params.skipBQSR || params.step == "calling") ? chFilteredBam : bqsrFlow.out.bam
    }

  /*
  ================================================================================
   DESIGN / PAIRED ANALYSIS / TUMOR ONLY / GERMLINE ONLY
  ================================================================================
  */

  chSingleBam = Channel.empty()

  if (params.design && (params.step == "mapping" || params.step == "filtering" || params.step == "calling")){

    //******************
    // PAIRED BAMS

    //[meta], tumor_bam, tumor_bai, normal_bam, normal_bai
    chProcBam
      .combine(chProcBam)
      .combine(chDesign.paired)
      .filter { it[0].id == it[6] && it[3].id == it[7] }
      .map{ it ->
        def meta = [tumor_id:it[6], normal_id:it[7], pair_id:it[8], id:it[6]+"_vs_"+it[7], status:"pair", sex:it[9]]
        return [meta, it[1], it[2], it[4], it[5] ]
      }.set{ chPairBam }

    //[meta], tumor_bam, tumor_bai
    chPairBam
      .map{ it ->
        def meta = [tumor_id:it[0].tumor_id, id:it[0].tumor_id, status: "tumor", sex:it[0].sex]
        return [meta, it[1], it[2] ]
      }.set{ chTumorBam }

    //[meta], normal_bam, normal_bai
    chPairBam
      .map{ it ->
        def meta = [normal_id:it[0].normal_id, id:it[0].normal_id, status: "normal", sex:it[0].sex]
        return [meta, it[3], it[4] ]
      }.set{ chNormalBam }

    chSingleBam = chNormalBam.mix(chTumorBam)

    //******************  
    // TUMOR ONLY 

    chProcBam
      .combine(chDesign.tumorOnly)
      .filter { it[0].id == it[3] }
      .map{ it ->
        def meta = [id:it[3], status: "tumor", sex:it[6]]
        return [meta, it[1], it[2] ]
      }.set{ chTumorOnlyBam }

    chSingleBam = chSingleBam.mix(chTumorOnlyBam)

    //*******************
    // GERMLINE ONLY

    chProcBam
      .combine(chDesign.germlineOnly)
      .filter { it[0].id == it[4] }
      .map{ it ->
        def meta = [id:it[4], status: "normal", sex:it[6]]
        return [meta, it[1], it[2] ]
      }.set{ chGermlineOnlyBam }

    chSingleBam = chSingleBam.mix(chGermlineOnlyBam)
   }

  /*
  ================================================================================
   SNV VARIANT CALLING
  ================================================================================
  */

  chAllVcf = Channel.empty()

  if (params.step == "mapping" || params.step == "filtering" || params.step == "calling"){

    //*******************************************
    //SUB-WORKFLOW : HaplotypeCaller

    chSingleBam = chSingleBam.unique()

    if('haplotypecaller' in tools){
      haplotypeCallerFlow(
        chSingleBam,
        chIntervalBeds,
        chDbsnp,
        chDbsnpIndex,
        chFasta,
        chFastaFai,
        chDict
      )
      chVersions = chVersions.mix(haplotypeCallerFlow.out.versions)
      //chHaplotypecallerMetricsMqc = haplotypeCallerFlow.out.mqc
      //chTsTvMqc = haplotypeCallerFlow.out.transition
      chAllVcf = chAllVcf.mix(haplotypeCallerFlow.out.vcfNorm)
    }

    //*******************************************
    //SUB-WORKFLOW : Mutect2

    chMutect2MetricsMqc = Channel.empty()
    
    if('mutect2' in tools){

      mutect2PairsFlow(
        chPairBam,
        chIntervalBeds,
        chFasta,
        chFastaFai,
        chDict,
        chGermlineResource,
        chGermlineResourceIndex,
        chPileupSum,
        chPileupSumIndex,
        chPon,
        chPonIndex,
      )
      chVersions = chVersions.mix(mutect2PairsFlow.out.versions)
      //chMutect2MetricsMqc = chMutect2MetricsMqc.mix(mutect2PairsFlow.out.mqc)
      //chTsTvMqc = chTsTvMqc.mix(mutect2PairsFlow.out.transition)
      chAllVcf = chAllVcf.mix(mutect2PairsFlow.out.vcfFiltered)

      mutect2TumorOnlyFlow(
        chTumorOnlyBam,
        chIntervalBeds,
        chFasta,
        chFastaFai,
        chDict,
        chGermlineResource,
        chGermlineResourceIndex,
        chPileupSum,
        chPileupSumIndex,
        chPon,
        chPonIndex
      )
      chVersions = chVersions.mix(mutect2TumorOnlyFlow.out.versions)
      //chMutect2MetricsMqc = chMutect2MetricsMqc.mix(mutect2TumorOnlyFlow.out.mqc)
      //chTsTvMqc = chTsTvMqc.mix(mutect2TumorOnlyFlow.out.transition)
      chAllVcf = chAllVcf.mix(mutect2TumorOnlyFlow.out.vcfFiltered)
    }
  }

  /*
  ================================================================================
                                   VCF QC
  ================================================================================
  */

  vcfQcFlow(
    chAllVcf
  )
  chVcfMetricsMqc = vcfQcFlow.out.mqc
  chTsTvMqc = vcfQcFlow.out.transition

  /*
  ================================================================================
                                   VCF ANNOTATION
  ================================================================================
  */

  // Annotation somatic vcf
  if('snpeff' in tools || params.step == 'annotate'){
    annotateFlow(
      chAllVcf,
      chSnpeffDb,
      chSnpeffCache,
      chCosmicDb,
      chCosmicDbIndex,
      chIcgcDb,
      chIcgcDbIndex,
      chCancerhotspotsDb,
      chCancerhotspotsDbIndex,
      chGnomadDb,
      chGnomadDbIndex,
      chDbnsfp,
      chDbnsfpIndex
    )
    chSnpEffMqc = annotateFlow.out.snpEffReport
    chTMB = annotateFlow.out.vcf.filter{ it[0].status != "normal" }
  }else{
    chTMB = Channel.empty()
  }

  /*
  ================================================================================
                                   TABLE REPORT
  ================================================================================
  */

  if('snpeff' in tools || params.step == 'annotate'){
    tableReportFlow(
      annotateFlow.out.vcf
    )
  }

  /*
  ================================================================================
                                         TMB
  ================================================================================
  */

  if('tmb' in tools || params.step == 'annotate'){
    tmbFlow(
      chTMB,
      chTargetBed,
      chEffGenomeSize
    )
    chTMBMqc = tmbFlow.out.report.map{it->it[1]}
  }

  /*
  ================================================================================
                                        MSI
  ================================================================================
  */

  if('msisensor' in tools && params.step != 'annotate'){
    msiFlow(
      chPairBam,
      chFasta,
      chTargetBed
    )
    chMSIMqc = msiFlow.out.report.map{it->it[1]}
  }

  /*
  ================================================================================
                             SV VARIANT CALLING
  ================================================================================
  */

  // STEP MANTA

  if ('manta' in tools && params.step != 'annotate'){
    mantaFlow(
      chPairBam,
      chTargetBed,
      chFasta,
      chFastaFai
    )
    chVersions = chVersions.mix(mantaFlow.out.versions)
  }


  /*
  ================================================================================
                                        CNV calling
  ================================================================================
  */

  if('facets' in tools && params.step != 'annotate'){
    facetsFlow(
      chPairBam,
      chDbsnp
    )
  }

  if('ascat' in tools && params.step != 'annotate'){
    ascatFlow(
      chSingleBam,
      chAcLoci,
      chAcLociGC,
      chDict,
      chFasta,
      chFastaFai
    )
  }

  /*
  ================================================================================
                                 MULTIQC
  ================================================================================
  */

  // Warnings that will be printed in the mqc report
  chWarn = Channel.empty()

  if (!params.skipMultiQC && params.step != 'calling' && params.step != 'annotate' ){
    getSoftwareVersions(
      chVersions.unique().collectFile()
    )

    multiqc(
      customRunName,
      chSplan.collect(),
      chMetadata.ifEmpty([]),
      chMultiqcConfig.ifEmpty([]),
      chFastqcMqc.collect().ifEmpty([]),
      chMappingMqc.collect().ifEmpty([]),
      chMappingStats.collect().ifEmpty([]),
      chMappingFlagstat.collect().ifEmpty([]),
      chFragSizeMqc.collect().ifEmpty([]),
      chWgsMetricsMqc.collect().ifEmpty([]),
      chPreseqMqc.collect().ifEmpty([]),
      chMarkdupStatsMqc.collect().ifEmpty([]),
      chOnTargetStatsMqc.collect().ifEmpty([]),
      chFilteringStatsMqc.collect().ifEmpty([]),
      chGeneCovMqc.collect().ifEmpty([]),
      chMosdepthMqc.collect().ifEmpty([]),
      chIdentitoMqc.collect().ifEmpty([]),
      chVcfMetricsMqc.collect().ifEmpty([]),
      chTsTvMqc.collect().ifEmpty([]),
      chSnpEffMqc.collect().ifEmpty([]),
      chMSIMqc.collect().ifEmpty([]),
      chTMBMqc.collect().ifEmpty([]),
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
