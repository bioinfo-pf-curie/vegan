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
params.bwaMem2Index = NFTools.getGenomeAttribute(params, 'bwaMem2Index')
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

/* TODO - add additional controls on annotation files ? */

/*
=====================================
 INIT FILE CHANNELS BASED ON PARAMS
=====================================
*/

chMetatdata             = params.metadata              ? Channel.fromPath(params.metadata, checkIfExists: true).collect()               : Channel.empty()
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
chPon                   = params.pon                   ? Channel.fromPath(params.pon, checkIfExists: true).collect()                    : Channel.value([]) //optional
chPonIndex              = params.ponIndex              ? Channel.fromPath(params.ponIndex, checkIfExists: true).collect()               : Channel.value([]) //optional
chKnownIndels           = params.knownIndels           ? Channel.fromPath(params.knownIndels, checkIfExists: true).collect()            : Channel.value([]) //optional
chKnownIndelsIndex      = params.knownIndelsIndex      ? Channel.fromPath(params.knownIndelsIndex, checkIfExists: true).collect()       : Channel.value([]) //optional
chSnpeffDb              = params.snpeffDb              ? Channel.of(params.snpeffDb)                                                    : Channel.empty()
chSnpeffCache           = params.snpeffCache           ? Channel.fromPath(params.snpeffCache, checkIfExists: true).collect()            : Channel.value([])

chBed                   = params.targetBed             ? Channel.fromPath(params.targetBed, checkIfExists: true).collect()              : Channel.value([]) //optional
chIntervals             = params.intervals             ? Channel.fromPath(params.intervals, checkIfExists: true).collect()              : Channel.value([]) //optional

chBwaIndex              = params.bwaIndex              ? Channel.fromPath(params.bwaIndex)                                              : Channel.empty()
chBwaMem2Index          = params.bwaMem2Index          ? Channel.fromPath(params.bwaMem2Index)                                          : Channel.empty()
chDragmapIndex          = params.dragmapIndex          ? Channel.fromPath(params.dragmapIndex)                                          : Channel.empty()


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
  'Design' : params.design ?: null,
  'Genome' : params.genome,
  'Tools' : params.tools ?: null,
  'Target Bed' : params.targetBed ?: null,
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
include { mappingFlow } from './nf-modules/local/subworkflow/mapping'
include { bamFiltersFlow } from './nf-modules/local/subworkflow/bamFiltering'
include { bamQcFlow } from './nf-modules/local/subworkflow/bamQc'
include { identitoFlow } from './nf-modules/common/subworkflow/identito'
include { bqsrFlow } from './nf-modules/local/subworkflow/bqsr'
include { haplotypeCallerFlow } from './nf-modules/local/subworkflow/haplotypeCaller'
include { mutect2PairsFlow } from './nf-modules/local/subworkflow/mutect2Pairs'
include { annotateFlow } from './nf-modules/local/subworkflow/annotate'
include { mantaFlow } from './nf-modules/local/subworkflow/manta'

// Processes
include { getSoftwareVersions } from './nf-modules/common/process/utils/getSoftwareVersions'
include { outputDocumentation } from './nf-modules/common/process/utils/outputDocumentation'
include { fastqc } from './nf-modules/common/process/fastqc/fastqc'
include { multiqc } from './nf-modules/local/process/multiqc'
include { preseq } from './nf-modules/common/process/preseq/preseq'
include { tmb } from './nf-modules/common/process/tmb/tmb'
include { msisensorproScan } from './nf-modules/common/process/msisensorpro/msisensorproScan.nf'
include { msisensorproMsi } from './nf-modules/common/process/msisensorpro/msisensorproMsi.nf'

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
    // SUB-WORFKLOW : MAPPING WITH BWA-MEM/BWA-MEM2/DRAGMAP

    chAlignerIndex = params.aligner == 'bwa-mem' ? chBwaIndex :
      params.aligner == 'bwa-mem2' ? chBwaMem2Index :
      chdragmapIndex

    mappingFlow(
      chRawReads,
      chAlignerIndex
    )
    chAlignedBam = mappingFlow.out.bam
    chVersions = chVersions.mix(mappingFlow.out.versions)

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

    bamFiltersFlow(
      chAlignedBam,
      chBed
    )
    chVersions = chVersions.mix(bamFiltersFlow.out.versions)
    chFilteredBam = bamFiltersFlow.out.bam

    //*******************************************
    //SUB-WORKFLOW : bamQcFlow

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

    //*****************************
    // GATK4 - PRE-PROCESSING

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
      chProcBam = bqsrFlow.out.bqsrBam
    }else{
      chProcBam = chFilteredBam
    }

  /*
  ================================================================================
   DESIGN / PAIRED ANALYSIS
  ================================================================================
  */

  //[meta], tumor_bam, tumor_bai, normal_bam, normal_bai
  chProcBam
    .combine(chProcBam)
    .combine(chDesign.paired)
    .filter { it[0].id == it[6] && it[3].id == it[7] }
    .map{ it ->
      meta = [tumor_id:it[6], normal_id:it[7], status: "pair", id:it[8], sex:it[9]]
      return [meta, it[1], it[2], it[4], it[5] ]
    }.set{ chPairBam }

  //[meta], tumor_bam, tumor_bai
  chProcBam
    .combine(chProcBam)
    .combine(chDesign.paired)
    .filter { it[0].id == it[6] && it[3].id == it[7] }
    .map{ it ->
      meta = [tumor_id:it[6], normal_id:it[7], status: "tumor", id:it[8], sex:it[9]]
      return [meta, it[1], it[2] ]
    }.set{ chTumorBam }

  //chTumorBam.view()

  //[meta], normal_bam, normal_bai
  chProcBam
    .combine(chProcBam)
    .combine(chDesign.paired)
    .filter { it[0].id == it[6] && it[3].id == it[7] }
    .map{ it ->
      meta = [tumor_id:it[6], normal_id:it[7], status: "normal", id:it[8], sex:it[9]]
      return [meta, it[4], it[5] ]
    }.set{ chNormalBam }

  chSingleBam = chNormalBam.mix(chTumorBam)

  /*
  ================================================================================
   SNV VARIANT CALLING
  ================================================================================
  */

  chAllVcf = Channel.empty()

  //*******************************************
  //SUB-WORKFLOW : HaplotypeCaller


  if(params.tools && params.tools.contains('haplotypecaller')){
    haplotypeCallerFlow(
      chProcBam,
      chBed,
      chDbsnp,
      chDbsnpIndex,
      chFasta,
      chFastaFai,
      chDict
    )
    chVersions = chVersions.mix(haplotypeCallerFlow.out.versions)
    chAllVcf = chAllVcf.mix(haplotypeCallerFlow.out.vcf)
  }

  //*******************************************
  //SUB-WORKFLOW : Mutect2

  if(params.tools && params.tools.contains('mutect2')){
    mutect2PairsFlow(
      chPairBam,
      chBed,
      chFasta,
      chFastaFai,
      chDict,
      chGermlineResource,
      chGermlineResourceIndex,
      chPon,
      chPonIndex,
      chIntervals
    )
    chVersions = chVersions.mix(mutect2PairsFlow.out.versions)
    chAllVcf = chAllVcf.mix(mutect2PairsFlow.out.vcfFiltered)
  }


  /*
  ================================================================================
                                   VCF ANNOTATION
  ================================================================================
  */

  // Annotation somatic vcf
  annotateFlow(
    mutect2PairsFlow.out.vcfFiltered,
    chSnpeffDb,
    chSnpeffCache
  )

  /*
  ================================================================================
                                         TMB
  ================================================================================
  */

  dbConfig = file("$projectDir/assets/tmb/snpeff.yml", checkIfExists: true)
  varConfig = file("$projectDir/assets/tmb/mutect2.yml", checkIfExists: true)

  //tmb(
  //  annotateFlow.out.vcf.map{ it -> [it[0], it[1][0], dbConfig, varConfig] },
  //  chBed
  //)
  //chVersions = chVersions.mix(tmb.out.versions)

  /*
  ================================================================================
                                        MSI
  ================================================================================
  */

  msisensorproScan(
    chFasta
  )

  msisensorproMsi(
    chPairBam,
    chFasta.collect(),
    msisensorproScan.out.list.collect()
  )
  chVersions = chVersions.mix(msisensorproMsi.out.versions)

  /*
  ================================================================================
                             SV VARIANT CALLING
  ================================================================================
  */

    // STEP MANTA.1 - SINGLE MODE

    if ('manta' in params.tools){
    mantaFlow(
      chPairBam,
      chSingleBam,
      chBed,
      chFasta,
      chFastaFai
      )

    chVersions = chVersions.mix(mantaFlow.out.versions)
  }


  /*
  ================================================================================
                                 MULTIQC
  ================================================================================
  */

    // MULTIQC

    // Warnings that will be printed in the mqc report
    // chWarn = Channel.empty()
    //
    // if (!params.skipMultiQC){
    //
    //   getSoftwareVersions(
    //     chVersions.unique().collectFile()
    //   )
    //
    //   multiqc(
    //     customRunName,
    //     chSplan.collect(),
    //     chMetadata.ifEmpty([]),
    //     chMultiqcConfig.ifEmpty([]),
    //     chFastqcMqc.ifEmpty([]),
    //     getSoftwareVersions.out.versionsYaml.collect().ifEmpty([]),
    //     workflowSummaryCh.collectFile(name: "workflow_summary_mqc.yaml"),
    //     chWarn.collect().ifEmpty([])
    //   )
    //
    //   mqcReport = multiqc.out.report.toList()
    // }
}

workflow.onComplete {
  NFTools.makeReports(workflow, params, summary, customRunName, mqcReport)
}
