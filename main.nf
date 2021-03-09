#!/usr/bin/env nextflow
/*
Copyright Institut Curie 2020
This software is a computer program whose purpose is to analyze high-throughput sequencing data.
You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.
Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data.
The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.

This script is based on the nf-core guidelines. See https://nf-co.re/ for more information
*/

import common.VeganTools
@BaseScript(VeganTools)
import groovy.transform.BaseScript

/*
================================================================================
                        project : EUCANCAN/vegan
================================================================================
Started February 2020.
--------------------------------------------------------------------------------
vegan: Variant calling pipeline for whole Exome and whole Genome sequencing cANcer data Pipeline.
  An open-source analysis pipeline to detect germline or somatic variants
  from whole genome or targeted sequencing
--------------------------------------------------------------------------------
 @Homepage
 https://gitlab.curie.fr/data-analysis/vegan
--------------------------------------------------------------------------------
 @Documentation
 https://gitlab.curie.fr/data-analysis/vegan/README.md
--------------------------------------------------------------------------------
*/


// TODO: implicit call in baseScript by overwriting setup ?
// Welcome message and setup specific to the workflow
// Initialize lintedParams and paramsWithUsage
welcome()

/*
================================================================================
                               CONFIGURATION VARIABLES
================================================================================
*/

// Use lintedParams as default params object
colors = generateLogColors(params.get("monochromeLogs", false) as Boolean)
paramsWithUsage = readParamsFromJsonSettings("${projectDir}/parameters.settings.json")
def params = lint(params, paramsWithUsage)

tools = params.tools
skipQC = params.skipQC
SVFilters = params.SVFilters
SNVFilters = params.SNVFilters
annotateTools = params.annotateTools

customRunName = checkRunName(workflow.runName, params.name)
step = getStep(params.samplePlan, params.step)
samplePlanPath = getPath(step, params.samplePlan, params.outDir)
samplePlanCh = getSamplePlan(samplePlanPath, params.singleEnd, params.reads)
samplePlanCheckCh = params.samplePlan ? Channel.fromPath(samplePlanPath) : Channel.empty()
(designCh, designCheckCh) = params.design ? [getDesign(params.design), Channel.fromPath(params.design)] : [Channel.empty(), Channel.empty()]

if (params.design){
  (genderMap, statusMap, pairMap) = extractInfos(designCh)
} else {
  log.info """\
========================================================================
  ${colors.greenBold}INFO${colors.reset}: No design file detected
${tools && (('manta' in tools) || ('haplotypecaller' in tools) || ('mutect2' in tools))? "  ${colors.redBold}WARNING${colors.reset}: You need a design file in order to use VC tools\n" : ''}\
  ${colors.yellowBold}Variant detection (SV/SNV) will be skipped
  Please set up a design file '--design' to run these steps${colors.reset}
========================================================================"""
  tools = []
}
/*
================================================================================
                               CHECKING REFERENCES
================================================================================
*/

// Initialize each reference with default values in params.genomes, catch the genome defined on the command line first
// if it was defined
params << [
  fasta: params.genome && !('annotate' in step) ? params.genomes[params.genome].fasta ?: null : null,
  gtf: params.genome ? params.genomes[params.genome].gtf ?: null : null,
  polyms: params.genome && !('annotate' in step) ? params.genomes[params.genome].polyms ?: null : null,
  acLoci: params.genome && 'ascat' in tools ? params.genomes[params.genome].acLoci ?: null : null,
  acLociGC: params.genome && 'ascat' in tools ? params.genomes[params.genome].acLociGC ?: null : null,
  bwaIndex: params.genome && params.genomes[params.genome].fasta && 'mapping' in step ? params.genomes[params.genome].bwaIndex ?: null : null,
  dbsnp: params.genome && ('mapping' in step || 'haplotypecaller' in tools || 'mutect2' in tools) ? params.genomes[params.genome].dbsnp ?: null : null,
  dbsnpIndex: params.genome && params.genomes[params.genome].dbsnp ? params.genomes[params.genome].dbsnpIndex ?: null : null,
  dict: params.genome && params.genomes[params.genome].fasta ? params.genomes[params.genome].dict ?: null : null,
  fastaFai: params.genome && params.genomes[params.genome].fasta ? params.genomes[params.genome].fastaFai ?: null : null,
  germlineResource: params.genome && 'mutect2' in tools ? params.genomes[params.genome].germlineResource ?: null : null,
  germlineResourceIndex: params.genome && params.genomes[params.genome].germlineResource ? params.genomes[params.genome].germlineResourceIndex ?: null : null,
  intervals: params.genome && !('annotate' in step) ? params.genomes[params.genome].intervals ?: null : null,
  knownIndels: params.genome && 'mapping' in step ? params.genomes[params.genome].knownIndels ?: null : null,
  knownIndelsIndex: params.genome && params.genomes[params.genome].knownIndels ? params.genomes[params.genome].knownIndelsIndex ?: null : null,
  snpeffDb: params.genome && 'snpeff' in tools ? params.genomes[params.genome].snpeffDb ?: null : null,
]

/*
================================================================================
                                SUMMARY
================================================================================
*/

// Header info
summary = [
  'Pipeline Release': workflow.revision ?: null,
  'Run Name': customRunName,
  'Max Resources': "${params.maxMemory} memory, ${params.maxCpus} cpus, ${params.maxTime} time per job",
  'Container': workflow.containerEngine && workflow.container ? "${workflow.containerEngine} - ${workflow.container}" : null,
  'Target BED': params.targetBED ?: null,
  'Step': step ?: null,
  'Tools': params.tools ? params.tools instanceof Collection ? params.tools.join(', ') : params.tools: null,
  'QC tools skip': params.skipQC ? params.skipQC instanceof Collection  ? params.skipQC.join(', ') : skipQC : null,
  'SV filters': params.SVFilters ? params.SVFilters instanceof Collection  ? params.SVFilters.join(', ') : params.SVFilters : null,
  'SNV filters': params.SNVFilters ? params.SNVFilters instanceof Collection  ? params.SNVFilters.join(', ') : params.SNVFilters : null,
  'Intervals': params.noIntervals && step != 'annotate' ? 'Do not use' : null,
  'GVCF': 'haplotypecaller' in tools ? params.noGVCF ? 'No' : 'Yes' : null,
  'Sequenced by': params.sequencingCenter ? params.sequencingCenter: null,
  'Panel of normals': params.pon && 'mutect2' in tools ? params.pon: null,
  'Save Genome Index': params.saveGenomeIndex ? 'Yes' : 'No',
  'Output dir': params.outDir,
  'Launch dir': workflow.launchDir,
  'Working dir': workflow.workDir,
  'Script dir': workflow.projectDir,
  'User': workflow.userName,
  'Genome': params.genome,
  'Fasta': params.fasta ?: null,
  'FastaFai': params.fastaFai ?: null,
  'Polyms': params.polyms ?: null,
  'Dict': params.dict ?: null,
  'BwaIndex': params.bwaIndex ?: null,
  'GermlineResource': params.germlineResource ?: null,
  'GermlineResourceIndex': params.germlineResourceIndex ?: null,
  'acLoci': params.acLoci ?: null,
  'acLociGC': params.acLociGC ?: null,
  'dbsnp': params.dbsnp ?: null,
  'dbsnpIndex': params.dbsnpIndex ?: null,
  'knownIndels': params.knownIndels ?: null,
  'knownIndelsIndex': params.knownIndelsIndex ?: null,
  'snpeffDb': params.snpeffDb ?: null,
  'snpEffCache': params.snpEffCache ?: null,
  'Config Profile': workflow.profile,
  'Config Description': params.configProfileDescription ?: null,
  'Config Contact': params.configProfileContact ?: null,
  'Config URL': params.configProfileUrl ?: null,
  'E-mail Address': params.email ?: null,
  'MultiQC maxsize': params.email ? params.maxMultiqcEmailFileSize: null,
].findAll{ it.value != null }

// Check the hostnames against configured profiles
checkHostname(params, workflow)


/*
================================================================================
                               INIT CHANNELS
================================================================================
*/

// Initialize channels based on params
acLociCh = params.acLoci && 'ascat' in tools ? Channel.value(file(params.acLoci)) : "null"
acLociGCCh = params.acLociGC && 'ascat' in tools ? Channel.value(file(params.acLociGC)) : "null"
dbsnpCh = params.dbsnp && ('mapping' in step || 'haplotypecaller' in tools || 'mutect2' in tools) ? Channel.value(file(params.dbsnp)) : "null"
fastaCh = params.fasta && !('annotate' in step) ? Channel.value(file(params.fasta)) : "null"
fastaFaiCh = params.fastaFai && !('annotate' in step) ? Channel.value(file(params.fastaFai)) : "null"
gtfCh = params.gtf ? Channel.value(file(params.gtf)) : "null"
polymsCh = params.polyms ? Channel.value(file(params.polyms)) : "null"
germlineResourceCh = params.germlineResource && 'mutect2' in tools ? Channel.value(file(params.germlineResource)) : "null"
intervalsCh = params.intervals && !params.noIntervals && !('annotate' in step) ? Channel.value(file(params.intervals)) : "null"
knownIndelsCh = params.knownIndels ? Channel.value(file(params.knownIndels)) : "null"
snpEffCacheCh = params.snpEffCache ? Channel.value(file(params.snpEffCache)) : "null"
snpeffDbCh = params.snpeffDb ? Channel.value(params.snpeffDb) : "null"

// Optional files, not defined within the params.genomes[params.genome] scope
ponCh = params.pon ? Channel.value(file(params.pon)) : "null"
targetBedCh = params.targetBED ? Channel.value(file(params.targetBED)) : "null"
// Print summary and genareta summary channel
workflowSummaryCh = summarize(params, summary, workflow)
metadataCh = params.metadata ? Channel.fromPath(params.metadata) : "null"
outputDocsCh = file("$projectDir/docs/output.md", checkIfExists: true)
outputDocsImagesCh = file("$projectDir/docs/images/", checkIfExists: true)

/*
================================================================================
                                BUILDING INDEXES
================================================================================
*/

// And then initialize channels based on params or indexes that were just built

process buildBWAindexes {
  label 'bwa'

  tag "${fasta}"

  publishDir params.outDir, mode: params.publishDirMode,
    saveAs: {params.saveGenomeIndex ? "reference_genome/BWAIndex/${it}" : null }

  input:
  file(fasta) from fastaCh

  output:
  file("${fasta}.*") into bwaIndexesCh
  file("v_bwa.txt") into bwaVersionCh

  when: !(params.bwaIndex) && params.fasta && 'mapping' in step

  script:
  """
  bwa index ${fasta}
  bwa &> v_bwa.txt 2>&1 || true
  """
}

bwaIndexCh = params.bwaIndex ? Channel.value(file(params.bwaIndex)) : bwaIndexesCh

process buildDict {
  label 'gatk'

  tag "${fasta}"

  publishDir params.outDir, mode: params.publishDirMode,
    saveAs: {params.saveGenomeIndex ? "reference_genome/${it}" : null }

  input:
  file(fasta) from fastaCh

  output:
  file("${fasta.baseName}.dict") into dictBuiltCh

  when: !(params.dict) && params.fasta && !('annotate' in step)

  script:
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \
      CreateSequenceDictionary \
      --REFERENCE ${fasta} \
      --OUTPUT ${fasta.baseName}.dict
  """
}

dictCh = params.dict ? Channel.value(file(params.dict)) : dictBuiltCh

process buildFastaFai {
  label 'samtools'

  tag "${fasta}"

  publishDir params.outDir, mode: params.publishDirMode,
    saveAs: {params.saveGenomeIndex ? "reference_genome/${it}" : null }

  input:
  file(fasta) from fastaCh

  output:
  file("${fasta}.fai") into fastaFaiBuiltCh

  when: !(params.fastaFai) && params.fasta && !('annotate' in step)

  script:
  """
  samtools faidx ${fasta}
  """
}

fastaFaiCh = params.fastaFai ? Channel.value(file(params.fastaFai)) : fastaFaiBuiltCh

process buildDbsnpIndex {
  label 'tabix'

  tag "${dbsnp}"

  publishDir params.outDir, mode: params.publishDirMode,
    saveAs: {params.saveGenomeIndex ? "reference_genome/${it}" : null }

  input:
  file(dbsnp) from dbsnpCh

  output:
  file("${dbsnp}.tbi") into dbsnpIndexBuiltCh

  when: !(params.dbsnpIndex) && params.dbsnp && ('mapping' in step || 'haplotypecaller' in tools || 'mutect2' in tools)

  script:
  """
  tabix -p vcf ${dbsnp}
  """
}

dbsnpIndexCh = params.dbsnp ? params.dbsnpIndex ? Channel.value(file(params.dbsnpIndex)) : dbsnpIndexBuiltCh : "null"

process buildGermlineResourceIndex {
  label 'tabix'

  tag "${germlineResource}"

  publishDir params.outDir, mode: params.publishDirMode,
    saveAs: {params.saveGenomeIndex ? "reference_genome/${it}" : null }

  input:
  file(germlineResource) from germlineResourceCh

  output:
  file("${germlineResource}.tbi") into germlineResourceIndexBuiltCh

  when: !(params.germlineResourceIndex) && params.germlineResource && 'mutect2' in tools

  script:
  """
  tabix -p vcf ${germlineResource}
  """
}

germlineResourceIndexCh = params.germlineResource ? params.germlineResourceIndex ? Channel.value(file(params.germlineResourceIndex)) : germlineResourceIndexBuiltCh : "null"

process buildKnownIndelsIndex {
  label 'tabix'

  tag "${knownIndels}"

  publishDir params.outDir, mode: params.publishDirMode,
    saveAs: {params.saveGenomeIndex ? "reference_genome/${it}" : null }

  input:
  file(knownIndels) from knownIndelsCh

  output:
  file("${knownIndels}.tbi") into knownIndelsIndexBuiltCh

  when: !(params.knownIndelsIndex) && params.knownIndels && 'mapping' in step

  script:
  """
  tabix -p vcf ${knownIndels}
  """
}

knownIndelsIndexCh = params.knownIndels ? params.knownIndelsIndex ? Channel.value(file(params.knownIndelsIndex)) : knownIndelsIndexBuiltCh.collect() : "null"

process buildPonIndex {
  label 'tabix'

  tag "${pon}"

  publishDir params.outDir, mode: params.publishDirMode,
    saveAs: {params.saveGenomeIndex ? "reference_genome/${it}" : null }

  input:
  file(pon) from ponCh

  output:
  file("${pon}.tbi") into ponIndexBuiltCh

  when: !(params.ponIndex) && params.pon && 'mutect2' in tools

  script:
  """
  tabix -p vcf ${pon}
  """
}

process buildIntervals {
  label 'onlyLinux'

  tag "${fastaFai}"

  publishDir params.outDir, mode: params.publishDirMode,
    saveAs: {params.saveGenomeIndex ? "reference_genome/${it}" : null }

  input:
  file(fastaFai) from fastaFaiCh

  output:
  file("${fastaFai.baseName}.bed") into intervalBuiltCh

  when: !(params.intervals) && !('annotate' in step) && !(params.noIntervals)

  script:
  """
  awk -v FS='\t' -v OFS='\t' '{ print \$1, \"0\", \$2 }' ${fastaFai} > ${fastaFai.baseName}.bed
  """
}

intervalsCh = params.noIntervals ? "null" : params.intervals && !('annotate' in step) ? Channel.value(file(params.intervals)) : intervalBuiltCh


/*
================================================================================
                                  PREPROCESSING
================================================================================
*/

// STEP 0: CREATING INTERVALS FOR PARALLELIZATION (PREPROCESSING AND VARIANT CALLING)

process createIntervalBeds {
  label 'onlyLinux'

  tag "${intervals.fileName}"

  input:
  file(intervals) from intervalsCh

  output:
  file('*.bed') into bedIntervalsCh mode flatten

  when: (!params.noIntervals) && step != 'annotate'

  script:
  // If the interval file is BED format, the fifth column is interpreted to
  // contain runtime estimates, which is then used to combine short-running jobs
  if (hasExtension(intervals, "bed"))
    """
    awk -v FS="\t" '{
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
    grep -v '^@' ${intervals} | awk -v FS="\t" '{
      name = sprintf("%s_%d-%d", \$1, \$2, \$3);
      printf("%s\\t%d\\t%d\\n", \$1, \$2-1, \$3) > name ".bed"
    }'
    """
  else
    """
    awk -v FS="[:-]" '{
      name = sprintf("%s_%d-%d", \$1, \$2, \$3);
      printf("%s\\t%d\\t%d\\n", \$1, \$2-1, \$3) > name ".bed"
    }' ${intervals}
    """
}

bedIntervalsCh = bedIntervalsCh
  .map { intervalFile ->
    duration = 0.0
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

bedIntervalsCh = bedIntervalsCh.dump(tag:'bedintervals')

if (params.noIntervals && step != 'annotate') {file("${params.outDir}/noIntervals.bed").text = "noIntervals\n"; bedIntervalsCh = Channel.from(file("${params.outDir}/noIntervals.bed"))}

(intBaseRecalibratorCh, intApplyBQSRCh, intHaplotypeCallerCh, bedIntervalsCh) = bedIntervalsCh.into(5)

// PREPARING CHANNELS FOR PREPROCESSING AND QC

// TODO: not working
//(inputBamCh, inputPairReadsCh) = step == 'mapping' ? forkMappingSamplePlan(samplePlanCh) : [Channel.empty(), Channel.empty()]
if (step == "mapping") {
  runIds = [:]

  samplePlanCh.map {
    // Sample Plan
    // platformID - biologicalName - fastq1 - fastq2
    // Design
    // platformIdNormal - platformIdTumor - pairName - sex
    runIds[it[0]] = runIds.containsKey(it[0]) ? runIds[it[0]] + 1 : 0
    // sampleId, sampleName, runID, inputFiles
    // runId = platformId_XX
    return it[0,1] + [[it[0], runIds[it[0]].toString()].join("_")] + [it[2..-1]]
  }.branch {
     bamCh: it[3][0] =~ /.*bam$/
     pairCh: it[3][0] =~ /.*(fastq.gz|fq.gz|fastq|fq)$/
  }.set { samplePlanForks }
  (inputBamCh, inputPairReadsCh) = [samplePlanForks.bamCh, samplePlanForks.pairCh]
} else (inputBamCh, inputPairReadsCh) = [Channel.empty(), Channel.empty()]

// TODO: chek if splitFastq works
// I guess it does not work for uBam ?
if (params.splitFastq){
  inputPairReadsCh = inputPairReadsCh
  // newly splitfastq are named based on split, so the name is easier to catch
     .splitFastq(by: params.splitFastq, compress:true, file:"split", pe:true)
     .map {sampleId, sampleName, runId, reads ->
       // The split fastq read1 is the 4th element (indexed 3) its name is split_3
       // The split fastq read2's name is split_4
       // It's followed by which split it's acutally based on the mother fastq file
       // Index start at 1
       // Extracting the index to get a new IdRun
       splitIndex = reads[0].fileName.toString().minus("split_3.").minus(".gz")
       newIdRun = runId + "_" + splitIndex
       // Giving the files a new nice name
       newReads1 = file("${sampleId}_${newIdRun}_R1.fastq.gz")
       newReads2 = file("${sampleId}_${newIdRun}_R2.fastq.gz")
       [sampleId, sampleName, newIdRun, [newReads1, newReads2]]}
}


/*
================================================================================
                           RAW DATA QUALITY CHECK
================================================================================
*/


// Removing inputFile2 wich is null in case of uBAM
// TODO - should be the same for singleEnd data
inputBamCh = inputBamCh.map {
    sampleId, sampleName, runID, inputFiles ->
    [sampleId, sampleName, runID, inputFiles[0]]
}
(inputBamCh, inputBamFastQCCh) = inputBamCh.into(2)
(inputPairReadsCh, inputPairReadsFastQC) = inputPairReadsCh.into(2)

inputPairReadsCh = inputPairReadsCh.dump(tag:'inputPairReadsCh')

/*
 * FastQC
 */

process fastQC {
  label 'fastqc'
  label 'lowCpu'

  tag "${sampleId}"

  publishDir "${params.outDir}/fastqc", mode: params.publishDirMode, overwrite: true

  input:
  tuple val(sampleId),
    val(sampleName),
    val(runId),
    file(reads) from inputPairReadsFastQC.mix(inputBamFastQCCh)

  output:
  file("*.{html,zip}") into fastqcReportCh
  file("v_fastqc.txt") into fastqcVersionCh

  when: !('fastqc' in skipQC)

  script:
  """
  fastqc -t ${task.cpus} -q ${reads}
  fastqc --version > v_fastqc.txt
  """
}

/*
================================================================================
                              READS MAPPING
================================================================================
*/


/*
 * ALIGN READS TO REFERENCE GENOME WITH BWA-MEM
 */

inputPairReadsCh = inputPairReadsCh.mix(inputBamCh)
inputPairReadsCh = inputPairReadsCh.dump(tag:'INPUT MAP READS')

process mapReads {
  label 'gatkBwaSamtools'
  label 'highCpu'
  label 'extraMem'

  tag "${sampleId}"

  publishDir "${params.outDir}/mapping", mode: params.publishDirMode,
               saveAs: {filename ->  if (params.saveAlignedIntermediates) filename}

  input:
  tuple val(sampleId),
    val(sampleName),
    val(runId),
    file(inputFile) from inputPairReadsCh
  file(bwaIndex) from bwaIndexCh
  file(fasta) from fastaCh
  file(fastaFai) from fastaFaiCh

  output:
  tuple val(sampleId),
    val(sampleName),
    val(runId),
    file("${sampleId}.bam") into bamMappedCh
  tuple val(sampleId),
    val("${sampleId}_${runId}"),
    file("${sampleId}.bam") into bamMappedBamQCCh
  file('v_samtools.txt') into samtoolsMapReadsVersionCh

  script:
  // -K is an hidden option, used to fix the number of reads processed by bwa mem
  // Chunk size can affect bwa results, if not specified,
  // the number of threads can change which can give not deterministic result.
  // cf https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md
  // and https://github.com/gatk-workflows/gatk4-data-processing/blob/8ffa26ff4580df4ac3a5aa9e272a4ff6bab44ba2/processing-for-variant-discovery-gatk4.b37.wgs.inputs.json#L29
  CN = params.sequencingCenter ? "CN:${params.sequencingCenter}\\t" : ""
  readGroup = "@RG\\tID:${runId}\\t${CN}PU:${runId}\\tSM:${sampleId}\\tLB:${sampleId}\\tPL:illumina"
  // adjust mismatch penalty for tumor samples
  //status = statusMap[sampleId]
  //extra = status == 1 ? "-B 3" : ""
  convertToFastq = hasExtension(inputFile[0], "bam") ? "gatk --java-options -Xmx${task.memory.toGiga()}g SamToFastq --INPUT=${inputFile[0]} --FASTQ=/dev/stdout --INTERLEAVE=true --NON_PF=true | \\" : ""
  input = hasExtension(inputFile[0], "bam") ? "-p /dev/stdin - 2> >(tee ${inputFile[0]}.bwa.stderr.log >&2)" : "${inputFile[0]} ${inputFile[1]}"
  """
  ${convertToFastq}
  bwa mem ${params.bwaOptions} -R \"${readGroup}\" -t ${task.cpus} ${fasta} \
  ${input} | \
  samtools sort --threads ${task.cpus} - > ${sampleId}.bam
  samtools --version &> v_samtools.txt 2>&1 || true
  """
}

// Sort BAM whether they are standalone or should be merged
singleBamCh = Channel.create()
multipleBamCh = Channel.create()
bamMappedCh.groupTuple(by:[0, 1])
  .branch {
    singleCh: it[2].size() == 1
    multipleCh: it[2].size() > 1
  }.set { bamMappedForks }
(singleBamCh, multipleBamCh) = [bamMappedForks.singleCh, bamMappedForks.multipleCh]

singleBamCh = singleBamCh.map {
  sampleId, sampleName, runId, bam ->
    [sampleId, sampleName, bam]
}

singleBamCh = singleBamCh.dump(tag:'sbams')
multipleBamCh = multipleBamCh
  .map{it -> [it[0], it[1], it[3][0]]}
  .groupTuple()
  .dump(tag:'mbams')


/*
 * MERGING BAM FROM MULTIPLE LANES
 */

process mergeBamMapped {
  label 'samtools'
  label 'highCpu'

  tag "${sampleId}"

  publishDir "${params.outDir}/mapping", mode: params.publishDirMode,
               saveAs: {filename ->  if (params.saveAlignedIntermediates) filename}

  input:
  tuple val(sampleId),
    val(sampleName),
    val(bams) from multipleBamCh

  output:
  tuple val(sampleId),
    val(sampleName),
    file("*_merged.bam") into mergedBamCh
  file('v_samtools.txt') into samtoolsMergeBamMappedVersionCh

  script:
  """
  samtools merge --threads ${task.cpus} ${sampleId}_merged.bam ${bams[0]}
  samtools --version &> v_samtools.txt 2>&1 || true
  """
}

mergedBamCh = mergedBamCh.mix(singleBamCh).dump(tag:'bams')
(mergedBamCh, mergedBamPreseqCh, mergedBamToStatsCh, mergedBamToIndexCh) = mergedBamCh.into(4)


/*
 * INDEX ALIGNED BAM FILE
 */

process indexBamFile {
  label 'samtools'
  label 'minCpu'

  tag "${sampleId}"

  publishDir "${params.outDir}/mapping", mode: params.publishDirMode,
               saveAs: {filename ->  if (params.saveAlignedIntermediates) filename}

  input:
  tuple val(sampleId),
    val(sampleName),
    file(bam) from mergedBamToIndexCh

  output:
  tuple val(sampleId),
    val(sampleName),
    file(bam),
    file("*.bai") into indexedBamCh
  file('v_samtools.txt') into samtoolsIndexBamFileVersionCh

  script:
  """
  samtools index ${bam}
  samtools --version &> v_samtools.txt 2>&1 || true
  """
}


/*
 * BWA-MEM MAPPING STATISTICS
 */

process bamStats {
  label 'samtools'
  label 'lowCpu'

  tag {sampleId}

  publishDir "${params.outDir}/mapping/logs", mode: params.publishDirMode

  input:
  tuple val(sampleId),
    val(sampleName),
    file(bam) from mergedBamToStatsCh

  output:
  file("*_mappingstats.mqc") into bamStatsMqcCh
  file("*bwa.log") into bwaMqcCh
  file('v_samtools.txt') into samtoolsMappingStatsVersionCh

  script:
  """
  getBWAstats.sh -i ${bam} -p ${task.cpus} > ${sampleId}_bwa.log
  aligned="\$(samtools view -@ ${task.cpus} -F 0x100 -F 0x4 -F 0x800 -c ${bam})"
  hqbam="\$(samtools view -@ ${task.cpus} -F 0x100 -F 0x800 -F 0x4 -q 20 -c ${bam})"
  lqbam="\$((\$aligned - \$hqbam))"
  echo -e "Mapped,\${aligned}" > ${sampleId}_mappingstats.mqc
  echo -e "HighQual,\${hqbam}" >> ${sampleId}_mappingstats.mqc
  echo -e "LowQual,\${lqbam}" >> ${sampleId}_mappingstats.mqc
  samtools --version &> v_samtools.txt 2>&1 || true
  """
}


/*
 * Preseq (before alignment filtering)
 */

process preseq {
  label 'preseq'
  label 'lowCpu'
  label 'medMem'

  tag "${sampleID}"

  publishDir "${params.outDir}/preseq", mode: params.publishDirMode

  when:
  !params.skipPreseq

  input:
  tuple val(sampleID),
    val(sampleName),
    file(bam) from mergedBamPreseqCh

  output:
  file("*.ccurve.txt") into preseqStatsCh
  file("v_preseq.txt") into preseqVersionCh

  script:
  defectMode = params.preseqDefect ? '-D' : ''
  """
  preseq &> v_preseq.txt
  preseq lc_extrap -v $defectMode -output ${bam.baseName}.ccurve.txt -bam ${bam}
  """
}


/*
================================================================================
                                  FILTERING
================================================================================
*/


/*
 * Duplicates - sambamba
 */

process markDuplicates {
  label 'sambamba'
  label 'highCpu'
  label 'highMem'

  tag "${sampleId}"

  publishDir "${params.outDir}/markDuplicates", mode: params.publishDirMode,
              saveAs: {filename -> if (params.saveAlignedIntermediates) filename}

  input:
  tuple val(sampleId),
    val(sampleName),
    file(bam) from mergedBamCh

  output:
  tuple val(sampleId),
    val(sampleName),
    file("${sampleId}.md.bam"),
    file("${sampleId}.md.bam.bai") into duplicateMarkedBamsCh,mdBamPolymCh
  file ("${sampleId}.md.bam.metrics") into markDuplicatesReportCh

  script:
  """
  sambamba markdup --nthreads ${task.cpus} --tmpdir . ${bam} ${sampleId}.md.bam
  sambamba flagstat --nthreads ${task.cpus} ${sampleId}.md.bam > ${sampleId}.md.bam.metrics
  """
}

/*
 * BAM on Target
 */

bamsToTargetCh = params.targetBED ? duplicateMarkedBamsCh : Channel.empty()

process bamOnTarget {
  label 'bedtools'
  label 'minCpu'
  label 'medMem'

  tag {sampleId}

  publishDir "${params.outDir}/onTarget", mode: params.publishDirMode,
               saveAs: {filename ->  if (params.saveAlignedIntermediates) filename}

  when:
  params.targetBED

  input:
  tuple val(sampleId),
    val(sampleName),
    file(bam),
    file(bai) from bamsToTargetCh
  file(targetBED) from targetBedCh

  output:
  tuple val(sampleId),
    val(sampleName),
    file("*_onTarget.bam"),
    file("*_onTarget.bam.bai") into onTargetBamsCh
  file("${bam.baseName}_onTarget.bam.metrics") into onTargetReportCh

  script:
  """
  intersectBed -abam ${bam} -b ${targetBED} > ${bam.baseName}_onTarget.bam
  samtools index ${bam.baseName}_onTarget.bam
  samtools flagstat ${bam.baseName}_onTarget.bam > ${bam.baseName}_onTarget.bam.metrics
  """
}

procBamsCh = params.targetBED ? onTargetBamsCh : duplicateMarkedBamsCh

/*
 * FILTER ALIGNED BAM FILE FOR SNV/SV
 */

if (('manta' in tools) && ('ascat' in tools || 'haplotypecaller' in tools || 'mutect2' in tools)){
  // Duplicates the channel for SV and SNV filtering
  procBamsCh = procBamsCh.flatMap { it -> [it + 'SV', it + 'SNV']}
}else if ( ('manta' in tools) && !('ascat' in tools || 'haplotypecaller' in tools || 'mutect2' in tools)){
  // SV only
  procBamsCh = procBamsCh.flatMap { it -> [it + 'SV']}
}else{
  // SNV only if no design or just SNV
  procBamsCh = procBamsCh.flatMap { it -> [it + 'SNV']}
}
procBamsCh = procBamsCh.dump(tag:'procBamsCh')

process bamFiltering {
  label 'samtools'
  label 'medCpu'
  tag "${sampleId}-${vCType}"

  publishDir "${params.outDir}/mapping/filtering", mode: params.publishDirMode

  input:
  tuple val(sampleId),
    val(sampleName),
    file(bam),
    file(bai),
    val(vCType) from procBamsCh

  output:
  tuple val(sampleId),
    val(sampleName),
    val(vCType),
    file("${sampleId}.filtered.${vCType}.bam"),
    file("${sampleId}.filtered.${vCType}.bam.bai") into filteredBamCh, filteredBamQCCh
  file("${sampleId}.filtered.${vCType}.idxstats") into bamFilterReportCh
  file('v_samtools.txt') into samtoolsBamFilterVersionCh

  script:
  dupParams = (vCType == 'SNV' && 'markduplicates' in SNVFilters) | (vCType == 'SV' && 'markduplicates' in SVFilters) ? "-F 0x0400" : ""
  mapqParams = (vCType == 'SNV' && 'mapq' in SNVFilters) | (vCType == 'SV' && 'mapq' in SVFilters) && (params.mapQual > 0) ? "-q ${params.mapQual}" : ""
  // Remove singletons and keep paired reads + Delete secondary and not primary alignment (0x100)
  uniqParams =  (vCType == 'SNV' && 'uniq' in SNVFilters) | (vCType == 'SV' && 'uniq' in SVFilters) ? "-F 0x004 -F 0x0008 -f 0x001 -F 0x100 -F 0x800" :  ""
  uniqFilter = (vCType == 'SNV' && 'uniq' in SNVFilters) | (vCType == 'SV' && 'uniq' in SVFilters) ? "| grep -v -e \\\"XA:Z:\\\" -e \\\"SA:Z:\\\" | samtools view -b -" : "| samtools view -b -"
  """
  samtools view -h -@ ${task.cpus} ${uniqParams} ${dupParams} ${mapqParams} ${bam} ${uniqFilter} > ${sampleId}.filtered.${vCType}.bam
  samtools index ${sampleId}.filtered.${vCType}.bam
  samtools flagstat ${sampleId}.filtered.${vCType}.bam > ${sampleId}.filtered.${vCType}.flagstats
  samtools idxstats ${sampleId}.filtered.${vCType}.bam > ${sampleId}.filtered.${vCType}.idxstats
  samtools stats ${sampleId}.filtered.${vCType}.bam > ${sampleId}.filtered.${vCType}.stats
  samtools --version &> v_samtools.txt 2>&1 || true
  """
}


/*
================================================================================
                                  QUALITY CHECK
================================================================================
*/


filteredBamCh = filteredBamCh.dump(tag:'filteredBamCh')
// Run the QC on the SNV bam only if available - on the SV otherwise
if ( ('manta' in tools) && !('ascat' in tools || 'haplotypecaller' in tools || 'mutect2' in tools)){
  filteredBamQCCh
    .filter { it[2] == 'SV' }
    .dump(tag:'qcbams')
    .into {bamQualimapCh; bamInsertSizeCh; bamMosdepthCh; bamGeneCovCh; bamWGSmetricsCh }
} else {
  filteredBamQCCh
    .filter { it[2] == 'SNV' }
    .dump(tag:'qcbams')
    .into {bamQualimapCh; bamInsertSizeCh; bamMosdepthCh; bamGeneCovCh; bamWGSmetricsCh }
}

/*
 * INSERT SIZE
 */

process getFragmentSize {
  label 'picard'
  label 'medCpu'
  label 'medMem'

  tag "${sampleId}"

  publishDir path: "${params.outDir}/fragSize", mode: params.publishDirMode

  when:
  !params.singleEnd

  input:
  tuple val(sampleId),
    val(sampleName),
    val(vCType),
    file(bam),
    file(bai) from bamInsertSizeCh

  output:
  file("*.{pdf,txt}") into fragmentSizeCh

  script:
  """
  picard CollectInsertSizeMetrics \
      I=${bam} \
      O=${bam.baseName}_insert_size_metrics.txt \
      H=${bam.baseName}_insert_size_hist.pdf \
      M=0.5
  """
}

/*
 * SEQUENCING DEPTH
 */

process getSeqDepth {
  label 'mosdepth'
  label 'medCpu'
  label 'medMem'

  tag "${sampleId}"

  publishDir path: "${params.outDir}/depth", mode: params.publishDirMode

  input:
  tuple val(sampleId),
    val(sampleName),
    val(vCType),
    file(bam),
    file(bai) from bamMosdepthCh
  file(bed) from targetBedCh

  output:
  file("*.txt") into mosdepthOutputCh
  file("*{.bed.gz,.bed.gz.csi}") into mosdepthBedOutputCh

  script:
  bedCmd = params.targetBED ? "--by ${bed}" : ''
  """
  mosdepth -t ${task.cpus} -n --quantize 0:1:10:50:100: ${bedCmd} ${bam.baseName} ${bam}
  """
}

/*
 * GENES COVERAGE
 */

process prepareExonInfo {
  label 'bedtools'
  label 'medCpu'
  label 'medMem'

  input:
  file(gtf) from gtfCh
  file(bed) from targetBedCh

  output:
  file("*exon.bed") into exonBedCh

  script:
  targetCmd = params.targetBED ? " | intersectBed -a stdin -b ${bed} ": ''
  """
  awk -F"\t" -v type='gene_id' 'BEGIN{OFS="\t"} \$3=="exon" {split(\$9,annot,";");for(i=1;i<=length(annot);i++){if (annot[i]~type){anntype=annot[i]}} print \$1,\$4-1,\$5,anntype}' ${gtf} | sed -e 's/gene_id//' -e 's/"//g' | sort -u -k1,1V -k2,2n ${targetCmd} > ${gtf.baseName}_exon.bed
  """
}

process genesCoverage {
  label 'mosdepth'
  label 'medCpu'
  label 'medMem'
  publishDir path: "${params.outDir}/depth", mode: params.publishDirMode

  tag "${sampleId}"

  input:
  tuple val(sampleId),
    val(sampleName),
    val(vCType),
    file(bam),
    file(bai) from bamGeneCovCh
  file(exon) from exonBedCh

  output:
  file("*.mqc") into geneCovMqc
  file("*.pdf") into geneCovOutput

  script:
  """
  mosdepth -n -t ${task.cpus} --by ${exon} ${sampleId}.genecov ${bam}
  apGeneCov.r --cov ${sampleId}.genecov.regions.bed.gz --oprefix ${sampleId}_covdensity
  """
}

/*
 * MATES OVERLAP
 */

process getWGSmetrics {
  tag "${sampleId}"
  label 'picard'
  label 'medCpu'
  label 'medMem'

  publishDir path: "${params.outDir}/WGSmetrics", mode: params.publishDirMode

  input:
  tuple val(sampleId),
    val(sampleName),
    val(vCType),
    file(bam),
    file(bai) from bamWGSmetricsCh
  file(reference) from fastaCh
  file(dict) from dictCh
  file(bed) from targetBedCh

  output:
  file("*.txt") into wgsMetricsOutputCh

  script:
  memOption = "\"-Xms" +  (task.memory.toGiga() / 2).trunc() + "g -Xmx" + (task.memory.toGiga() - 1) + "g\""
  bedTointerCmd = params.targetBED ? "picard BedToIntervalList I=${bed} O=intervals.bed SD=${dict}":""
  bedCmd = params.targetBED ? "INTERVALS=intervals.bed" : ""

  """
  ${bedTointerCmd}
  picard ${memOption} CollectWgsMetrics \
       USE_FAST_ALGORITHM=true \
       I=${bam} \
       O=${bam.baseName}_collect_wgs_metrics.txt \
       R=${reference} \
       ${bedCmd}
  """
}

/*
 * Identito monitoring
 */

mdBamPolymCh = mdBamPolymCh.dump(tag:'polymCh')

process getPolym {
  label 'lowCpu'
  label 'lowMem'
  label 'polym'

  publishDir "${params.outDir}/identito", mode: params.publishDirMode

  when: !(params.skipIdentito)

  input:
  file(fasta) from fastaCh
  file(fastaFai) from fastaFaiCh
  file(polyms) from polymsCh
  tuple val(sampleId),
    val(sampleName),
    file("${sampleId}.md.bam"),
    file("${sampleId}.md.bam.bai") from mdBamPolymCh

  output:
  file("${sampleId}_matrix.tsv") into clustPolymCh

  script:
  """
  bcftools mpileup -R ${polyms} -f ${fasta} -x -A -B -q 20 -I -Q 0 -d 1000 --annotate FORMAT/DP,FORMAT/AD ${sampleId}.md.bam > ${sampleId}_bcftools.vcf
  SnpSift extractFields -e "."  -s ";" ${sampleId}_bcftools.vcf CHROM POS REF ALT GEN[*].DP GEN[*].AD > ${sampleId}_bcftools.tsv
  apComputePolym.R ${sampleId}_bcftools.tsv ${sampleId}_matrix.tsv ${sampleId} ${polyms}
  """
}

clustPolymCh = clustPolymCh.dump(tag:'clustPolymCh')

process computePolym {
  label 'lowCpu'
  label 'lowMem'
  label 'polym'

  publishDir "${params.outDir}/identito", mode: params.publishDirMode

  when: !(params.skipIdentito)

  input:
  file(matrix) from clustPolymCh.collect()

  output:
  file("*") into clustPolymResultsCh

  script:
  """
  cat *matrix.tsv |awk 'NR==1{print \$0}' > clust_mat.tsv

  for in_matrix in ${matrix}
  do
  awk 'NR>1{print \$0}' \$in_matrix >> clust_mat.tsv
  done

  apComputeClust.R clust_mat.tsv . clustering_plot
  """
}


/*
================================================================================
                          RECALIBRATING
================================================================================
*/

filteredBamCh
  .branch { it ->
    snvCh: it[2] == 'SNV'
    svCh: it[2] == 'SV'
    otherCh: true
  }.set { filteredBamForks }
(filteredSNVBamsCh, filteredSVBamsCh, filteredOtherBamsCh) = [filteredBamForks.snvCh, filteredBamForks.svCh, filteredBamForks.otherCh]

filteredSNVBamsCh
  .into {bamBaseRecalibratorCh; bamBaseRecalibratorToJoinCh}
bamBaseRecalibratorCh = bamBaseRecalibratorCh.combine(intBaseRecalibratorCh)


/*
 * CREATING RECALIBRATION TABLES
 */

process baseRecalibrator {
  label 'gatk'
  label 'minCpu'

  tag "${sampleId}"

  input:
  tuple val(sampleId),
    val(sampleName),
    val(vCType),
    file(bam),
    file(bai),
    file(intervalBed) from bamBaseRecalibratorCh.dump(tag:'bamBaseRecalibratorCh')
  file(dbsnp) from dbsnpCh
  file(dbsnpIndex) from dbsnpIndexCh
  file(fasta) from fastaCh
  file(fastaFai) from fastaFaiCh
  file(dict) from dictCh
  file(knownIndels) from knownIndelsCh
  file(knownIndelsIndex) from knownIndelsIndexCh

  output:
  tuple val(sampleId),
    val(sampleName),
    val(vCType),
    file("${prefix}${sampleId}.recal.table") into tableGatherBQSRReportsCh
  tuple val(sampleId),
    val(sampleName),
    val(vCType) into recalTableTSVnoIntCh

  when: 'haplotypecaller' in tools || 'mutect2' in tools

  script:
  dbsnpOptions = params.dbsnp ? "--known-sites ${dbsnp}" : ""
  knownOptions = params.knownIndels ? knownIndels.collect{"--known-sites ${it}"}.join(' ') : ""
  prefix = params.noIntervals ? "${vCType}_" : "${intervalBed.baseName}_${vCType}"
  intervalsOptions = params.noIntervals ? "" : "-L ${intervalBed}"
  // TODO: --use-original-qualities ???
  """
  gatk --java-options -Xmx${task.memory.toGiga()}g \
      BaseRecalibrator \
      -I ${bam} \
      -O ${prefix}${sampleId}.recal.table \
      --tmp-dir ${params.baseRecalibratorOpts} \
      -R ${fasta} \
      ${intervalsOptions} \
      ${dbsnpOptions} \
      ${knownOptions} \
      --verbosity INFO
  """
}

/*
 * MERGE BQSR TABLES PER INERVALS
 */

if (!params.noIntervals) {
  tableGatherBQSRReportsCh = tableGatherBQSRReportsCh.groupTuple(by:[0, 1, 2])
}else{
  (tableGatherBQSRReportsCh, recalTableCh) = tableGatherBQSRReportsCh.into(2)
}

process gatherBQSRReports {
  label 'gatk'
  label 'lowMem'
  label 'lowCpu'

  tag "${sampleId}"

  publishDir "${params.outDir}/Preprocessing/${sampleId}/DuplicateMarked", mode: params.publishDirMode, overwrite: false

  input:
  tuple val(sampleId),
    val(sampleName),
    val(vCType),
    file(recal) from tableGatherBQSRReportsCh

  output:
  tuple val(sampleId),
    val(sampleName),
    val(vCType),
    file("${prefix}${sampleId}.recal.table") into recalTableCh
  tuple val(sampleId),
    val(sampleName),
    val(vCType) into recalTableCSVCh

  when: !(params.noIntervals)

  script:
  input = recal.collect{"-I ${it}"}.join(' ')
  prefix = "${vCType}_"
  """
  gatk --java-options -Xmx${task.memory.toGiga()}g \
      GatherBQSRReports \
      ${input} \
      -O ${prefix}${sampleId}.recal.table \
  """
}

// Create TSV files to restart from this step
recalTableCSVCh = recalTableCSVCh.mix(recalTableTSVnoIntCh)
recalTableCSVCh = recalTableCSVCh.dump(tag: 'recalTableCSVCh')
recalTableCSVCh.map { sampleId, sampleName, vCType ->
  bam = "${params.outDir}/Preprocessing/${sampleId}/DuplicateMarked/${sampleId}.md.bam"
  bai = "${params.outDir}/Preprocessing/${sampleId}/DuplicateMarked/${sampleId}.md.bai"
  recalTable = "${params.outDir}/Preprocessing/${sampleId}/DuplicateMarked/${vCType}_${sampleId}.recal.table"
  "${sampleId},${sampleName},${vCType},${bam},${bai},${recalTable}\n"
}.collectFile(
  name: 'samplePlan.recal.csv', sort: true, storeDir: "${params.outDir}/Preprocessing/CSV"
)

// TODO: find if generated files below are useful or not
//recalTableSampleTSVCh
//    .collectFile(storeDir: "${params.outDir}/Preprocessing/TSV/") {
//        sampleId, sampleName, vCType ->
//        status = statusMap[sampleId]
//        gender = genderMap[sampleId]
//        bam = "${params.outDir}/Preprocessing/${sampleName}/DuplicateMarked/${sampleName}.md.bam"
//        bai = "${params.outDir}/Preprocessing/${sampleName}/DuplicateMarked/${sampleName}.md.bai"
//        recalTable = "${params.outDir}/Preprocessing/${sampleName}/DuplicateMarked/${sampleName}.recal.table"
//        ["duplicateMarked_${sampleName}.tsv", "${sampleId}\t${vCType}\t${sampleName}\t${bam}\t${bai}\t${recalTable}\n"]
//}

bamApplyBQSRCh = step in 'recalibrate' ? samplePlanCh : bamBaseRecalibratorToJoinCh.join(recalTableCh, by:[0,1,2])
bamApplyBQSRCh = bamApplyBQSRCh.combine(intApplyBQSRCh)


/*
 * RECALIBRATING
 */

process applyBQSR {
  label 'gatk'
  label 'lowMem'
  label 'lowCpu'

  tag "${sampleId}-${sampleName}-${vCType}-${intervalBed.baseName}"

  input:
  tuple val(sampleId),
    val(sampleName),
    val(vCType),
    file(bam),
    file(bai),
    file(recalibrationReport),
    file(intervalBed) from bamApplyBQSRCh
  file(dict) from dictCh
  file(fasta) from fastaCh
  file(fastaFai) from fastaFaiCh

  output:
  tuple val(sampleId),
    val(sampleName),
    val(vCType),
    file("${prefix}${sampleId}.recal.bam") into bamMergeBamRecalCh
  file("v_gatk.txt") into gatkVersionCh

  script:
  prefix = params.noIntervals ? "${vCType}_noInterval_" : "${vCType}_${intervalBed.baseName}_"
  intervalsOptions = params.noIntervals ? "" : "-L ${intervalBed}"
  """
  gatk --java-options -Xmx${task.memory.toGiga()}g \
      ApplyBQSR \
      -R ${fasta} \
      --input ${bam} \
      --output ${prefix}${sampleId}.recal.bam \
      ${intervalsOptions} \
      --bqsr-recal-file ${recalibrationReport}
  gatk ApplyBQSR --help &> v_gatk.txt 2>&1 || true
  """
}

bamMergeBamRecalCh = bamMergeBamRecalCh.groupTuple(by:[0, 1, 2])
(bamMergeBamRecalCh, bamMergeBamRecalNoIntCh) = bamMergeBamRecalCh.into(2)


/*
 * MERGING THE RECALIBRATED BAM FILES
 */

process mergeAndIndexBamRecal {
  label 'samtools'
  label 'medCpu'

  tag "${sampleId}-${sampleName}-${vCType}"

  publishDir "${params.outDir}/Preprocessing/${sampleId}/Recalibrated", mode: params.publishDirMode

  input:
  tuple val(sampleId),
    val(sampleName),
    val(vCType),
    file(bam) from bamMergeBamRecalCh

  output:
  tuple val(sampleId),
    val(sampleName),
    val(vCType),
    file("*recal.bam"),
    file("*recal.bam.bai") into bamRecalCh
  tuple val(sampleId),
    val(sampleName),
    val(vCType) into bamRecalTSVCh
  file('v_samtools.txt') into samtoolsMergeBamRecalVersionCh

  when: !(params.noIntervals)

  script:
  """
  samtools merge --threads ${task.cpus} ${sampleId}_${vCType}_recal.bam ${bam}
  samtools index ${sampleId}_${vCType}_recal.bam
  samtools --version &> v_samtools.txt 2>&1 || true
  """
}


/*
 * INDEXING THE MERGED RECALIBRATED BAM FILES
 */

process indexBamRecal {
  label 'samtools'
  label 'lowCpu'

  tag "${sampleId}-${sampleName}-${vCType}"

  publishDir "${params.outDir}/Preprocessing/${sampleId}/Recalibrated", mode: params.publishDirMode

  input:
  tuple val(sampleId),
    val(sampleName),
    val(vCType),
    file(bam) from bamMergeBamRecalNoIntCh

  output:
  tuple val(sampleId),
    val(sampleName),
    val(vCType),
    file(bam),
    file("*bam.bai") into bamRecalNoIntCh
  tuple val(sampleId),
    val(sampleName),
    val(vCType) into bamRecalTSVnoIntCh
  file('v_samtools.txt') into samtoolsIndexBamRecalVersionCh

  when: params.noIntervals

  script:
  """
  samtools index ${bam}
  samtools --version &> v_samtools.txt 2>&1 || true
  """
}

bamRecalCh = bamRecalCh.mix(bamRecalNoIntCh)
//bamRecalQCCh = bamRecalQCCh.mix(bamRecalQCnoIntCh)
bamRecalTSVCh = bamRecalTSVCh.mix(bamRecalTSVnoIntCh)

//(bamRecalQualimapCh, bamRecalSamToolsStatsCh, bamRecalMosdepthCh, bamRecalInsertSizeCh) = bamRecalCh.into(4)
(bamRecalTSVCh, bamRecalSampleTSVCh) = bamRecalTSVCh.into(2)

// Creating a TSV file to restart from this step
bamRecalTSVCh.map { sampleId, sampleName, vCType ->
  bam = "${params.outDir}/Preprocessing/${sampleId}/Recalibrated/${sampleId}.recal.bam"
  bai = "${params.outDir}/Preprocessing/${sampleId}/Recalibrated/${sampleId}.recal.bam.bai"
  "${sampleId}\t${sampleName}\t${vCType}\t${bam}\t${bai}\n"
}.collectFile(
  name: 'recal.samplePlan.tsv', sort: true, storeDir: "${params.outDir}/Preprocessing/TSV"
)

// TODO: find if generated files below are useful or not
//bamRecalSampleTSVCh
//    .collectFile(storeDir: "${params.outDir}/Preprocessing/TSV") {
//        sampleId, sampleName ->
//        status = statusMap[sampleId]
//        gender = genderMap[sampleId]
//        bam = "${params.outDir}/Preprocessing/${sampleName}/Recalibrated/${sampleName}.recal.bam"
//        bai = "${params.outDir}/Preprocessing/${sampleName}/Recalibrated/${sampleName}.recal.bam.bai"
//        ["recalibrated_${sampleName}.tsv", "${sampleId}\t${gender}\t${status}\t${sampleName}\t${bam}\t${bai}\n"]
//}

// When no knownIndels for mapping, Channel bamRecalCh is indexedBamCh
// TODO: seems not suited to the actual layout, have to refactor line below (indexed bam are not filtered)
// bamRecalCh = (params.knownIndels && step == 'mapping') ? bamRecalCh : indexedBamCh.flatMap { it -> [it.plus(2, 'SV'), it.plus(2, 'SNV')]}


/*
================================================================================
                            VARIANT CALLING
================================================================================
*/

if (params.design){
  // When starting with variant calling, Channel bamRecalCh is samplePlanCh
  bamRecalCh = step in 'variantcalling' ? samplePlanCh : bamRecalCh

  // Here we have a recalibrated bam set
  // The TSV file is formatted like: "sampleId status sampleName vCType bamFile baiFile"
  // Manta will be run in Germline mode, or in Tumor mode depending on status
  // HaplotypeCaller will be run for Normal and Tumor samples

  // SV Bams
  (filteredSVBamsCh, bamMantaSingleCh) = filteredSVBamsCh.into(2)

  // SNV Bams
  (bamAscatCh, bamRecalAllCh, bamRecalAllTempCh) = bamRecalCh.into(4)
  bamHaplotypeCallerCh = bamRecalAllTempCh.combine(intHaplotypeCallerCh)
  //(bamAscatCh, bamRecalAllCh) = bamRecalAllCh.into(2)

  // separate BAM by status for somatic variant calling
  bamRecalAllCh = bamRecalAllCh.dump(tag: 'bamRecalAllCh')

  bamRecalAllCh.branch{
    normalCh: statusMap[it[0]] == 0
    tumorCh: statusMap[it[0]] == 1
  }.set { bamRecalAllForks }

  filteredSVBamsCh.branch{
    normalCh: statusMap[it[0]] == 0
    tumorCh: statusMap[it[0]] == 1
  }.set { filteredSVBamsForks }

  (bamRecalNormalCh, bamRecalTumorCh) = [bamRecalAllForks.normalCh, bamRecalAllForks.tumorCh]
  (bamSVNormalCh, bamSVTumorCh) = [filteredSVBamsForks.normalCh, filteredSVBamsForks.tumorCh]

  // Crossing Normal and Tumor to get a T/N pair for Somatic Variant Calling
  // Remapping channel to remove common key sampleId

  pairBamsSNVCh = bamRecalNormalCh.combine(bamRecalTumorCh)
  pairBamsSVCh = bamSVNormalCh.combine(bamSVTumorCh)
  pairBamsSNVCh = pairBamsSNVCh.dump(tag: 'pairBamsSNVCh')

  pairBamsSNVCh = pairBamsSNVCh.filter{ pairMap.containsKey([it[0], it[5]]) && it[2] == it[7] }
  pairBamsSVCh = pairBamsSVCh.filter{ pairMap.containsKey([it[0], it[5]]) && it[2] == it[7] }

  pairBamsSNVCh = pairBamsSNVCh.dump(tag: 'pairBamsSNVCh')
  pairBamsSVCh = pairBamsSVCh.dump(tag: 'pairBamsSVCh')

  // Manta,  Mutect2
  (pairBamCalculateContaminationCh, pairBamsSNVCh) = pairBamsSNVCh.into(3)
  (pairBamMantaCh, pairBamsSVCh) = pairBamsSVCh.into(3)

  intervalPairBamCh = pairBamsSNVCh.combine(bedIntervalsCh)

  // intervals for Mutect2 calls and pileups for Mutect2 filtering
  (pairBamMutect2Ch, pairBamPileupSummariesCh) = intervalPairBamCh.into(2)
  pairBamMutect2Ch = pairBamMutect2Ch.dump(tag: 'pairBamMutect2Ch')
} else {
  bamHaplotypeCallerCh = Channel.empty()
  pairBamPileupSummariesCh = Channel.empty()
  pairBamCalculateContaminationCh = Channel.empty()
  pairBamMutect2Ch = Channel.empty()
  bamMantaSingleCh = Channel.empty()
  pairBamMantaCh = Channel.empty()
  bamAscatCh = Channel.empty()
}

/*
================================================================================
                            SNV VARIANT CALLING
================================================================================
*/

// To speed Variant Callers up we are chopping the reference into smaller pieces
// Do variant calling by this intervals, and re-merge the VCFs


// STEP GATK HAPLOTYPECALLER.1
process haplotypeCaller {
  label 'gatk'
  label 'medMemSq'
  label 'lowCpu'

  tag "${sampleId}-${intervalBed.baseName}"

  input:
  tuple val(sampleId),
    val(sampleName),
    val(vCType),
    file(bam),
    file(bai),
    file(intervalBed) from bamHaplotypeCallerCh
  file(dbsnp) from dbsnpCh
  file(dbsnpIndex) from dbsnpIndexCh
  file(dict) from dictCh
  file(fasta) from fastaCh
  file(fastaFai) from fastaFaiCh

  output:
  tuple val("HaplotypeCallerGVCF"),
    val(sampleId),
    val(sampleName),
    file("${intervalBed.baseName}_${sampleId}.g.vcf") into gvcfHaplotypeCallerCh
  tuple val(sampleId),
    val(sampleName),
    file(intervalBed),
    file("${intervalBed.baseName}_${sampleId}.g.vcf") into gvcfGenotypeGVCFsCh

  when: 'haplotypecaller' in tools && vCType == 'SNV'

  script:
  intervalOpts = params.noIntervals ? "" : "-L ${intervalBed}"
  dbsnpOpts = params.dbsnp ? "--D ${dbsnp}" : ""
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g -Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
    HaplotypeCaller \
    -R ${fasta} \
    -I ${bam} \
    ${intervalOpts} \
    ${dbsnpOpts} \
    -O ${intervalBed.baseName}_${sampleId}.g.vcf \
    -ERC GVCF
  """
}

gvcfHaplotypeCallerCh = params.noGVCF ? gvcfHaplotypeCallerCh.close() :  gvcfHaplotypeCallerCh.groupTuple(by:[0, 1, 2]).dump(tag:'GVCF HaplotypeCaller')

// STEP GATK HAPLOTYPECALLER.2

process genotypeGVCFs {
  label 'gatk'
  tag "${sampleId}-${sampleName}-${intervalBed.baseName}"

  input:
  tuple val(sampleId),
    val(sampleName),
    file(intervalBed),
    file(gvcf) from gvcfGenotypeGVCFsCh
  file(dbsnp) from dbsnpCh
  file(dbsnpIndex) from dbsnpIndexCh
  file(dict) from dictCh
  file(fasta) from fastaCh
  file(fastaFai) from fastaFaiCh

  output:
  tuple val("HaplotypeCaller"),
    val(sampleId),
    val(sampleName),
    file("${intervalBed.baseName}_${sampleId}.vcf") into vcfGenotypeGVCFsCh

  when: 'haplotypecaller' in tools

  script:
  // Using -L is important for speed and we have to index the interval files also
  intervalOpts = params.noIntervals ? "" : "-L ${intervalBed}"
  dbsnpOpts = params.dbsnp ? "--D ${dbsnp}" : ""
  """
  gatk --java-options -Xmx${task.memory.toGiga()}g \
    IndexFeatureFile -I ${gvcf}

  gatk --java-options -Xmx${task.memory.toGiga()}g \
    GenotypeGVCFs \
    -R ${fasta} \
    ${intervalOpts} \
    ${dbsnpOpts} \
    -V ${gvcf} \
    -O ${intervalBed.baseName}_${sampleId}.vcf
  """
}
vcfGenotypeGVCFsCh = vcfGenotypeGVCFsCh.groupTuple(by:[0, 1, 2])

// STEP GATK MUTECT2.1 - RAW CALLS

ponIndexCh = Channel.value(params.ponIndex ? file(params.ponIndex) : ponIndexBuiltCh)

process mutect2 {
  tag "${sampleIdNormal}_vs_${sampleIdTumor}-${intervalBed.baseName}"
  label 'gatk'
  label 'minCpu'

  input:
  tuple val(sampleIdNormal),
    val(sampleNameNormal),
    val(vCTypeNormal),
    file(bamNormal),
    file(baiNormal),
    val(sampleIdTumor),
    val(sampleNameTumor),
    val(vCTypeTumor),
    file(bamTumor),
    file(baiTumor),
    file(intervalBed) from pairBamMutect2Ch
  file(dict) from dictCh
  file(fasta) from fastaCh
  file(fastaFai) from fastaFaiCh
  file(germlineResource) from germlineResourceCh
  file(germlineResourceIndex) from germlineResourceIndexCh
  file(intervals) from intervalsCh
  file(ponIndex) from ponIndexCh

  output:
  tuple val("Mutect2"),
    val(pairName),
    val("${sampleIdTumor}_vs_${sampleIdNormal}"),
    file("${intervalBed.baseName}_${sampleIdTumor}_vs_${sampleIdNormal}.vcf") into mutect2OutputCh
  tuple val(pairName),
    val(sampleIdTumor),
    val(sampleIdNormal),
    file("${intervalBed.baseName}_${sampleIdTumor}_vs_${sampleIdNormal}.vcf.stats") optional true into mutect2StatsCh

  when: 'mutect2' in tools && vCTypeNormal == 'SNV' && vCTypeTumor == 'SNV'

  script:
  pairName = pairMap[[sampleIdNormal, sampleIdTumor]]
  // please make a panel-of-normals, using at least 40 samples
  // https://gatkforums.broadinstitute.org/gatk/discussion/11136/how-to-call-somatic-mutations-using-gatk4-mutect2
  PON = params.pon ? "--panel-of-normals ${pon}" : ""
  intervalOpts = params.noIntervals ? "" : "-L ${intervalBed}"
  """
  # Get raw calls
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    Mutect2 \
    -R ${fasta}\
    -I ${bamTumor}  -tumor ${sampleIdTumor} \
    -I ${bamNormal} -normal ${sampleIdNormal} \
    ${intervalOpts} \
    --germline-resource ${germlineResource} \
    ${PON} \
    -O ${intervalBed.baseName}_${sampleIdTumor}_vs_${sampleIdNormal}.vcf
  """
}

mutect2OutputCh = mutect2OutputCh.groupTuple(by:[0,1,2])
(mutect2OutputCh, mutect2OutForStatsCh) = mutect2OutputCh.into(2)

mutect2StatsCh = mutect2StatsCh.groupTuple(by:[0,1,2])
(mutect2StatsCh, intervalStatsFilesCh) = mutect2StatsCh.into(2)

// STEP GATK MUTECT2.2 - MERGING STATS

process mergeMutect2Stats {
  label 'gatk'

  tag "${sampleIdTumor}_vs_${sampleIdNormal}"

  publishDir "${params.outDir}/VariantCalling/${sampleIdTumor}_vs_${sampleIdNormal}/Mutect2", mode: params.publishDirMode

  input:
  tuple val(caller),
    val(pairName),
    val(sampleIdTumor_vs_sampleIdNormal),
    file(vcfFiles) from mutect2OutForStatsCh // corresponding small VCF chunks
  tuple val(pairName),
    val(sampleIdTumor),
    val(sampleIdNormal),
    file(statsFiles) from mutect2StatsCh               // the actual stats files
  file(dict) from dictCh
  file(fasta) from fastaCh
  file(fastaFai) from fastaFaiCh
  file(germlineResource) from germlineResourceCh
  file(germlineResourceIndex) from germlineResourceIndexCh
  file(intervals) from intervalsCh

  output:
  file("${sampleIdTumor}_vs_${sampleIdNormal}.vcf.gz.stats") into mergedStatsFileCh

  when: 'mutect2' in tools

  script:
  stats = statsFiles.collect{ "-stats ${it} " }.join(' ')
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    MergeMutectStats \
    ${stats} \
    -O ${sampleIdTumor}_vs_${sampleIdNormal}.vcf.gz.stats
  """
}

// we are merging the VCFs that are called separatelly for different intervals
// so we can have a single sorted VCF containing all the calls for a given caller

// STEP MERGING VCF - GATK HAPLOTYPECALLER & GATK MUTECT2 (UNFILTERED)

vcfConcatenateVCFsCh = mutect2OutputCh.mix(vcfGenotypeGVCFsCh, gvcfHaplotypeCallerCh)
vcfConcatenateVCFsCh = vcfConcatenateVCFsCh.dump(tag:'VCF to merge')

process concatVCF {
  label 'bcftools'
  label 'higCpu'

  tag "${variantCaller}-${sampleId}"

  publishDir "${params.outDir}/VariantCalling/${sampleId}/${variantCaller}", mode: params.publishDirMode

  input:
  tuple val(variantCaller),
    val(sampleId),
    val(sampleIdTN),
    file(vcFiles) from vcfConcatenateVCFsCh
  file(fastaFai) from fastaFaiCh
  file(targetBED) from targetBedCh

  output:
  // we have this funny *_* pattern to avoid copying the raw calls to publishdir
  tuple val(variantCaller),
    val(sampleId),
    val(sampleIdTN),
    file("*_*.vcf.gz"),
    file("*_*.vcf.gz.tbi") into vcfConcatenatedCh
  file("v_bcftools.txt") into bcftoolsVersionCh

  when: ('haplotypecaller' in tools || 'mutect2' in tools)

  script:
  if (variantCaller == 'HaplotypeCallerGVCF')
    outputFile = "HaplotypeCaller_${sampleId}.g.vcf"
  else if (variantCaller == "Mutect2")
    outputFile = "unfiltered_${variantCaller}_${sampleId}.vcf"
  else
    outputFile = "${variantCaller}_${sampleId}.vcf"
  options = params.targetBED ? "-t ${targetBED}" : ""
  intervalsOptions = params.noIntervals ? "-n" : ""
  """
  apConcatenateVCFs.sh -i ${fastaFai} -c ${task.cpus} -o ${outputFile} ${options} ${intervalsOptions}
  bcftools --version &> v_bcftools.txt 2>&1 || true
  """
}

vcfConcatenatedCh
  .branch {
    vcfConcatForFilterCh: it[0] == "Mutect2"
    otherCh: true
  }.set { vcfConcatenatedForks }
(vcfConcatenatedForFilterCh, vcfConcatenatedCh) = [vcfConcatenatedForks.vcfConcatForFilterCh, vcfConcatenatedForks.otherCh]

vcfConcatenatedCh = vcfConcatenatedCh.dump(tag:'VCF')

// STEP GATK MUTECT2.3 - GENERATING PILEUP SUMMARIES

process pileupSummariesForMutect2 {
  label 'gatk'
  label 'minCpu'

  tag "${sampleIdTumor}_vs_${sampleIdNormal}_${intervalBed.baseName}"

  input:
  tuple val(sampleIdNormal),
    val(sampleNameNormal),
    val(vCType),
    file(bamNormal),
    file(baiNormal),
    val(sampleIdTumor),
    val(sampleNameTumor),
    val(vCType),
    file(bamTumor),
    file(baiTumor),
    file(intervalBed) from pairBamPileupSummariesCh
  tuple val(sampleId),
    val(sampleNameTumor),
    val(sampleNameNormal),
    file(statsFile) from intervalStatsFilesCh
  file(germlineResource) from germlineResourceCh
  file(germlineResourceIndex) from germlineResourceIndexCh

  output:
  tuple val(pairName),
    val(sampleIdTumor),
    file("${intervalBed.baseName}_${sampleIdTumor}_pileupsummaries.table") into pileupSummariesCh

  when: 'mutect2' in tools && vCType == 'SNV'

  script:
  pairName = pairMap[[sampleIdNormal, sampleIdTumor]]
  intervalOpts = params.noIntervals ? "-L ${germlineResource}" : "-L ${intervalBed}"
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    GetPileupSummaries \
    -I ${bamTumor} \
    -V ${germlineResource} \
    ${intervalOpts} \
    -O ${intervalBed.baseName}_${sampleIdTumor}_pileupsummaries.table
  """
}

pileupSummariesCh = pileupSummariesCh.groupTuple(by:[0,1])

// STEP GATK MUTECT2.4 - MERGING PILEUP SUMMARIES

process mergePileupSummaries {
  label 'gatk'
  label 'minCpu'

  tag "${pairName}_${sampleIdTumor}"

  publishDir "${params.outDir}/VariantCalling/${sampleIdTumor}/Mutect2", mode: params.publishDirMode

  input:
  tuple val(pairName),
    val(sampleIdTumor),
    file(pileupSums) from pileupSummariesCh
  file(dict) from dictCh

  output:
  file("${sampleIdTumor}_pileupsummaries.table.tsv") into mergedPileupFileCh

  when: 'mutect2' in tools

  script:
  allPileups = pileupSums.collect{ "-I ${it} " }.join(' ')
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    GatherPileupSummaries \
    --sequence-dictionary ${dict} \
    ${allPileups} \
    -O ${sampleIdTumor}_pileupsummaries.table.tsv
  """
}

// STEP GATK MUTECT2.5 - CALCULATING CONTAMINATION

process calculateContamination {
  label 'gatk'
  label 'minCpu'

  tag "${sampleIdTumor}_vs_${sampleIdNormal}"

  publishDir "${params.outDir}/VariantCalling/${sampleIdTumor}/Mutect2", mode: params.publishDirMode

  input:
  tuple val(sampleIdNormal),
    val(sampleNameNormal),
    val(vCType),
    file(bamNormal),
    file(baiNormal),
    val(sampleIdTumor),
    val(sampleNameTumor),
    val(vCType),
    file(bamTumor),
    file(baiTumor) from pairBamCalculateContaminationCh
  file("${sampleIdTumor}_pileupsummaries.table") from mergedPileupFileCh

  output:
  file("${sampleIdTumor}_contamination.table") into contaminationTableCh

  when: 'mutect2' in tools && vCType == 'SNV'

  script:
  """
  # calculate contamination
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    CalculateContamination \
    -I ${sampleIdTumor}_pileupsummaries.table \
    -O ${sampleIdTumor}_contamination.table
  """
}

// STEP GATK MUTECT2.6 - FILTERING CALLS

process filterMutect2Calls {
  label 'gatk'
  label 'medCpu'
  label 'medMem'

  tag "${sampleIdTN}"

  publishDir "${params.outDir}/VariantCalling/${sampleNameTN}/${"$variantCaller"}", mode: params.publishDirMode

  input:
  tuple val(variantCaller),
    val(sampleId),
    val(sampleIdTN),
    file(unfiltered),
    file(unfilteredIndex) from vcfConcatenatedForFilterCh
  file("${sampleIdTN}.vcf.gz.stats") from mergedStatsFileCh
  file("${sampleIdTN}_contamination.table") from contaminationTableCh
  file(dict) from dictCh
  file(fasta) from fastaCh
  file(fastaFai) from fastaFaiCh
  file(germlineResource) from germlineResourceCh
  file(germlineResourceIndex) from germlineResourceIndexCh
  file(intervals) from intervalsCh

  output:
  tuple val("Mutect2"),
    val(sampleId),
    val(sampleIdTN),
    file("filtered_${variantCaller}_${sampleIdTN}.vcf.gz"),
    file("filtered_${variantCaller}_${sampleIdTN}.vcf.gz.tbi"),
    file("filtered_${variantCaller}_${sampleIdTN}.vcf.gz.filteringStats.tsv") into filteredMutect2OutputCh

  when: 'mutect2' in tools

  script:
  """
  # do the actual filtering
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    FilterMutectCalls \
    -V ${unfiltered} \
    --contamination-table ${sampleIdTN}_contamination.table \
    --stats ${sampleIdTN}.vcf.gz.stats \
    -R ${fasta} \
    -O filtered_${variantCaller}_${sampleIdTN}.vcf.gz
  """
}


/*
================================================================================
                            SV VARIANT CALLING
================================================================================
*/

// STEP MANTA.1 - SINGLE MODE

process mantaSingle {
  label 'manta'
  label 'highCpu'
  label 'highMem'

  tag "${sampleId}"

  publishDir "${params.outDir}/VariantCalling/${sampleId}/Manta", mode: params.publishDirMode

  input:
  tuple val(sampleId),
    val(sampleName),
    val(vCType),
    file(bam),
    file(bai) from bamMantaSingleCh
  file(fasta) from fastaCh
  file(fastaFai) from fastaFaiCh
  file(targetBED) from targetBedCh

  output:
  tuple val("Manta"),
    val(sampleId),
    val(sampleName),
    file("*.vcf.gz"),
    file("*.vcf.gz.tbi") into vcfMantaSingleCh
  file('v_manta.txt') into mantaSingleVersionCh

  when: 'manta' in tools && vCType == 'SV'

  script:
  beforeScript = params.targetBED ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
  options = params.targetBED ? "--exome --callRegions call_targets.bed.gz" : ""
  status = statusMap[sampleId]
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
    Manta_${sampleId}.candidateSmallIndels.vcf.gz
  mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
    Manta_${sampleId}.candidateSmallIndels.vcf.gz.tbi
  mv Manta/results/variants/candidateSV.vcf.gz \
    Manta_${sampleId}.candidateSV.vcf.gz
  mv Manta/results/variants/candidateSV.vcf.gz.tbi \
    Manta_${sampleId}.candidateSV.vcf.gz.tbi
  mv Manta/results/variants/${vcftype}SV.vcf.gz \
    Manta_${sampleId}.${vcftype}SV.vcf.gz
  mv Manta/results/variants/${vcftype}SV.vcf.gz.tbi \
    Manta_${sampleId}.${vcftype}SV.vcf.gz.tbi
  configManta.py --version &> v_manta.txt 2>&1 || true
  """
}

vcfMantaSingleCh = vcfMantaSingleCh.dump(tag:'Single Manta')

// STEP MANTA.2 - SOMATIC PAIR

process manta {
  label 'manta'
  label 'highCpu'
  label 'highMem'

  tag "${sampleIdTumor}_vs_${sampleIdNormal}_${vCType}"

  publishDir "${params.outDir}/VariantCalling/${sampleIdTumor}_vs_${sampleIdNormal}/Manta", mode: params.publishDirMode

  input:
  tuple val(sampleIdNormal),
    val(sampleNameNormal),
    val(vCType),
    file(bamNormal),
    file(baiNormal),
    val(sampleIdTumor),
    val(sampleNameTumor),
    val(vCType),
    file(bamTumor),
    file(baiTumor) from pairBamMantaCh
  file(fasta) from fastaCh
  file(fastaFai) from fastaFaiCh
  file(targetBED) from targetBedCh

  output:
  tuple val("Manta"),
    val(pairName),
    val("${sampleIdTumor}_vs_${sampleIdNormal}"),
    file("*.vcf.gz"),
    file("*.vcf.gz.tbi") into vcfMantaCh
  file('v_manta.txt') into mantaVersionCh

  when: 'manta' in tools && vCType == 'SV'

  script:
  pairName = pairMap[[sampleIdNormal, sampleIdTumor]]
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
        Manta_${sampleIdTumor}_vs_${sampleIdNormal}.candidateSmallIndels.vcf.gz
    mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
        Manta_${sampleIdTumor}_vs_${sampleIdNormal}.candidateSmallIndels.vcf.gz.tbi
    mv Manta/results/variants/candidateSV.vcf.gz \
        Manta_${sampleIdTumor}_vs_${sampleIdNormal}.candidateSV.vcf.gz
    mv Manta/results/variants/candidateSV.vcf.gz.tbi \
        Manta_${sampleIdTumor}_vs_${sampleIdNormal}.candidateSV.vcf.gz.tbi
    mv Manta/results/variants/diploidSV.vcf.gz \
        Manta_${sampleIdTumor}_vs_${sampleIdNormal}.diploidSV.vcf.gz
    mv Manta/results/variants/diploidSV.vcf.gz.tbi \
        Manta_${sampleIdTumor}_vs_${sampleIdNormal}.diploidSV.vcf.gz.tbi
    mv Manta/results/variants/somaticSV.vcf.gz \
        Manta_${sampleIdTumor}_vs_${sampleIdNormal}.somaticSV.vcf.gz
    mv Manta/results/variants/somaticSV.vcf.gz.tbi \
        Manta_${sampleIdTumor}_vs_${sampleIdNormal}.somaticSV.vcf.gz.tbi
    configManta.py --version &> v_manta.txt 2>&1 || true
    """
}

vcfMantaCh = vcfMantaCh.dump(tag:'Manta')
(vcfMantaSomaticSVCh, vcfMantaDiploidSVCh) = vcfMantaCh.into(2)


/*
================================================================================
                            CNV VARIANT CALLING
================================================================================
*/

// STEP ASCAT.1 - ALLELECOUNTER

// Run commands and code from Malin Larsson
// Based on Jesper Eisfeldt's code
process alleleCounter {
  label 'canceritAllelecount'
  label 'lowMem'

  tag "${sampleId}"

  input:
  tuple val(sampleId),
    val(sampleName),
    val(vCType),
    file(bam),
    file(bai) from bamAscatCh
  file(acLoci) from acLociCh
  file(dict) from dictCh
  file(fasta) from fastaCh
  file(fastaFai) from fastaFaiCh

  output:
  tuple val(sampleId),
    val(sampleName),
    file("${sampleId}.alleleCount") into alleleCounterOutCh
  file("v_allelecount.txt") into alleleCountsVersionCh

  when: 'ascat' in tools && vCType == 'SNV'

  script:
  """
  alleleCounter \
    -l ${acLoci} \
    -r ${fasta} \
    -b ${bam} \
    -o ${sampleId}.alleleCount;
  alleleCounter --version &> v_allelecount.txt 2>&1 || true
  """
}

alleleCounterOutCh
  .branch {
    normalCh: statusMap[it[0]] == 0
    tumorCh: statusMap[it[0]] == 1
  }.set { alleleCountOutForks }
(alleleCounterOutNormalCh, alleleCounterOutTumorCh) = [alleleCountOutForks.normalCh, alleleCountOutForks.tumorCh]
alleleCounterOutCh = alleleCounterOutNormalCh.combine(alleleCounterOutTumorCh).filter{ pairMap.containsKey([it[0], it[3]]) }

alleleCounterOutCh = alleleCounterOutCh.map {
  sampleIdNormal, sampleNameNormal, alleleCountOutNormal, sampleIdTumor, sampleNameTumor, alleleCountOutTumor ->
    [sampleIdNormal, sampleIdTumor, alleleCountOutNormal, alleleCountOutTumor]
}
// STEP ASCAT.2 - CONVERTALLELECOUNTS

// R script from Malin Larssons bitbucket repo:
// https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
process convertAlleleCounts {
  label 'ascat'
  label 'lowMem'

  tag "${sampleIdTumor}_vs_${sampleIdNormal}"

  publishDir "${params.outDir}/VariantCalling/${sampleIdTumor}_vs_${sampleIdNormal}/ASCAT", mode: params.publishDirMode

  input:
  tuple val(sampleIdNormal),
    val(sampleIdTumor),
    file(alleleCountNormal),
    file(alleleCountTumor) from alleleCounterOutCh

  output:
  tuple val(sampleIdNormal),
    val(sampleIdTumor),
    file("${sampleIdNormal}.BAF"),
    file("${sampleIdNormal}.LogR"),
    file("${sampleIdTumor}.BAF"),
    file("${sampleIdTumor}.LogR") into convertAlleleCountsOutCh
  file("v_ascat.txt") into convertAlleleCountsVersionCh

  when: 'ascat' in tools

  script:
  pairName = pairMap[[sampleIdNormal, sampleIdTumor]]
  gender = genderMap[sampleIdNormal]
  """
  Rscript ${workflow.projectDir}/bin/apConvertAlleleCounts.r ${sampleIdTumor} ${alleleCountTumor} ${sampleIdNormal} ${alleleCountNormal} ${gender}
  R -e "packageVersion('ASCAT')" > v_ascat.txt
  """
}

// STEP ASCAT.3 - ASCAT

// R scripts from Malin Larssons bitbucket repo:
// https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
process ascat {
  label 'ascat'
  label 'lowMem'

  tag "${sampleIdTumor}_vs_${sampleIdNormal}"

  publishDir "${params.outDir}/VariantCalling/${sampleIdTumor}_vs_${sampleIdNormal}/ASCAT", mode: params.publishDirMode

  input:
  tuple val(sampleIdNormal),
    val(sampleIdTumor),
    file(bafNormal),
    file(logrNormal),
    file(bafTumor),
    file(logrTumor) from convertAlleleCountsOutCh
  file(acLociGC) from acLociGCCh

  output:
  tuple val("ASCAT"),
    val(sampleIdNormal),
    val(sampleIdTumor),
    file("${sampleIdTumor}.*.{png,txt}") into ascatOutCh
  file("v_ascat.txt") into ascatVersionCh

  when: 'ascat' in tools

  script:
  gender = genderMap[sampleIdNormal]
  purity_ploidy = (params.ascat_purity && params.ascat_ploidy) ? "--purity ${params.ascat_purity} --ploidy ${params.ascat_ploidy}" : ""
  """
  for f in *BAF *LogR; do sed 's/chr//g' \$f > tmpFile; mv tmpFile \$f;done
  apRunAscat.r \
    --tumorbaf ${bafTumor} \
    --tumorlogr ${logrTumor} \
    --normalbaf ${bafNormal} \
    --normallogr ${logrNormal} \
    --tumorname ${sampleIdTumor} \
    --basedir ${workflow.projectDir} \
    --gcfile ${acLociGC} \
    --gender ${gender} \
    ${purity_ploidy}
  R -e "packageVersion('ASCAT')" > v_ascat.txt
  """
}

ascatOutCh.dump(tag:'ASCAT')


/*
================================================================================
                                   ANNOTATION
================================================================================
*/

// Remapping channels for QC and annotation

vcfAnnotationCh = Channel.empty().mix(
  filteredMutect2OutputCh.map{
    variantcaller, sampleId, sampleName, vcf, tbi, tsv ->
      [variantcaller, sampleId, vcf]
  },
  vcfConcatenatedCh.map{
    variantcaller, sampleId, sampleName, vcf, tbi ->
      [variantcaller, sampleId, vcf]
  },
  vcfMantaSingleCh.map {
    variantcaller, sampleId, sampleName, vcf, tbi ->
      [variantcaller, sampleId, vcf[2]]
  },
  vcfMantaDiploidSVCh.map {
    variantcaller, sampleId, sampleName, vcf, tbi ->
      [variantcaller, sampleId, vcf[2]]
  },
  vcfMantaSomaticSVCh.map {
    variantcaller, sampleId, sampleName, vcf, tbi ->
      [variantcaller, sampleId, vcf[3]]
  })

if (step == 'annotate') {
  vcfToAnnotateCh = Channel.create()
  vcfNoAnnotateCh = Channel.create()

  if (samplePlanPath == []) {
    // By default, annotates all available vcfs that it can find in the VariantCalling directory
    // Excluding vcfs from and g.vcf from HaplotypeCaller
    // Basically it's: results/VariantCalling/*/{HaplotypeCaller,Manta,Mutect2}/*.vcf.gz
    // Without *SmallIndels.vcf.gz from Manta
    // The small snippet `vcf.minus(vcf.fileName)[-2]` catches sampleId
    // This field is used to output final annotated VCFs in the correct directory
    Channel.empty().mix(
      Channel.fromPath("${params.outDir}/VariantCalling/*/HaplotypeCaller/*.vcf.gz")
        .flatten().map{vcf -> ['HaplotypeCaller', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
      Channel.fromPath("${params.outDir}/VariantCalling/*/Manta/*[!candidate]SV.vcf.gz")
        .flatten().map{vcf -> ['Manta', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
      Channel.fromPath("${params.outDir}/VariantCalling/*/Mutect2/*.vcf.gz")
        .flatten().map{vcf -> ['Mutect2', vcf.minus(vcf.fileName)[-2].toString(), vcf]}
    ).choice(vcfToAnnotateCh, vcfNoAnnotateCh) {
      annotateTools == [] || (annotateTools != [] && it[0] in annotateTools) ? 0 : 1
    }
  } else if (annotateTools == []) {
    // Annotate user-submitted VCFs
    // If user-submitted, assume that the sampleId should be assumed automatically
    vcfToAnnotateCh = Channel.fromPath(samplePlanPath)
      .map{vcf -> ['userspecified', vcf.minus(vcf.fileName)[-2].toString(), vcf]}
  } else exit 1, "specify only tools or files to annotate, not both"

  vcfNoAnnotateCh.close()
  vcfAnnotationCh = vcfAnnotationCh.mix(vcfToAnnotateCh)
}
// as now have the list of VCFs to annotate, the first step is to annotate with allele frequencies, if there are any

// STEP SNPEFF

process snpEff {
  label 'snpeff'

  tag "${sampleId} - ${variantCaller} - ${vcf}"

  publishDir params.outDir, mode: params.publishDirMode, saveAs: {
    if (it == "${reducedVCF}_snpEff.ann.vcf") null
    else "Reports/${sampleId}/snpEff/${it}"
  }

  input:
  tuple val(variantCaller),
    val(sampleId),
    file(vcf) from vcfAnnotationCh
  file(dataDir) from snpEffCacheCh
  val(snpeffDb) from snpeffDbCh

  output:
  tuple file("${reducedVCF}_snpEff.txt"),
    file("${reducedVCF}_snpEff.html"),
    file("${reducedVCF}_snpEff.csv") into snpeffReportCh
  tuple val(variantCaller),
    val(sampleId),
    file("${reducedVCF}_snpEff.ann.vcf") into snpeffVCFCh
  file('v_snpeff.txt') into snpeffVersionCh

  when: 'snpeff' in tools

  script:
  reducedVCF = reduceVCF(vcf.fileName)
  cache = params.snpEffCache ? "-dataDir \${PWD}/${dataDir}" : ""
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
  snpEff -version &> v_snpeff.txt 2>&1 || true
  """
}

snpeffReportCh = snpeffReportCh.dump(tag:'snpEff report')

// STEP COMPRESS AND INDEX VCF.1 - SNPEFF

process compressVCFsnpEff {
  label 'tabix'

  tag "${sampleId} - ${vcf}"

  publishDir "${params.outDir}/Annotation/${sampleId}/snpEff", mode: params.publishDirMode

  input:
  tuple val(variantCaller),
    val(sampleId),
    file(vcf) from snpeffVCFCh

  output:
  tuple val(variantCaller),
    val(sampleId),
    file("*.vcf.gz"),
    file("*.vcf.gz.tbi") into compressVCFsnpEffOutCh

  script:
  """
  bgzip < ${vcf} > ${vcf}.gz
  tabix ${vcf}.gz
  """
}

compressVCFsnpEffOutCh = compressVCFsnpEffOutCh.dump(tag:'VCF')


/*
================================================================================
                                     MultiQC
================================================================================
*/

/**
 * Parse software version numbers
 *
 * @output software_versions_mqc.yaml
 */
// TODO: find a way to get multiqc version ?
process getSoftwareVersions {
  label 'python'

  publishDir path:"${params.outDir}/softwareVersions", mode: params.publishDirMode

  input:
  file('v_ascat.txt') from ascatVersionCh.mix(convertAlleleCountsVersionCh).first().ifEmpty('')
  file('v_allelecount.txt') from alleleCountsVersionCh.first().ifEmpty('')
  file('v_bcftools.txt') from bcftoolsVersionCh.first().ifEmpty('')
  file('v_bwa.txt') from bwaVersionCh.ifEmpty('')
  file('v_fastqc.txt') from fastqcVersionCh.first().ifEmpty('')
  file('v_gatk.txt') from gatkVersionCh.first().ifEmpty('')
  file 'v_preseq.txt' from preseqVersionCh.first().ifEmpty([])
  file('v_manta.txt') from mantaVersionCh.mix(mantaSingleVersionCh).first().ifEmpty('')
  //file('v_qualimap.txt') from qualimapVersionCh.first().ifEmpty('')
  file('v_samtools.txt') from samtoolsIndexBamFileVersionCh.mix(samtoolsIndexBamRecalVersionCh).mix(samtoolsMapReadsVersionCh).mix(samtoolsMergeBamMappedVersionCh).mix(samtoolsMergeBamRecalVersionCh).mix(samtoolsBamFilterVersionCh).first().ifEmpty('')
  file('v_snpeff.txt') from snpeffVersionCh.first().ifEmpty('')

  output:
  file('software_versions_mqc.yaml') into yamlSoftwareVersionCh

  when: !('versions' in skipQC)

  script:
  """
  echo "${workflow.manifest.version}" &> v_pipeline.txt 2>&1 || true
  echo "${workflow.nextflow.version}" &> v_nextflow.txt 2>&1 || true

  apScrapeSoftwareVersions.py &> software_versions_mqc.yaml
  """
}

yamlSoftwareVersionCh = yamlSoftwareVersionCh.dump(tag:'SOFTWARE VERSIONS')

process multiQC {
  label 'multiqc'

  publishDir "${params.outDir}/MultiQC", mode: params.publishDirMode

  when: !(params.skipMultiqc)

  input:
  file(splan) from Channel.value(samplePlanPath ? file(samplePlanPath) : "")
  file(metadata) from metadataCh.ifEmpty([])
  file(multiqcConfig) from Channel.value(params.multiqcConfig ? file(params.multiqcConfig) : "")
  file(workflowSummary) from workflowSummaryCh.collectFile(name: "workflow_summary_mqc.yaml")
  file(versions) from yamlSoftwareVersionCh
  file('Mapping/*') from bwaMqcCh.collect().ifEmpty([])
  file('Mapping/*') from bamStatsMqcCh.collect().ifEmpty([])
  file('Mapping/*') from onTargetReportCh.collect().ifEmpty([])
  file('preseq/*') from preseqStatsCh.collect().ifEmpty([])
  //file ('BamQC/*') from bamQCReportCh.collect().ifEmpty([])
  file('coverage/*') from mosdepthOutputCh.collect().ifEmpty([])
  file('coverage/*') from geneCovMqc.collect().ifEmpty([])
  file('BamQC/*') from fragmentSizeCh.collect().ifEmpty([])
  file('BamQC/*') from wgsMetricsOutputCh.collect().ifEmpty([])
  file('FastQC/*') from fastqcReportCh.collect().ifEmpty([])
  file('MarkDuplicates/*') from markDuplicatesReportCh.collect().ifEmpty([])
  //file ('SamToolsStats/*') from samtoolsStatsReportCh.collect().ifEmpty([])
  file('SnpEff/*') from snpeffReportCh.collect().ifEmpty([])
  file('Identito/*') from clustPolymResultsCh.collect().ifEmpty([])

  output:
  file("*multiqc_report.html") into multiQCOutCh
  file("*_data")

  when: !('multiqc' in skipQC)

  script:
  rtitle = customRunName ? "--title \"$customRunName\"" : ''
  rfilename = customRunName ? "--filename " + customRunName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
  metadataOpts = params.metadata ? "--metadata ${metadata}" : ""
  designOpts= params.design ? "-d ${params.design}" : ""
  isPE = params.singleEnd ? "" : "-p"
  modules_list = "-m custom_content -m fastqc -m preseq -m picard -m gatk -m bcftools -m snpeff -m qualimap -m picard -m mosdepth"
  """
  apStats2MultiQC.sh -s ${splan} ${designOpts} ${isPE}
  mqcHeader.py --splan ${splan} --name "VEGAN" --version ${workflow.manifest.version} ${metadataOpts} > multiqc-config-header.yaml
  multiqc . -f ${rtitle} ${rfilename} -c multiqc-config-header.yaml -c $multiqcConfig $modules_list
  """
}

multiQCOutCh.dump(tag:'MultiQC')

/****************
 * Sub-routines *
 ****************/

process checkDesign{
  label 'python'
  label 'minCpu'
  label 'minMem'
  publishDir "${params.summaryDir}/", mode: 'copy'

  when:
  params.design

  input:
  file(design) from designCheckCh
  file(samplePlan) from samplePlanCheckCh

  script:
  optSE = params.singleEnd ? "--singleEnd" : ""
  """
  apCheckDesign.py -d $design -s $samplePlan ${optSE}
  """
}

process outputDocumentation {
  label 'python'
  label 'minCpu'
  label 'minMem'

  publishDir "${params.summaryDir}/", mode: 'copy'

  input:
  file(outputDocs) from outputDocsCh
  file(images) from outputDocsImagesCh

  output:
  file "resultsDescription.html"

  script:
  """
  markdownToHtml.py $outputDocs -o resultsDescription.html
  """
}

workflow.onComplete {

  reportFields = params + [
    workflowName: workflow.manifest.name,
    version : workflow.manifest.version,
    runName: customRunName,
    success: workflow.success,
    dateComplete: workflow.complete,
    duration:workflow.duration,
    exitStatus: workflow.exitStatus,
    errorMessage: (workflow.errorMessage ?: 'None'),
    errorReport: (workflow.errorReport ?: 'None'),
    commandLine: workflow.commandLine,
    projectDir: workflow.projectDir,
    summary: summary + [
      'Date Started': workflow.start,
      'Date Completed': workflow.complete,
      'Pipeline script file path': workflow.scriptFile,
      'Pipeline script hash ID': workflow.scriptId,
      'Pipeline repository Git URL': workflow.repository ?: null,
      'Pipeline repository Git Commit': workflow.commitId ?: null,
      'Pipeline Git branch/tag': workflow.revision ?: null,
      'Docker image': workflow.container ?: null,
      'Nextflow Version': workflow.nextflow.version,
      'Nextflow Build': workflow.nextflow.build,
      'Nextflow Compile Timestamp': workflow.nextflow.timestamp,
    ].findAll{ it.value != null },
  ]

  // Write reports and send email notifications
  makeReports(workflow, params, reportFields, multiQCOutCh)
}
