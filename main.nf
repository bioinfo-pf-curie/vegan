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
params.putAll(lint(params, paramsWithUsage))

tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase()} : []
SVFilters = params.SVFilters ? params.SVFilters.split(',').collect{it.trim().toLowerCase()} : []
SNVFilters = params.SNVFilters ? params.SNVFilters.split(',').collect{it.trim().toLowerCase()} : []
annotateTools = params.annotateTools ? params.annotateTools.split(',').collect{it.trim().toLowerCase()} : []

customRunName = checkRunName(workflow.runName, params.name)
step = getStep(params.samplePlan, params.step)
samplePlanPath = getPath(step, params.samplePlan, params.outDir)
samplePlanCh = getSamplePlan(samplePlanPath, step, params.singleEnd, params.reads, params.readPaths).dump(tag: 'samplePlanCh')
samplePlanCheckCh = params.samplePlan ? Channel.fromPath(samplePlanPath) : Channel.empty()
(designCh, designCheckCh) = params.design ? [getDesign(params.design), Channel.fromPath(params.design)] : [Channel.empty(), Channel.empty()]

if (params.design){
  (genderMap, statusMap, pairMap) = extractInfos(designCh)
} else {
  log.info """\
========================================================================
  ${colors.greenBold}INFO${colors.reset}: No design file detected
  ${tools && (('ascat' in tools) || ('mutect2' in tools))? "  ${colors.redBold}WARNING${colors.reset}: You need a design file in order to run somatic analysis\n" : ''}\
  ${colors.yellowBold}Somatic variant detection (SNV/CNV/SV) will be skipped
  Please set up a design file '--design' to run these steps${colors.reset}
========================================================================"""
  tools.removeAll(["mutect2", "ascat", "manta"])
}

if (tools && ('manta' in tools) && params.singleEnd) {
  log.info """\
========================================================================
${colors.redBold}WARNING${colors.reset}: Manta is not compatible with singleEnd option 
========================================================================"""
}
/*
================================================================================
                               CHECKING REFERENCES
================================================================================
*/

// Initialize each reference with default values in params.genomes, catch the genome defined on the command line first
// if it was defined
params.putAll([
  fasta: params.genome ? params.genomes[params.genome].fasta ?: false : false,
  gtf: params.genome ? params.genomes[params.genome].gtf ?: false : false,
  polyms: params.genome ? params.genomes[params.genome].polyms ?: false : false,
  acLoci: params.genome ? params.genomes[params.genome].acLoci ?: false : false,
  acLociGC: params.genome ? params.genomes[params.genome].acLociGC ?: false : false,
  bwaIndex: params.genome ? params.genomes[params.genome].bwaIndex ?: false : false,
  dbsnp: params.genome ? params.genomes[params.genome].dbsnp ?: false : false,
  dbsnpIndex: params.genome && params.genomes[params.genome].dbsnp ? params.genomes[params.genome].dbsnpIndex ?: false : false,
  dict: params.genome && params.genomes[params.genome].fasta ? params.genomes[params.genome].dict ?: false : false,
  fastaFai: params.genome && params.genomes[params.genome].fasta ? params.genomes[params.genome].fastaFai ?: false : false,
  germlineResource: params.genome ? params.genomes[params.genome].germlineResource ?: false : false,
  germlineResourceIndex: params.genome && params.genomes[params.genome].germlineResource ? params.genomes[params.genome].germlineResourceIndex ?: false : false,
  intervals: params.genome ? params.genomes[params.genome].intervals ?: false : false,
  knownIndels: params.genome ? params.genomes[params.genome].knownIndels ?: false : false,
  knownIndelsIndex: params.genome && params.genomes[params.genome].knownIndels ? params.genomes[params.genome].knownIndelsIndex ?: false : false,
  snpeffDb: params.genome ? params.genomes[params.genome].snpeffDb ?: false : false,
  snpeffCache: params.genome ? params.genomes[params.genome].snpeffCache ?: false : false,
  outDir: params.outDir ? file(params.outDir).toAbsolutePath() : './results',
  bamDir: params.bamDir ? file(params.bamDir).toAbsolutePath() : "${params.outDir}/preprocessing/bams",
  filteredBamDir: params.filteredBamDir ? file(params.filteredBamDir).toAbsolutePath() : "${params.bamDir}/filtering",
  bqsrBamDir: params.bqsrBamDir ? file(params.bqsrBamDir).toAbsolutePath() : "${params.bamDir}/bqsr",
  summaryDir: params.summaryDir ? file(params.summaryDir).toAbsolutePath() : "${params.outDir}/summary",
])

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
  'Genome': params.genome,
  'Fasta': params.fasta ?: null,
  'Target BED': params.targetBED ?: null,
  'Intervals': params.noIntervals ? 'Do not use' : params.intervals,
  'Step': step ?: null,
  'Tools': params.tools ? tools instanceof Collection ? tools.join(', ') : tools: null,
  'QC tools skip': params.skipQC ? 'Yes' : 'No',
  'SV filters': params.SVFilters ? SVFilters instanceof Collection  ? SVFilters.join(', ') : SVFilters : null,
  'SNV filters': params.SNVFilters ? SNVFilters instanceof Collection  ? SNVFilters.join(', ') : SNVFilters : null,
  'Polyms': params.polyms ?: null,
  'BwaIndex': params.bwaIndex ?: null,
  'GermlineResource': params.germlineResource ?: null,
  'acLoci': params.acLoci ?: null,
  'acLociGC': params.acLociGC ?: null,
  'dbsnp': params.dbsnp ?: null,
  'snpeffDb': params.snpeffDb ?: null,
  'snpeffCache': params.snpeffCache ?: null,
  'Save GVCF': 'haplotypecaller' in tools ? params.saveGVCF ? 'Yes' : 'No' : null,
  'Sequenced by': params.sequencingCenter ? params.sequencingCenter: null,
  'Panel of normals': params.pon && 'mutect2' in tools ? params.pon: null,
  'Save Genome Index': params.saveGenomeIndex ? 'Yes' : 'No',
  'Output dir': params.outDir,
  'Launch dir': workflow.launchDir,
  'Working dir': workflow.workDir,
  'Script dir': workflow.projectDir,
  'User': workflow.userName,
  'Config Profile': workflow.profile,
  'Config Description': params.configProfileDescription ?: null,
  'Config Contact': params.configProfileContact ?: null,
  'Config URL': params.configProfileUrl ?: null,
].findAll{ it.value != null }

// Check the hostnames against configured profiles
checkHostname(params, workflow)

/*
================================================================================
                               INIT CHANNELS
================================================================================
*/

// Initialize channels based on params
fastaCh = params.fasta && !('annotate' in step) ? Channel.value(file(params.fasta)) : "null"
fastaFaiCh = params.fastaFai ? Channel.value(file(params.fastaFai)) : fastaFaiBuiltCh
dictCh = params.dict ? Channel.value(file(params.dict)) : dictBuiltCh
polymsCh = params.polyms ? Channel.value(file(params.polyms)) : "null"
gtfCh = params.gtf ? Channel.value(file(params.gtf)) : "null"

// databases
acLociCh = params.acLoci ? Channel.value(file(params.acLoci)) : "null"
acLociGCCh = params.acLociGC ? Channel.value(file(params.acLociGC)) : "null"
dbsnpCh = params.dbsnp ? Channel.value(file(params.dbsnp)) : "null"
dbsnpIndexCh = params.dbsnp ? params.dbsnpIndex ? Channel.value(file(params.dbsnpIndex)) : "null" : "null"
germlineResourceCh = params.germlineResource ? Channel.value(file(params.germlineResource)) : "null"
germlineResourceIndexCh = params.germlineResource ? params.germlineResourceIndex ? Channel.value(file(params.germlineResourceIndex)) : "null" : "null"
knownIndelsCh = params.knownIndels ? Channel.value(file(params.knownIndels)) : "null"
knownIndelsIndexCh = params.knownIndels ? params.knownIndelsIndex ? Channel.value(file(params.knownIndelsIndex)) : "null" : "null"
snpeffCacheCh = params.snpeffCache ? Channel.value(file(params.snpeffCache)) : "null"
snpeffDbCh = params.snpeffDb ? Channel.value(params.snpeffDb) : "null"
polymHeaderCh = Channel.fromPath("$baseDir/assets/polym_header.txt")

// Optional files, not defined within the params.genomes[params.genome] scope
intervalsCh = params.intervals && !params.noIntervals ? Channel.value(file(params.intervals)) : "null"
ponCh = params.pon ? Channel.value(file(params.pon)) : "null"
targetBedCh = params.targetBED ? Channel.value(file(params.targetBED)) : "null"
ponIndexCh = Channel.value(params.ponIndex ? file(params.ponIndex) : "null")

// Print summary and genareta summary channel
workflowSummaryCh = summarize(params, summary, workflow)
metadataCh = params.metadata ? Channel.fromPath(params.metadata) : "null"
outputDocsCh = file("$projectDir/docs/output.md", checkIfExists: true)
outputDocsImagesCh = file("$projectDir/docs/images/", checkIfExists: true)

// Bwa Indexes
if (params.bwaIndex){
  lastPath = params.bwaIndex.lastIndexOf(File.separator)
  bwaDir =  params.bwaIndex.substring(0,lastPath+1)
  bwaBase = params.bwaIndex.substring(lastPath+1)
  Channel
    .fromPath(bwaDir, checkIfExists: true)
    .ifEmpty {exit 1, "BWA index file not found: ${params.bwaIndex}"}
    .combine( [ bwaBase ] )
    .dump(tag :'bwaindexch')
    .set { bwaIndexCh }
} else {
  exit 1, "BWA index file not found: ${params.bwaIndex}"
}

/*
================================================================================
                                  INTERVALS
================================================================================
*/

// STEP 0: CREATING INTERVALS FOR PARALLELIZATION (PREPROCESSING AND VARIANT CALLING)

process createIntervalBeds {
  label 'onlyLinux'

  input:
  file(intervals) from intervalsCh

  output:
  file('*.bed') into bedIntervalsCh mode flatten

  //when: (!params.noIntervals) && step != 'annotate'

  script:
  // If the interval file is BED format, the fifth column is interpreted to
  // contain runtime estimates, which is then used to combine short-running jobs
  if (params.noIntervals)
    """
    echo "noIntervals\n" > noIntervals.bed
    """
  else if (hasExtension(intervals, "bed"))
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

if (!params.noIntervals){
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
}

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
    return it[0,1] + [[it[0], runIds[it[0]].toString()].join("_")] + it[2..-1]
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
// TODO - should be the same for singleEnd
inputBamCh = inputBamCh.map {
    sampleId, sampleName, runID, inputFiles ->
    [sampleId, sampleName, runID, inputFiles[0]]
}
(inputBamCh, inputBamFastQCCh) = inputBamCh.into(2)
(inputPairReadsCh, inputPairReadsFastQC) = inputPairReadsCh.dump(tag: "inputPairReadsCh").into(2)
inputPairReadsCh = inputPairReadsCh.dump(tag:'inputPairReadsCh')

/*
 * FastQC
 */

process fastQC {
  label 'fastqc'
  label 'medCpu'
  label 'lowMem'

  tag "${sampleId}"

  publishDir "${params.outDir}/preprocessing/metrics/fastqc", mode: params.publishDirMode, overwrite: true

  input:
  tuple val(sampleId), val(sampleName), val(runId), file(reads) from inputPairReadsFastQC.mix(inputBamFastQCCh).dump(tag: 'inputFastQC')

  output:
  file("*.{html,zip}") into fastqcReportCh
  file("v_fastqc.txt") into fastqcVersionCh

  when:
  !params.skipQC

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

process bwaMem {
  label 'bwa'
  label 'highCpu'
  label 'extraMem'

  tag "${sampleId}"

  publishDir "${params.outDir}/preprocessing/bams/bwa", mode: params.publishDirMode,
               saveAs: {filename ->  if (params.saveAlignedIntermediates) filename}

  input:
  tuple val(sampleId), val(sampleName), val(runId), file(reads), file(index), val(genomeBase) from inputPairReadsCh.combine(bwaIndexCh)

  output:
  tuple val(sampleId), val(sampleName), val(runId), file("${sampleId}.bam") into bamMappedCh
  tuple val(sampleId), val("${sampleId}_${runId}"), file("${sampleId}.bam") into bamMappedBamQCCh
  file('v_samtools.txt') into samtoolsMapReadsVersionCh
  file("v_bwa.txt") into bwaVersionCh

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
  convertToFastq = hasExtension(reads[0], "bam") ? "gatk --java-options -Xmx${task.memory.toGiga()}g SamToFastq --INPUT=${reads[0]} --FASTQ=/dev/stdout --INTERLEAVE=true --NON_PF=true | \\" : ""
  input = hasExtension(reads[0], "bam") ? "-p /dev/stdin - 2> >(tee ${reads[0]}.bwa.stderr.log >&2)" : "$reads"
  """
  ${convertToFastq}
  bwa mem ${params.bwaOptions} -R \"${readGroup}\" -t ${task.cpus} ${index}/${genomeBase} \
  ${input} | \
  samtools sort --threads ${task.cpus} - > ${sampleId}.bam
  samtools --version &> v_samtools.txt 2>&1 || true
  bwa &> v_bwa.txt 2>&1 || true
  """
}

// Sort BAM whether they are standalone or should be merged
singleBamCh = Channel.create()
multipleBamCh = Channel.create()
bamMappedCh.groupTuple(by:[0, 1])
  .map{it -> [it[0], it[1], it[3]]}
  .branch {
    singleCh: it[2].size() == 1
    multipleCh: it[2].size() > 1
  }.set { bamMappedForks }
(singleBamCh, multipleBamCh) = [bamMappedForks.singleCh, bamMappedForks.multipleCh]
singleBamCh = singleBamCh.dump(tag:'sbams')
multipleBamCh = multipleBamCh.dump(tag:'mbams')

/*
 * MERGING BAM FROM MULTIPLE LANES
 */

process mergeBamMapped {
  label 'samtools'
  label 'highCpu'
  label 'lowMem'

  tag "${sampleId}"

  publishDir "${params.outDir}/preprocessing/bams/bwa/", mode: params.publishDirMode,
               saveAs: {filename ->  if (params.saveAlignedIntermediates) filename}

  input:
  tuple val(sampleId), val(sampleName), val(bams) from multipleBamCh

  output:
  tuple val(sampleId), val(sampleName), file("*_merged.bam") into mergedBamCh
  file('v_samtools.txt') into samtoolsMergeBamMappedVersionCh

  script:
  inputs = bams.collect{"${it}"}.join(' ')
  """
  samtools merge --threads ${task.cpus} ${sampleId}_merged.bam ${inputs}
  samtools --version &> v_samtools.txt 2>&1 || true
  """
}

mergedBamCh = mergedBamCh.mix(singleBamCh).dump(tag:'mergedBamCh')
(mergedBamCh, mergedBamPreseqCh, mergedBamToStatsCh, mergedBamToIndexCh) = mergedBamCh.into(4)


/*
 * INDEX ALIGNED BAM FILE
 */

process indexBamFile {
  label 'samtools'
  label 'minCpu'
  label 'minMem'

  tag "${sampleId}"

  publishDir "${params.outDir}/preprocessing/bams/bwa", mode: params.publishDirMode,
               saveAs: {filename ->  if (params.saveAlignedIntermediates) filename}

  input:
  tuple val(sampleId), val(sampleName), file(bam) from mergedBamToIndexCh

  output:
  tuple val(sampleId), val(sampleName), file(bam), file("*.bai") into indexedBamCh
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
  label 'minMem'

  tag {sampleId}

  publishDir "${params.outDir}/preprocessing/bams/bwa/logs", mode: params.publishDirMode

  input:
  tuple val(sampleId), val(sampleName), file(bam) from mergedBamToStatsCh

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
 * PRESEQ
 */

process preseq {
  label 'preseq'
  label 'lowCpu'
  label 'medMem'

  tag "${sampleID}"

  publishDir "${params.outDir}/preprocessing/metrics/preseq", mode: params.publishDirMode

  when:
  !params.skipPreseq && !params.skipQC

  input:
  tuple val(sampleID), val(sampleName), file(bam) from mergedBamPreseqCh

  output:
  file("*.ccurve.txt") into preseqStatsCh
  file("v_preseq.txt") into preseqVersionCh

  script:
  defectMode = params.preseqDefect ? '-D' : ''
  """
  preseq &> v_preseq.txt
  preseq lc_extrap -v $defectMode -output ${bam.baseName}.ccurve.txt -bam ${bam} -e 500e+06
  """
}

/*
================================================================================
                                  FILTERING
================================================================================
*/

/*
 * DUPLICATES - SAMBAMBA
 */

process markDuplicates {
  label 'sambamba'
  label 'highCpu'
  label 'highMem'

  tag "${sampleId}"

  publishDir "${params.outDir}/preprocessing/bams/markDuplicates", mode: params.publishDirMode,
              saveAs: {filename -> if ( filename.endsWith("md.flagstats")) "stats/$filename"
                                   else if (params.saveAlignedIntermediates) "$filename" }

  input:
  tuple val(sampleId), val(sampleName), file(bam) from mergedBamCh

  output:
  tuple val(sampleId), val(sampleName), file("${sampleId}.md.bam"), file("${sampleId}.md.bam.bai") into duplicateMarkedBamsCh,mdBamPolymCh
  file ("${sampleId}.md.bam.metrics") into markDuplicatesReportCh

  script:
  """
  sambamba markdup --nthreads ${task.cpus} --tmpdir . ${bam} ${sampleId}.md.bam
  sambamba flagstat --nthreads ${task.cpus} ${sampleId}.md.bam > ${sampleId}.md.flagstats
  """
}

/*
 * BAM on Target
 */

bamsToTargetCh = params.targetBED ? duplicateMarkedBamsCh : Channel.empty()

process bamOnTarget {
  label 'bedtools'
  label 'minCpu'
  label 'lowMem'

  tag {sampleId}

  publishDir "${params.outDir}/preprocessing/bams/onTarget", mode: params.publishDirMode,
             saveAs: {filename -> if ( filename.endsWith("_onTarget.flagstats") ) "stats/$filename"
	                          else if (params.saveAlignedIntermediates) "$filename"}

  when:
  params.targetBED

  input:
  tuple val(sampleId), val(sampleName), file(bam), file(bai) from bamsToTargetCh
  file(targetBED) from targetBedCh

  output:
  tuple val(sampleId), val(sampleName), file("*_onTarget.bam"), file("*_onTarget.bam.bai") into onTargetBamsCh
  file("*_onTarget.flagstats") into onTargetReportCh

  script:
  """
  intersectBed -abam ${bam} -b ${targetBED} > ${bam.baseName}_onTarget.bam
  samtools index ${bam.baseName}_onTarget.bam
  samtools flagstat ${bam.baseName}_onTarget.bam > ${bam.baseName}_onTarget.flagstats
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

process bamFiltering {
  label 'samtools'
  label 'medCpu'
  label 'minMem'
  tag "${sampleId}-${vCType}"

  publishDir "${params.filteredBamDir}", mode: params.publishDirMode

  input:
  tuple val(sampleId), val(sampleName), file(bam), file(bai), val(vCType) from procBamsCh

  output:
  tuple val(sampleId), val(sampleName), val(vCType), file("${sampleId}.filtered.${vCType}.bam"), file("${sampleId}.filtered.${vCType}.bam.bai") into filteredBamCh
  file("${sampleId}.filtered.${vCType}.idxstats") into bamFilterReportCh
  file("*.flagstats") into filteringReportCh
  file('v_samtools.txt') into samtoolsBamFilterVersionCh

  script:
  dupParams = (vCType == 'SNV' && 'duplicates' in SNVFilters) | (vCType == 'SV' && 'duplicates' in SVFilters) ? "-F 0x0400" : ""
  mapqParams = (vCType == 'SNV' && 'mapq' in SNVFilters) | (vCType == 'SV' && 'mapq' in SVFilters) && (params.mapQual > 0) ? "-q ${params.mapQual}" : ""
  singleParams = (vCType == 'SNV' && 'singleton' in SNVFilters) | (vCType == 'SV' && 'single' in SVFilters) ? "-F 0x004 -F 0x008 -f 0x001": "-F 0x004"
  uniqParams =  (vCType == 'SNV' && 'multihits' in SNVFilters) | (vCType == 'SV' && 'multi' in SVFilters) ? "-F 0x100 -F 0x800" :  ""
  uniqFilter = (vCType == 'SNV' && 'multihits' in SNVFilters) | (vCType == 'SV' && 'multi' in SVFilters) ? "| grep -v -e \\\"XA:Z:\\\" -e \\\"SA:Z:\\\" | samtools view -b -" : "| samtools view -b -"
  """
  samtools view -h -@ ${task.cpus} ${uniqParams} ${singleParams} ${dupParams} ${mapqParams} ${bam} ${uniqFilter} > ${sampleId}.filtered.${vCType}.bam
  samtools index ${sampleId}.filtered.${vCType}.bam
  samtools flagstat ${sampleId}.filtered.${vCType}.bam > ${sampleId}.filtered.${vCType}.flagstats
  samtools idxstats ${sampleId}.filtered.${vCType}.bam > ${sampleId}.filtered.${vCType}.idxstats
  samtools stats ${sampleId}.filtered.${vCType}.bam > ${sampleId}.filtered.${vCType}.stats
  samtools --version &> v_samtools.txt 2>&1 || true
  """
}

/*
 * Separate SVN/SV/CNV bams after filtering
 */
                                                                                                                                                                                      

(filteredBamCh, filteredBamCSVCh) = filteredBamCh.dump(tag: "filteredBamCh").into(2)

// Save in a csv file in order to restart after filtering step
filteredBamCSVCh.map { sampleId, sampleName, vCType, bam, bai ->
   "${sampleId},${sampleName},${vCType},${params.filteredBamDir}/${bam.getName()},${params.filteredBamDir}/${bai.getName()}\n"
}.collectFile(
   name: 'samplePlan.filtered.csv', sort: true, storeDir: "${params.outDir}/resume/CSV"
)

// Allow to restart here with a sample plan with filteredBams
filteredBamCh = step in 'recalibrate' ? samplePlanCh : filteredBamCh

filteredBamCh
  .branch { sampleId, sampleName, vCType, bam, bai ->
    snvCh: vCType == 'SNV'
      return [sampleId, sampleName, bam, bai]
    svCh: vCType == 'SV'
      return [sampleId, sampleName, bam, bai]
    otherCh: true
      return [sampleId, sampleName, bam, bai]
  }.set { filteredBamForks }

// Remove VCType from all channels
(filteredSNVBamsCh, filteredSVBamsCh, filteredOtherBamsCh) = [filteredBamForks.snvCh, filteredBamForks.svCh, filteredBamForks.otherCh]

// Copy SNV bams for CNV calling
(filteredSNVBamsCh, filteredCNVBamsCh) = filteredSNVBamsCh.into(2)

/*
================================================================================
                                  QUALITY CHECK
================================================================================
*/

// Run the QC on the SNV bam only if available - on the SV otherwise
if ( ('manta' in tools) && !('ascat' in tools || 'haplotypecaller' in tools || 'mutect2' in tools)){
  (filteredSVBamsCh, filteredBamQCCh) = filteredSVBamsCh.into(2)
} else {
  (filteredSNVBamsCh, filteredBamQCCh) = filteredSNVBamsCh.into(2)
}

filteredBamQCCh
  .dump(tag:'qcbams')
  .into {bamInsertSizeCh; bamMosdepthCh; bamGeneCovCh; bamWGSmetricsCh }

/*
 * INSERT SIZE
 */

process getFragmentSize {
  label 'picard'
  label 'minCpu'
  label 'medMem'

  tag "${sampleId}"

  publishDir path: "${params.outDir}/preprocessing/metrics/fragSize", mode: params.publishDirMode

  when:
  !params.singleEnd && !params.skipQC

  input:
  tuple val(sampleId), val(sampleName), file(bam), file(bai) from bamInsertSizeCh

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

  publishDir path: "${params.outDir}/preprocessing/metrics/depth", mode: params.publishDirMode

  when:
  !params.skipQC

  input:
  tuple val(sampleId), val(sampleName), file(bam), file(bai) from bamMosdepthCh
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
  label 'minCpu'
  label 'minMem'

  when:
  !params.skipQC

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
  label 'minCpu'
  label 'lowMem'
  publishDir path: "${params.outDir}/preprocessing/metrics/depth", mode: params.publishDirMode

  tag "${sampleId}"

  when:
  !params.skipQC

  input:
  tuple val(sampleId), val(sampleName), file(bam), file(bai) from bamGeneCovCh
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
  label 'minCpu'
  label 'lowMem'

  publishDir path: "${params.outDir}/preprocessing/metrics/WGSmetrics", mode: params.publishDirMode

  when:
  !params.skipQC

  input:
  tuple val(sampleId), val(sampleName), file(bam), file(bai) from bamWGSmetricsCh
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
       I=${bam} \
       O=${bam.baseName}_collect_wgs_metrics.txt \
       R=${reference} \
       ${bedCmd}
  """
}

/*
 * IDENTITO MONITORING
 */

process getPolym {
  label 'lowCpu'
  label 'medMem'
  label 'polym'

  publishDir "${params.outDir}/preprocessing/metrics/identito", mode: params.publishDirMode

  when:
  !params.skipIdentito && !params.skipQC

  input:
  file(fasta) from fastaCh
  file(fastaFai) from fastaFaiCh
  file(polyms) from polymsCh
  tuple val(sampleId), val(sampleName), file("${sampleId}.md.bam"), file("${sampleId}.md.bam.bai") from mdBamPolymCh

  output:
  file("${sampleId}_matrix.tsv") into clustPolymCh

  script:
  """
  bcftools mpileup -R ${polyms} -f ${fasta} -x -A -B -q 20 -I -Q 0 -d 1000 --annotate FORMAT/DP,FORMAT/AD ${sampleId}.md.bam > ${sampleId}_bcftools.vcf
  SnpSift extractFields -e "."  -s ";" ${sampleId}_bcftools.vcf CHROM POS REF ALT GEN[*].DP GEN[*].AD > ${sampleId}_bcftools.tsv
  apComputePolym.R ${sampleId}_bcftools.tsv ${sampleId}_matrix.tsv ${sampleId} ${polyms}
  """
}

process computePolym {
  label 'lowCpu'
  label 'lowMem'
  label 'polym'

  publishDir "${params.outDir}/preprocessing/metrics/identito", mode: params.publishDirMode

  when:
  !params.skipIdentito && !params.skipQC

  input:
  file(matrix) from clustPolymCh.collect()
  file(header) from polymHeaderCh

  output:
  file("*.csv") into clustPolymResultsCh

  script:
  """
  cat *matrix.tsv |awk 'NR==1{print \$0}' > clust_mat.tsv

  for in_matrix in ${matrix}
  do
  awk 'NR>1{print \$0}' \$in_matrix >> clust_mat.tsv
  done

  apComputeClust.R clust_mat.tsv . clustering_plot
  cat polym_header.txt > temp_clust_mat.csv
  cat clustering_plot_identito.csv >> temp_clust_mat.csv
  mv temp_clust_mat.csv clustering_plot_identito.csv
  """
}

/*
================================================================================
                          RECALIBRATING
================================================================================
Recalibration is performed for SNV bams only
SNVs filtered Bams are used for both SNV and CNV calling
*/

// Channels for baseRecalibration
(bamBaseRecalibratorCh, bamBaseRecalibratorToJoinCh) = params.skipBQSR ? [Channel.empty(), Channel.empty()] : filteredSNVBamsCh.into(2)
bamBaseRecalibratorCh = params.skipBQSR ? bamBaseRecalibratorCh : bamBaseRecalibratorCh.combine(intBaseRecalibratorCh).dump(tag: "bamBaseRecalibratorCh")

/*
 * CREATING RECALIBRATION TABLES
 */

process baseRecalibrator {
  label 'gatk'
  label 'minCpu'
  label 'medMem'

  tag "${sampleId}"

  publishDir "${params.bqsrBamDir}", mode: params.publishDirMode,
             saveAs: {filename ->  if (params.noIntervals) filename}

  input:
  tuple val(sampleId), val(sampleName), file(bam), file(bai), file(intervalBed) from bamBaseRecalibratorCh.dump(tag:'bambqsr')
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
    file("${prefix}.recal.table") into tableGatherBQSRReportsCh

  when: 'haplotypecaller' in tools || 'mutect2' in tools

  script:
  dbsnpOptions = dbsnp.collect{"--known-sites ${it}"}.join(' ')
  knownOptions = knownIndels.collect{"--known-sites ${it}"}.join(' ')
  prefix = params.noIntervals ? "${sampleId}" : "${sampleId}_${intervalBed.baseName}"
  intervalsOptions = params.noIntervals ? "" : "-L ${intervalBed}"
  """
  gatk --java-options -Xmx${task.memory.toGiga()}g \
      BaseRecalibrator \
      -I ${bam} \
      -O ${prefix}.recal.table \
      --tmp-dir ${params.baseRecalibratorTmpDir} \
      -R ${fasta} \
      ${intervalsOptions} \
      ${dbsnpOptions} \
      ${knownOptions} \
      --verbosity INFO
  """
}

/*
 * MERGE BQSR TABLES PER INTERVALS
 */

if (!params.noIntervals) {
  tableGatherBQSRReportsCh = tableGatherBQSRReportsCh.groupTuple(by:[0, 1])
}else{
  (tableGatherBQSRReportsCh, recalTableCh) = tableGatherBQSRReportsCh.into(2)
}

process gatherBQSRReports {
  label 'gatk'
  label 'lowMem'
  label 'lowCpu'

  tag "${sampleId}"

  publishDir "${params.bqsrBamDir}", mode: params.publishDirMode, overwrite: false

  input:
  tuple val(sampleId),
    val(sampleName),
    file(recal) from tableGatherBQSRReportsCh

  output:
  tuple val(sampleId),
    val(sampleName),
    file("${sampleId}.recal.table") into recalTableCh

  when: !(params.noIntervals)

  script:
  input = recal.collect{"-I ${it}"}.join(' ')
  """
  gatk --java-options -Xmx${task.memory.toGiga()}g \
      GatherBQSRReports \
      ${input} \
      -O ${sampleId}.recal.table \
  """
}

bamApplyBQSRCh = bamBaseRecalibratorToJoinCh
  .join(recalTableCh, by:[0, 1])
  .combine(intApplyBQSRCh)
  .dump(tag: 'bamApplyBQSRCh')


/*
 * RECALIBRATING
 */

process applyBQSR {
  label 'gatk'
  label 'lowMem'
  label 'lowCpu'

  tag "${sampleId}-${sampleName}-${intervalBed.baseName}"

  input:
  tuple val(sampleId), val(sampleName), file(bam), file(bai), file(recalibrationReport), file(intervalBed) from bamApplyBQSRCh.dump(tag:'bqsr')
  file(dict) from dictCh
  file(fasta) from fastaCh
  file(fastaFai) from fastaFaiCh

  output:
  tuple val(sampleId), val(sampleName), file("${prefix}${sampleId}.recal.bam") into bamMergeBamRecalCh
  file("v_gatk.txt") into gatkVersionCh

  script:
  prefix = params.noIntervals ? "noInterval_" : "${intervalBed.baseName}_"
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

//bamMergeBamRecalCh = bamMergeBamRecalCh.groupTuple(by:[0, 1])
(bamMergeBamRecalCh, bamMergeBamRecalNoIntCh) = bamMergeBamRecalCh.groupTuple(by:[0, 1]).into(2)


/*
 * MERGING THE RECALIBRATED BAM FILES
 */

process mergeAndIndexBamRecal {
  label 'samtools'
  label 'medCpu'
  label 'medMem'

  tag "${sampleId}-${sampleName}"

  publishDir "${params.bqsrBamDir}", mode: params.publishDirMode

  input:
  tuple val(sampleId), val(sampleName), file(bam) from bamMergeBamRecalCh

  output:
  tuple val(sampleId),
    val(sampleName),
    file("*recal.bam"),
    file("*recal.bam.bai") into bamRecalCh
  file('v_samtools.txt') into samtoolsMergeBamRecalVersionCh

  when: !(params.noIntervals)

  script:
  """
  samtools merge --threads ${task.cpus} ${sampleId}_recal.bam ${bam}
  samtools index ${sampleId}_recal.bam
  samtools --version &> v_samtools.txt 2>&1 || true
  """
}

/*
 * INDEXING THE MERGED RECALIBRATED BAM FILES
 */

process indexBamRecal {
  label 'samtools'
  label 'minCpu'
  label 'minMem'

  tag "${sampleId}-${sampleName}"

  publishDir "${params.bqsrBamDir}", mode: params.publishDirMode

  input:
  tuple val(sampleId), val(sampleName), file(bam) from bamMergeBamRecalNoIntCh

  output:
  tuple val(sampleId),
    val(sampleName),
    file(bam),
    file("*bam.bai") into bamRecalNoIntCh
  file('v_samtools.txt') into samtoolsIndexBamRecalVersionCh

  when: params.noIntervals

  script:
  """
  samtools index ${bam}
  samtools --version &> v_samtools.txt 2>&1 || true
  """
}

// Since bamRecalCh or bamRecalNoIntCh are empty according to noIntervals flag, we merge them
(bamRecalCh, bamRecalCSVCh) = bamRecalCh.mix(bamRecalNoIntCh).into(2)

// When no knownIndels for mapping, Channel bamRecalCh is indexedBamCh
// TODO: seems not suited to the actual layout, have to refactor line below (indexed bam are not filtered)
// bamRecalCh = (params.knownIndels && step == 'mapping') ? bamRecalCh : indexedBamCh.flatMap { it -> [it.plus(2, 'SV'), it.plus(2, 'SNV')]}

// Save here filteredSVBamsCh, filteredCNVBamsCh and bamRecalCh or filteredSNVBamsCh if skipBQSR
// into a single samplePlan file
(filteredSVBamsCh, filteredSVBamsCSVCh) = filteredSVBamsCh.into(2)
filteredSVBamsCSVCh = filteredSVBamsCSVCh.map{ it + ["SV"]}
(filteredCNVBamsCh, filteredCNVBamsCSVCh) = filteredCNVBamsCh.into(2)
filteredCNVBamsCSVCh = filteredCNVBamsCSVCh.map{ it + ["CNV"]}

if (params.skipBQSR) {
  (filteredSNVBamsCh, recalOrSNVBamsCSVCh) = filteredSNVBamsCh.into(2)
  recalOrSNVBamsCSVCh = recalOrSNVBamsCSVCh.map{ it + ["SNV"]}
} else {
  (bamRecalCh, recalOrSNVBamsCSVCh) = bamRecalCh.into(2)
  recalOrSNVBamsCSVCh = recalOrSNVBamsCSVCh.map{ it + ["RECAL"]}
}

bamvCCSVCh = filteredSVBamsCSVCh.mix(filteredCNVBamsCSVCh, recalOrSNVBamsCSVCh)
bamvCCSVCh.map { sampleId, sampleName, bam, bai, vCType ->
   (vCType == 'RECAL') ? "${sampleId},${sampleName},${vCType},${params.bqsrBamDir}/${bam.getName()},${params.bqsrBamDir}/${bai.getName()}\n" : "${sampleId},${sampleName},${vCType},${params.filteredBamDir}/${bam.getName()},${params.filteredBamDir}/${bai.getName()}\n" 
}.collectFile(
   name: 'samplePlan.recal.csv', sort: true, storeDir: "${params.outDir}/resume/CSV"
)

/*
================================================================================
                            DESIGN / PAIRED ANALYSIS
================================================================================
*/

if (step in 'variantcalling') {
  // if step is variantcalling, start from samplePlanCh in order to init 
  // filteredSVBamsCh, filteredCNVBamsCh and filteredSNVBamsCh or bamRecalCh 
  // according to skipBQSR option
  samplePlanCh.branch{ sampleId, sampleName, vCType, bam, bai ->
    snvCh:  vCType == 'SNV'
      return [sampleId, sampleName, bam, bai]
    recalCh: vCType == 'RECAL'
      return [sampleId, sampleName, bam, bai]
    sVCh: vCType == 'SV'
      return [sampleId, sampleName, bam, bai]
    cNVCh: vCType == 'CNV'
      return [sampleId, sampleName, bam, bai]
    otherCh: true
      return [sampleId, sampleName, bam, bai]
  }.set{ samplePlanForks }
  bamRecalCh = params.skipBQSR ? samplePlanForks.snvCh : samplePlanForks.recalCh
  (filteredSVBamsCh, filteredCNVBamsCh) = [samplePlanForks.sVCh, samplePlanForks.cNVCh]
} else {
  bamRecalCh = params.skipBQSR ? filteredSNVBamsCh : bamRecalCh
}

// By default, HaplotypeCaller can be run without design in germline mode
(bamRecalAllCh, bamRecalAllTempCh) = bamRecalCh.into(2)
bamHaplotypeCallerCh = bamRecalAllTempCh.combine(intHaplotypeCallerCh) 

if (params.design) {

  // Manta requires some annotation, even in Single mode
  (filteredSVBamsCh, bamMantaSingleCh) = filteredSVBamsCh.into(2)
                                                                                                                                           
  // CNV Bams
  bamAscatCh = filteredCNVBamsCh

  // separate BAM by status for somatic variant calling
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

  // Crossing Normal and Tumor to get a T/N pair
  // Remapping channel to remove common key sampleId
  pairBamsSNVCh = bamRecalNormalCh.combine(bamRecalTumorCh)
  pairBamsSNVCh = pairBamsSNVCh.filter{ pairMap.containsKey([it[0], it[4]]) }
  pairBamsSVCh = bamSVNormalCh.combine(bamSVTumorCh)
  pairBamsSVCh = pairBamsSVCh.filter{ pairMap.containsKey([it[0], it[4]]) }
  
  // Manta,  Mutect2
  (pairBamCalculateContaminationCh, pairBamsSNVCh) = pairBamsSNVCh.into(2)
  (pairBamMantaCh, pairBamsSVCh) = pairBamsSVCh.into(2)
  intervalPairBamCh = pairBamsSNVCh.combine(bedIntervalsCh)

  // intervals for Mutect2 calls and pileups for Mutect2 filtering
  (pairBamMutect2Ch, pairBamPileupSummariesCh) = intervalPairBamCh.into(2)

} else {
  pairBamPileupSummariesCh = Channel.empty()
  pairBamCalculateContaminationCh = Channel.empty()
  pairBamMutect2Ch = Channel.empty()
  pairBamMantaCh = Channel.empty()
  bamMantaSingleCh = Channel.empty()
  bamAscatCh = Channel.empty()
}

/*
================================================================================
                            SNV VARIANT CALLING
================================================================================
*/

// To speed Variant Callers up we are chopping the reference into smaller pieces
// Do variant calling by this intervals, and re-merge the VCFs

/*
 * HAPLOTYPECALLER
 */

process haplotypeCaller {
  label 'gatk'
  label 'medMemSq'
  label 'lowCpu'

  tag "${sampleId}-${intervalBed.baseName}"

  input:
  tuple val(sampleId), val(sampleName), file(bam), file(bai), file(intervalBed) from bamHaplotypeCallerCh
  file(dbsnp) from dbsnpCh
  file(dbsnpIndex) from dbsnpIndexCh
  file(dict) from dictCh
  file(fasta) from fastaCh
  file(fastaFai) from fastaFaiCh

  output:
  tuple val("HaplotypeCallerGVCF"), val(sampleId), val(sampleName), file("${intervalBed.baseName}_${sampleId}.g.vcf") into gvcfHaplotypeCallerCh
  tuple val(sampleId), val(sampleName), file(intervalBed), file("${intervalBed.baseName}_${sampleId}.g.vcf") into gvcfGenotypeGVCFsCh

  when: 'haplotypecaller' in tools

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

gvcfHaplotypeCallerCh = !params.saveGVCF ? gvcfHaplotypeCallerCh.close() :  gvcfHaplotypeCallerCh.groupTuple(by:[0, 1, 2]).dump(tag:'GVCF')

process genotypeGVCFs {
  label 'gatk'
  tag "${sampleId}-${sampleName}-${intervalBed.baseName}"

  input:
  tuple val(sampleId), val(sampleName), file(intervalBed), file(gvcf) from gvcfGenotypeGVCFsCh
  file(dbsnp) from dbsnpCh
  file(dbsnpIndex) from dbsnpIndexCh
  file(dict) from dictCh
  file(fasta) from fastaCh
  file(fastaFai) from fastaFaiCh

  output:
  tuple val("HaplotypeCaller"), val(sampleId), val(sampleName), file("${intervalBed.baseName}_${sampleId}.vcf") into vcfGenotypeGVCFsCh

  when: 'haplotypecaller' in tools

  script:
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

/*
 * MUTECT2
 */

// STEP GATK MUTECT2.1 - RAW CALLS

process mutect2 {
  tag "${sampleIdNormal}_vs_${sampleIdTumor}-${intervalBed.baseName}"
  label 'gatk'
  label 'minCpu'
  label 'medMem'

  input:
  tuple val(sampleIdNormal), val(sampleNameNormal), file(bamNormal), file(baiNormal),
        val(sampleIdTumor), val(sampleNameTumor), file(bamTumor), file(baiTumor),
        file(intervalBed) from pairBamMutect2Ch
  file(dict) from dictCh
  file(fasta) from fastaCh
  file(fastaFai) from fastaFaiCh
  file(germlineResource) from germlineResourceCh
  file(germlineResourceIndex) from germlineResourceIndexCh
  file(intervals) from intervalsCh
  file(ponIndex) from ponIndexCh

  output:
  tuple val("Mutect2"), val(pairName), val("${sampleIdTumor}_vs_${sampleIdNormal}"), file("${intervalBed.baseName}_${sampleIdTumor}_vs_${sampleIdNormal}.vcf") into mutect2OutputCh
  tuple val(pairName), val(sampleIdTumor), val(sampleIdNormal), file("${intervalBed.baseName}_${sampleIdTumor}_vs_${sampleIdNormal}.vcf.stats") optional true into mutect2StatsCh

  when: 'mutect2' in tools

  script:
  pairName = pairMap[[sampleIdNormal, sampleIdTumor]]
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
mutect2StatsCh = mutect2StatsCh.groupTuple(by:[0,1,2])

// STEP GATK MUTECT2.2 - MERGING STATS

process mergeMutect2Stats {
  label 'gatk'
  label 'minCpu'
  label 'medMem'

  tag "${sampleIdTumor}_vs_${sampleIdNormal}"

  publishDir "${params.outDir}/Mutect2/stats/", mode: params.publishDirMode

  input:
  tuple val(pairName), val(sampleIdTumor), val(sampleIdNormal), file(statsFiles) from mutect2StatsCh
  file(dict) from dictCh
  file(fasta) from fastaCh
  file(fastaFai) from fastaFaiCh
  file(germlineResource) from germlineResourceCh
  file(germlineResourceIndex) from germlineResourceIndexCh
  file(intervals) from intervalsCh

  output:
  tuple val(pairName), val("${sampleIdTumor}_vs_${sampleIdNormal}"), file("${sampleIdTumor}_vs_${sampleIdNormal}.vcf.gz.stats") into mergedStatsFileCh

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

// STEP MERGING VCF - /!\ GATK HAPLOTYPECALLER & GATK MUTECT2 (UNFILTERED) /!\

// we are merging the VCFs that are called separatelly for different intervals
// so we can have a single sorted VCF containing all the calls for a given caller

vcfConcatenateVCFsCh = mutect2OutputCh.mix(vcfGenotypeGVCFsCh, gvcfHaplotypeCallerCh)

process concatVCF {
  label 'bcftools'
  label 'highCpu'
  label 'medMem'

  tag "${variantCaller}-${sampleId}"

  publishDir "${params.outDir}/", mode: params.publishDirMode,
             saveAs: { filename -> if ("${variantCaller}"=="Mutect2") "Mutect2/$filename"
                       else "HaplotypeCaller/$filename" }

  input:
  tuple val(variantCaller), val(sampleId), val(sampleIdTN), file(vcFiles) from vcfConcatenateVCFsCh
  file(fastaFai) from fastaFaiCh
  file(targetBED) from targetBedCh

  output:
  tuple val(variantCaller), val(sampleId), val(sampleIdTN), file("*_*.vcf.gz"), file("*_*.vcf.gz.tbi") into vcfConcatenatedCh
  file("v_bcftools.txt") into bcftoolsVersionCh

  when: ('haplotypecaller' in tools || 'mutect2' in tools)

  script:
  if (variantCaller == 'HaplotypeCallerGVCF')
    outputFile = "${sampleId}_HaplotypeCaller.g.vcf"
  else if (variantCaller == 'Mutect2')
    outputFile = "${sampleIdTN}_Mutect2_unfiltered.vcf"
  else
    outputFile = "${sampleId}_${variantCaller}.vcf"

  options = params.targetBED ? "-t ${targetBED}" : ""
  intervalsOptions = params.noIntervals ? "-n" : ""
  """
  apConcatenateVCFs.sh -i ${fastaFai} -c ${task.cpus} -o ${outputFile} ${options} ${intervalsOptions}
  bcftools --version &> v_bcftools.txt 2>&1 || true
  """
}

// Seperate Mutect2 vs HC and remove GVCF from annotation
vcfConcatenatedCh
  .branch {
    vcfMutect2: it[0] == "Mutect2"
    gvcf: it[0] == "HaplotypeCallerGVCF"
    other: true
  }.set { vcfConcatenatedForks }
(vcfConcatenatedForMutect2FilterCh, vcfConcatenatedHaplotypeCallerGVCFCh, vcfConcatenatedCh) = [vcfConcatenatedForks.vcfMutect2, vcfConcatenatedForks.gvcf, vcfConcatenatedForks.other]
(transitionCh, vcfForMqcStatsCh, vcfForAnnotationCh) = vcfConcatenatedCh.into(3)


/*
 * VCF STATS BEFORE VARIANTS FILTERING
 */

process collectVCFmetrics {
  label 'minCpu'
  label 'minMem'
  label 'onlyLinux'

  input:
  tuple val(variantCaller), val(sampleId), val(sampleName), file(vcf), file(tbi) from vcfForMqcStatsCh

  output:
  file("*.mqc") into callingMetricsMqcCh

  """
  getCallingMetrics.sh -i ${vcf} \
                       -n ${sampleId} > ${sampleId}_${variantCaller}_callingMetrics.mqc
  """
}

/*
 * TRANSITION/TRANSVERSION RATIO
 */

process computeTransition {
  label 'minCpu'
  label 'lowMem'
  label 'transition'

  publishDir "${params.outDir}/${variantCaller}/transition", mode: params.publishDirMode

  input:
  tuple val(variantCaller),
        val(sampleId),
        val(sampleName),
        file(vcf),
        file(tbi) from transitionCh.filter{it[0] == "HaplotypeCaller"}.dump(tag:'transi')

  output:
   file("*table.tsv") into transitionPerSampleCh

  when: 'haplotypecaller' in tools

  script:
  """
  apParseTransition.py -i $vcf -o ${vcf.baseName}.transi.tsv
  apTransition.R ${vcf.baseName}.transi.tsv ${vcf.baseName}.table.tsv
  """
}

// STEP GATK MUTECT2.3 - GENERATING PILEUP SUMMARIES

process pileupSummariesForMutect2 {
  label 'gatk'
  label 'minCpu'
  label 'extraMem'

  tag "${sampleIdTumor}_vs_${sampleIdNormal}_${intervalBed.baseName}"

  input:
  tuple val(sampleIdNormal), val(sampleNameNormal), file(bamNormal), file(baiNormal),
        val(sampleIdTumor), val(sampleNameTumor), file(bamTumor), file(baiTumor),
        file(intervalBed) from pairBamPileupSummariesCh
  file(germlineResource) from germlineResourceCh
  file(germlineResourceIndex) from germlineResourceIndexCh

  output:
  tuple val(pairName), val(sampleIdNormal), val(sampleIdTumor), file("*_pileupsummaries.table") into pileupSummariesCh

  when: 'mutect2' in tools && !params.skipMutectContamination

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

pileupSummariesCh = pileupSummariesCh.groupTuple(by:[0,1,2])

// STEP GATK MUTECT2.4 - MERGING PILEUP SUMMARIES

process mergePileupSummaries {
  label 'gatk'
  label 'minCpu'
  label 'medMem'

  tag "${pairName}_${sampleIdTumor}"

  publishDir "${params.outDir}/Mutect2/contamination", mode: params.publishDirMode

  input:
  tuple val(pairName), val(sampleIdNormal), val(sampleIdTumor), file(pileupSums) from pileupSummariesCh
  file(dict) from dictCh

  output:
  tuple val(sampleIdNormal), val(sampleIdTumor), file("${sampleIdTumor}_pileupsummaries.table.tsv") into mergedPileupFileCh

  when: 'mutect2' in tools  && !params.skipMutectContamination

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

pairBamCalculateContaminationCh = pairBamCalculateContaminationCh.map{
  sampleIdNormal,sampleNameNormal,bamNormal,baiNormal,sampleIdTumor,sampleNameTumor,bamTumor,baiTumor ->
  [sampleIdNormal,sampleIdTumor,sampleNameNormal,sampleNameTumor,bamNormal,baiNormal,bamTumor,baiTumor]
}.join(mergedPileupFileCh, by:[0,1])

process calculateContamination {
  label 'gatk'
  label 'minCpu'
  label 'lowMem'

  tag "${sampleIdTumor}_vs_${sampleIdNormal}"

  publishDir "${params.outDir}/Mutect2/contamination", mode: params.publishDirMode

  input:
  tuple val(sampleIdNormal), val(sampleIdTumor), 
        val(sampleNameNormal), val(sampleNameTumor),
        file(bamNormal), file(baiNormal),
        file(bamTumor), file(baiTumor),
        file(mergedPileup) from pairBamCalculateContaminationCh.dump(tag:'debug5')

  output:
  tuple val(pairName), val("${sampleIdTumor}_vs_${sampleIdNormal}"), file("${sampleIdTumor}_contamination.table.tsv") into contaminationTableCh

  when: 'mutect2' in tools && !params.skipMutectContamination

  script:
  pairName = pairMap[[sampleIdNormal, sampleIdTumor]]
  """
  # calculate contamination
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    CalculateContamination \
    -I ${sampleIdTumor}_pileupsummaries.table.tsv \
    -O ${sampleIdTumor}_contamination.table.tsv
  """
}

// STEP GATK MUTECT2.6 - FILTERING CALLS

if (params.skipMutectContamination){
  contaminationTableCh=file('NO_FILE')
}
mutect2CallsToFilterCh = vcfConcatenatedForMutect2FilterCh.map{
    variantCaller, sampleId, sampleIdTN, vcf, index ->
    [sampleId, sampleIdTN, vcf, index]
}.join(mergedStatsFileCh, by:[0,1]).join(contaminationTableCh, by:[0,1])

process filterMutect2Calls {
  label 'gatk'
  label 'medCpu'
  label 'medMem'

  tag "${sampleIdTN}"

  publishDir "${params.outDir}/Mutect2/", mode: params.publishDirMode,
              saveAs: {filename -> if ( filename.endsWith("callingMetrics.mqc") || filename.endsWith("filteringStats.tsv")) "stats/$filename"
                                   else "$filename"}

  input:
  tuple val(sampleId), val(sampleIdTN),
        file(unfiltered), file(unfilteredIndex), 
        file(statsFile),
        file(contaminationTable) from mutect2CallsToFilterCh
  file(dict) from dictCh
  file(fasta) from fastaCh
  file(fastaFai) from fastaFaiCh
  file(germlineResource) from germlineResourceCh
  file(germlineResourceIndex) from germlineResourceIndexCh
  file(intervals) from intervalsCh

  output:
  tuple val("Mutect2"), val(sampleId), val(sampleIdTN), file("*_filtered.vcf.gz"), file("*_filtered.vcf.gz.tbi"), file("*filteringStats.tsv") into filteredMutect2OutputCh
  file("*.mqc") into mutect2CallingMetricsMqcCh

  when: 'mutect2' in tools

  script:
  contaOpts = contaminationTable.name != 'NO_FILE' ? "--contamination-table ${contaminationTable}" : ""
  contaMetricsOpts = contaminationTable.name != 'NO_FILE' ? "-c ${contaminationTable}" : ""
  """
  # do the actual filtering
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    FilterMutectCalls \
    -V ${unfiltered} \
    ${contaOpts} \
    --stats ${sampleIdTN}.vcf.gz.stats \
    -R ${fasta} \
    -O ${sampleIdTN}_Mutect2_filtered.vcf.gz

  getCallingMetrics.sh -i ${unfiltered} \
                       -f ${sampleIdTN}_Mutect2_filtered.vcf.gz \
                       ${contaMetricsOpts} \
                       -n ${sampleIdTN} > ${sampleIdTN}_Mutect2_callingMetrics.mqc
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

  publishDir "${params.outDir}/Manta", mode: params.publishDirMode

  input:
  tuple val(sampleId), val(sampleName), file(bam), file(bai) from bamMantaSingleCh
  file(fasta) from fastaCh
  file(fastaFai) from fastaFaiCh
  file(targetBED) from targetBedCh

  output:
  tuple val("Manta"), val(sampleId), val(sampleName), file("*.vcf.gz"), file("*.vcf.gz.tbi") into vcfMantaSingleCh
  file('v_manta.txt') into mantaSingleVersionCh

  when: 'manta' in tools && !params.singleEnd

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

  tag "${sampleIdTumor}_vs_${sampleIdNormal}"

  publishDir "${params.outDir}/Manta", mode: params.publishDirMode

  input:
  tuple val(sampleIdNormal), val(sampleNameNormal), file(bamNormal), file(baiNormal),
        val(sampleIdTumor), val(sampleNameTumor), file(bamTumor), file(baiTumor) from pairBamMantaCh
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

  when: 'manta' in tools && !params.singleEnd

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
  label 'ascat'
  label 'lowMem'
  label 'minCpu'

  tag "${sampleId}"

  input:
  tuple val(sampleId), val(sampleName), file(bam), file(bai) from bamAscatCh
  file(acLoci) from acLociCh
  file(dict) from dictCh
  file(fasta) from fastaCh
  file(fastaFai) from fastaFaiCh

  output:
  tuple val(sampleId), val(sampleName), file("${sampleId}.alleleCount") into alleleCounterOutCh
  file("v_allelecount.txt") into alleleCountsVersionCh

  when: 'ascat' in tools

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
  label 'minCpu'

  tag "${sampleIdTumor}_vs_${sampleIdNormal}"

  publishDir "${params.outDir}/ASCAT", mode: params.publishDirMode

  input:
  tuple val(sampleIdNormal), val(sampleIdTumor), file(alleleCountNormal), file(alleleCountTumor) from alleleCounterOutCh

  output:
  tuple val(sampleIdNormal), val(sampleIdTumor), 
        file("${sampleIdNormal}.BAF"), file("${sampleIdNormal}.LogR"),
        file("${sampleIdTumor}.BAF"), file("${sampleIdTumor}.LogR") into convertAlleleCountsOutCh
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
  label 'minCpu'

  tag "${sampleIdTumor}_vs_${sampleIdNormal}"

  publishDir "${params.outDir}/ASCAT", mode: params.publishDirMode

  input:
  tuple val(sampleIdNormal), val(sampleIdTumor),
        file(bafNormal), file(logrNormal),
        file(bafTumor), file(logrTumor) from convertAlleleCountsOutCh
  file(acLociGC) from acLociGCCh

  output:
  tuple val("ASCAT"), val(sampleIdNormal), val(sampleIdTumor), file("${sampleIdTumor}.*.{png,txt}") into ascatOutCh
  file("v_ascat.txt") into ascatVersionCh

  when: 'ascat' in tools

  script:
  gender = genderMap[sampleIdNormal]
  purityPloidy = (params.ascatPurity && params.ascatPloidy) ? "--purity ${params.ascatPurity} --ploidy ${params.ascatPloidy}" : ""
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
    ${purityPloidy}
  R -e "packageVersion('ASCAT')" > v_ascat.txt
  """
}


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
  vcfForAnnotationCh.map{
    variantcaller, sampleId, sampleName, vcf, tbi ->
      [variantcaller, sampleId, vcf]
  }
//  vcfMantaSingleCh.map {
//    variantcaller, sampleId, sampleName, vcf, tbi ->
//      [variantcaller, sampleId, vcf[2]]
//  },
//  vcfMantaDiploidSVCh.map {
//    variantcaller, sampleId, sampleName, vcf, tbi ->
//      [variantcaller, sampleId, vcf[2]]
//  },
//  vcfMantaSomaticSVCh.map {
//    variantcaller, sampleId, sampleName, vcf, tbi ->
//      [variantcaller, sampleId, vcf[3]]
//  }
)

if (step == 'annotate') {
  vcfToAnnotateCh = Channel.create()
  vcfNoAnnotateCh = Channel.create()

  if (samplePlanPath == []) {
    // By default, annotates all available vcfs that it can find in the VariantCalling directory
    // Excluding vcfs from and g.vcf from HaplotypeCaller
    // Basically it's: results/{HaplotypeCaller,Manta,Mutect2}/*.vcf.gz
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


/*
 * SNPEFF
 */

process snpEff {
  label 'snpeff'
  label 'lowMem'
  label 'lowCpu'

  tag "${sampleId} - ${variantCaller} - ${vcf}"

  publishDir "${params.outDir}/${variantCaller}/snpEff/", mode: params.publishDirMode,
             saveAs: { if (it == "${reducedVCF}_snpEff.ann.vcf") null
                       else "reports/${it}" }

  input:
  tuple val(variantCaller), val(sampleId), file(vcf) from vcfAnnotationCh
  file(dataDir) from snpeffCacheCh
  val(snpeffDb) from snpeffDbCh

  output:
  tuple file("${reducedVCF}_snpEff.txt"), file("${reducedVCF}_snpEff.html"), file("${reducedVCF}_snpEff.csv") into snpeffReportCh
  tuple val(variantCaller), val(sampleId), file("${reducedVCF}_snpEff.ann.vcf") into snpeffVCFCh
  file('v_snpeff.txt') into snpeffVersionCh

  when: 'snpeff' in tools

  script:
  reducedVCF = reduceVCF(vcf.fileName)
  cache = params.snpeffCache ? "-dataDir \${PWD}/${dataDir}" : ""
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

// STEP COMPRESS AND INDEX VCF.1 - SNPEFF

process compressVCFsnpEff {
  label 'tabix'
  label 'lowMem'
  label 'lowCpu'

  tag "${sampleId} - ${vcf}"

  publishDir "${params.outDir}/${variantCaller}/snpEff", mode: params.publishDirMode

  input:
  tuple val(variantCaller), val(sampleId), file(vcf) from snpeffVCFCh

  output:
  tuple val(variantCaller), val(sampleId), file("*.vcf.gz"), file("*.vcf.gz.tbi") into compressVCFsnpEffOutCh

  script:
  """
  bgzip < ${vcf} > ${vcf}.gz
  tabix ${vcf}.gz
  """
}

/*
================================================================================
                                     MultiQC
================================================================================
*/

/*
 * Parse software version numbers
 * @output software_versions_mqc.yaml
 */

process getSoftwareVersions {
  label 'python'
  label 'minCpu'
  label 'minMem'

  publishDir path:"${params.outDir}/softwareVersions", mode: params.publishDirMode

  input:
  file('v_ascat.txt') from ascatVersionCh.mix(convertAlleleCountsVersionCh).first().ifEmpty('')
  file('v_allelecount.txt') from alleleCountsVersionCh.first().ifEmpty('')
  file('v_bcftools.txt') from bcftoolsVersionCh.first().ifEmpty('')
  file('v_bwa.txt') from bwaVersionCh.first().ifEmpty('')
  file('v_fastqc.txt') from fastqcVersionCh.first().ifEmpty('')
  file('v_gatk.txt') from gatkVersionCh.first().ifEmpty('')
  file 'v_preseq.txt' from preseqVersionCh.first().ifEmpty([])
  file('v_manta.txt') from mantaVersionCh.mix(mantaSingleVersionCh).first().ifEmpty('')
  file('v_samtools.txt') from samtoolsIndexBamFileVersionCh.mix(samtoolsIndexBamRecalVersionCh).mix(samtoolsMapReadsVersionCh)
                              .mix(samtoolsMergeBamMappedVersionCh).mix(samtoolsMergeBamRecalVersionCh).mix(samtoolsBamFilterVersionCh).first().ifEmpty('')
  file('v_snpeff.txt') from snpeffVersionCh.first().ifEmpty('')

  output:
  file('software_versions_mqc.yaml') into yamlSoftwareVersionCh

  script:
  """
  echo "${workflow.manifest.version}" &> v_pipeline.txt 2>&1 || true
  echo "${workflow.nextflow.version}" &> v_nextflow.txt 2>&1 || true
  apScrapeSoftwareVersions.py &> software_versions_mqc.yaml
  """
}

process multiQC {
  label 'multiqc'
  label 'lowMem'
  label 'minCpu'

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
  file('Mapping/*') from filteringReportCh.collect().ifEmpty([])
  file('preseq/*') from preseqStatsCh.collect().ifEmpty([])
  file('coverage/*') from mosdepthOutputCh.collect().ifEmpty([])
  file('coverage/*') from geneCovMqc.collect().ifEmpty([])
  file('BamQC/*') from fragmentSizeCh.collect().ifEmpty([])
  file('BamQC/*') from wgsMetricsOutputCh.collect().ifEmpty([])
  file('FastQC/*') from fastqcReportCh.collect().ifEmpty([])
  file('MarkDuplicates/*') from markDuplicatesReportCh.collect().ifEmpty([])
  file('SnpEff/*') from snpeffReportCh.collect().ifEmpty([])
  file('Identito/*') from clustPolymResultsCh.collect().ifEmpty([])
  file('vcfMetrics/*') from mutect2CallingMetricsMqcCh.collect().ifEmpty([])
  file('vcfMetrics/*') from callingMetricsMqcCh.collect().ifEmpty([])
  file('Transition/*') from transitionPerSampleCh.collect().ifEmpty([])

  output:
  file("*multiqc_report.html") into multiQCOutCh
  file("*_data")

  script:
  rtitle = customRunName ? "--title \"$customRunName\"" : ''
  rfilename = customRunName ? "--filename " + customRunName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
  metadataOpts = params.metadata ? "--metadata ${metadata}" : ""
  designOpts= params.design ? "-d ${params.design}" : ""
  isPE = params.singleEnd ? "" : "-p"
  modules_list = "-m custom_content -m fastqc -m preseq -m picard -m gatk -m bcftools -m snpeff -m picard -m mosdepth"
  """
  apStats2MultiQC.sh -s ${splan} ${designOpts} ${isPE}
  medianReadNb="\$(sort -t, -k3,3n mqc.stats | awk -F, '{a[i++]=\$3;} END{x=int((i+1)/2); if (x<(i+1)/2) printf "%.0f", (a[x-1]+a[x])/2; else printf "%.0f",a[x-1];}')"
  mqcHeader.py --splan ${splan} --name "VEGAN" --version ${workflow.manifest.version} ${metadataOpts} --nbreads \${medianReadNb} > multiqc-config-header.yaml
  multiqc . -f ${rtitle} ${rfilename} -c multiqc-config-header.yaml -c $multiqcConfig $modules_list
  """
}



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
