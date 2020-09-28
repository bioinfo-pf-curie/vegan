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
                        project : EUCANCAN/nf-vegan
================================================================================
Started February 2020.
--------------------------------------------------------------------------------
nf-vegan: Variant calling pipeline for whole Exome and whole Genome sequencing cANcer data Pipeline.
  An open-source analysis pipeline to detect germline or somatic variants
  from whole genome or targeted sequencing
--------------------------------------------------------------------------------
 @Homepage
 https://gitlab.curie.fr/data-analysis/nf-vegan
--------------------------------------------------------------------------------
 @Documentation
 https://gitlab.curie.fr/data-analysis/nf-vegan/README.md
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
//def params = lintedParams
def paramsWithUsage = readParamsFromJsonSettings("$baseDir/parameters.settings.json")
def params = lint(params, paramsWithUsage)

tools = params.tools
skipQC = params.skipQC
skipFilterSV = params.skipFilterSV
skipFilterSNV = params.skipFilterSNV
annotateTools = params.annotateTools

customRunName = checkRunName(workflow.runName, params.runName)
step    = getStep(params.samplePlan, params.step)
inputPath = getPath(step, params.samplePlan, params.outputDir)
samplePlanCh = getSamplePlan(inputPath)

(genderMap, statusMap, pairMap) = extractInfos(getDesign(params.design))
// genderMap[sampleId] = "XX"/"XY"
// Statusmap[sampleId] = O/1
// pairMap[normalSampleId, tumorSampleId] = "pairName"

/*
================================================================================
                               CHECKING REFERENCES
================================================================================
*/

// Initialize each reference with params.genomes, catch the command line first if it was defined
params << [
    fasta: params.genome && !('annotate' in step) ? params.genomes[params.genome].fasta ?: null : null,
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
def summary = [
        'Pipeline Release': workflow.revision ?: null,
        'Run Name': customRunName,
        'Max Resources': "${params.maxMemory} memory, ${params.maxCpus} cpus, ${params.maxTime} time per job",
        'Container': workflow.containerEngine && workflow.container ? "${workflow.containerEngine} - ${workflow.container}" : null,
        'Input': params.input ?: null,
        'Target BED': params.targetBED ?: null,
        'Step': step ?: null,
        'Tools': params.tools ? params.tools.join(', '): null,
        'QC tools skip': params.skipQC ? skipQC.join(', ') : null,
        'SV tools skip': params.skipFilterSV ? skipFilterSV.join(', ') : null,
        'SNV tools skip': params.skipFilterSNV ? skipFilterSNV.join(', ') : null,
        'Intervals': params.noIntervals && step != 'annotate' ? 'Do not use' : null,
        'GVCF': 'haplotypecaller' in tools ? params.noGVCF ? 'No' : 'Yes' : null,
        'Sequenced by': params.sequencingCenter ? params.sequencingCenter: null,
        'Panel of normals': params.pon && 'mutect2' in tools ? params.pon: null,
        'Save Genome Index': params.saveGenomeIndex ? 'Yes' : 'No',
        'Output dir': params.outputDir,
        'Launch dir': workflow.launchDir,
        'Working dir': workflow.workDir,
        'Script dir': workflow.projectDir,
        'User': workflow.userName,
        'Genome': params.genome,
        'Fasta': params.fasta ?: null,
        'FastaFai': params.fastaFai ?: null,
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
germlineResourceCh = params.germlineResource && 'mutect2' in tools ? Channel.value(file(params.germlineResource)) : "null"
intervalsCh = params.intervals && !params.noIntervals && !('annotate' in step) ? Channel.value(file(params.intervals)) : "null"

// knownIndels is currently a list of file for smallGRCh37, so transform it in a channel
knownIndelsList = params.knownIndels && ('mapping' in step) ? params.knownIndels.collect{file(it)} : []
knownIndelsCh = params.knownIndels && params.genome == 'smallGRCh37' ? Channel.value(knownIndelsList.collect()) : params.knownIndels ? Channel.value(file(params.knownIndels)) : "null"

snpEffCacheCh = params.snpEffCache ? Channel.value(file(params.snpEffCache)) : "null"
snpeffDbCh = params.snpeffDb ? Channel.value(params.snpeffDb) : "null"

// Optional files, not defined within the params.genomes[params.genome] scope
ponCh = params.pon ? Channel.value(file(params.pon)) : "null"
targetBEDCh = params.targetBED ? Channel.value(file(params.targetBED)) : "null"

// Print summary and genareta summary channel
workflowSummaryCh = summarize(params, summary, workflow)


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

process BuildDict {
    label 'gatk'
    tag {fasta}

    publishDir params.outputDir, mode: params.publishDirMode,
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

process BuildFastaFai {
    label 'samtools'
    tag {fasta}

    publishDir params.outputDir, mode: params.publishDirMode,
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

process BuildDbsnpIndex {
    label 'tabix'
    tag {dbsnp}

    publishDir params.outputDir, mode: params.publishDirMode,
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

process BuildGermlineResourceIndex {
    label 'tabix'
    tag {germlineResource}

    publishDir params.outputDir, mode: params.publishDirMode,
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

process BuildKnownIndelsIndex {
    label 'tabix'
    tag {knownIndels}

    publishDir params.outputDir, mode: params.publishDirMode,
        saveAs: {params.saveGenomeIndex ? "reference_genome/${it}" : null }

    input:
        each file(knownIndels) from knownIndelsCh

    output:
        file("${knownIndels}.tbi") into knownIndelsIndexBuiltCh

    when: !(params.knownIndelsIndex) && params.knownIndels && 'mapping' in step

    script:
    """
    tabix -p vcf ${knownIndels}
    """
}

knownIndelsIndexCh = params.knownIndels ? params.knownIndelsIndex ? Channel.value(file(params.knownIndelsIndex)) : knownIndelsIndexBuiltCh.collect() : "null"

process BuildPonIndex {
    label 'tabix'
    tag {pon}

    publishDir params.outputDir, mode: params.publishDirMode,
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

process BuildIntervals {
  label 'onlyLinux'
  tag {fastaFai}

  publishDir params.outputDir, mode: params.publishDirMode,
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

process CreateIntervalBeds {
    label 'onlyLinux'
    tag {intervals.fileName}

    input:
        file(intervals) from intervalsCh

    output:
        file '*.bed' into bedIntervalsCh mode flatten

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

bedIntervalsCh = bedIntervalsCh.dump(tag:'bedintervals')

if (params.noIntervals && step != 'annotate') bedIntervalsCh = Channel.from(file("noIntervals.bed"))

(intBaseRecalibratorCh, intApplyBQSRCh, intHaplotypeCallerCh, bedIntervalsCh) = bedIntervalsCh.into(4)

// PREPARING CHANNELS FOR PREPROCESSING AND QC

// TODO: not working
//(inputBamCh, inputPairReadsCh) = step == 'mapping' ? forkMappingSamplePlan(samplePlanCh) : [Channel.empty(), Channel.empty()]
if (step == "mapping") {
    def runIds = [:]
    samplePlanCh.map {
        runIds[it[0]] = runIds.containsKey(it[0]) ? runIds[it[0]] + 1 : 0
        return it[0,1] + [[it[0], runIds[it[0]].toString()].join("_")] + it[2..-1]
    }.branch {
        bamCh: it[3] =~ /.*bam$/
        pairCh: it[3] =~ /.*(fastq.gz|fq.gz|fastq|fq)$/
    }.set { samplePlanForks }
    (inputBamCh, inputPairReadsCh) = [samplePlanForks.bamCh, samplePlanForks.pairCh]
} else (inputBamCh, inputPairReadsCh) = [Channel.empty(), Channel.empty()]


// TODO: chek if splitFastq works
if (params.splitFastq){
    inputPairReadsCh = inputPairReadsCh
    // newly splitfastq are named based on split, so the name is easier to catch
            .splitFastq(by: params.splitFastq, compress:true, file:"split", pe:true)
            .map {sampleId, sampleName, runId, reads1, reads2 ->
                // The split fastq read1 is the 4th element (indexed 3) its name is split_3
                // The split fastq read2's name is split_4
                // It's followed by which split it's acutally based on the mother fastq file
                // Index start at 1
                // Extracting the index to get a new IdRun
                splitIndex = reads1.fileName.toString().minus("split_3.").minus(".gz")
                newIdRun = runId + "_" + splitIndex
                // Giving the files a new nice name
                newReads1 = file("${sampleName}_${newIdRun}_R1.fastq.gz")
                newReads2 = file("${sampleName}_${newIdRun}_R2.fastq.gz")
                [sampleId, sampleName, newIdRun, reads1, reads2]}
}

inputPairReadsCh = inputPairReadsCh.dump(tag:'INPUT')


/*
================================================================================
                                  QUALITY CHECK
================================================================================
*/

(inputBamCh, inputBamFastQCCh) = inputBamCh.into(2)

// Removing inputFile2 wich is null in case of uBAM
inputBamFastQCCh = inputBamFastQCCh.map {
    sampleId, sampleName, runID, inputFile1, inputFile2 ->
    [sampleId, sampleName, runID, inputFile1]
}

(inputPairReadsCh, inputPairReadsFastQC) = inputPairReadsCh.into(2)

// STEP 0.5: QC ON READS

// TODO: Use only one FastQC process for FASTQ and uBAM files ?
// FASTQ and uBAM files are renamed based on the sample name

process FastQCFQ {
    label 'fastqc'
    label 'cpus2'

    tag {sampleId + "-" + runId}

    publishDir "${params.outputDir}/Reports/${sampleName}/FastQC/${sampleName}_${runId}", mode: params.publishDirMode

    input:
        set sampleId, sampleName, runId, file("${sampleName}_${runId}_R1.fastq.gz"), file("${sampleName}_${runId}_R2.fastq.gz") from inputPairReadsFastQC

    output:
        file("*.{html,zip}") into fastQCFQReportCh
        file("v_fastqc.txt") into fastqcVersionCh

    when: !('fastqc' in skipQC)

    script:
    """
    fastqc -t 2 -q ${sampleName}_${runId}_R1.fastq.gz ${sampleName}_${runId}_R2.fastq.gz
    fastqc --version > v_fastqc.txt
    """
}

process FastQCBAM {
    label 'fastqc'
    label 'cpus2'

    tag {sampleId + "-" + runId}

    publishDir "${params.outputDir}/Reports/${sampleName}/FastQC/${sampleName}_${runId}", mode: params.publishDirMode

    input:
        set sampleId, sampleName, runId, file("${sampleName}_${runId}.bam") from inputBamFastQCCh

    output:
        file("*.{html,zip}") into fastQCBAMReportCh
        file("v_fastqc.txt") into fastqcBamVersionCh

    when: !('fastqc' in skipQC)

    script:
    """
    fastqc -t 2 -q ${sampleName}_${runId}.bam
    fastqc --version > v_fastqc.txt
    """
}

fastQCReportCh = fastQCFQReportCh.mix(fastQCBAMReportCh)

fastQCReportCh = fastQCReportCh.dump(tag:'FastQC')


/*
================================================================================
                                  MAPPING
================================================================================
*/

// STEP 1: MAPPING READS TO REFERENCE GENOME WITH BWA MEM

inputPairReadsCh = inputPairReadsCh.dump(tag:'INPUT')

inputPairReadsCh = inputPairReadsCh.mix(inputBamCh)

process MapReads {
    label 'gatk_bwa_samtools'
    label 'cpusMax'
    label 'memoryMax'

    tag {sampleId + "-" + runId}

    input:
//        set idPatient, sampleName, runId, file (inputFile1), file(inputFile2) from inputPairReadsCh
        set sampleId, sampleName, runId, file(inputFile1), file(inputFile2) from inputPairReadsCh
        file(bwaIndex) from bwaIndexCh
        file(fasta) from fastaCh
        file(fastaFai) from fastaFaiCh

    output:
        set sampleId, sampleName, runId, file("${sampleName}_${runId}.bam") into bamMappedCh
        set sampleId, val("${sampleName}_${runId}"), file("${sampleName}_${runId}.bam") into bamMappedBamQCCh
        file 'v_samtools.txt' into samtoolsMapReadsVersionCh

    script:
    // -K is an hidden option, used to fix the number of reads processed by bwa mem
    // Chunk size can affect bwa results, if not specified,
    // the number of threads can change which can give not deterministic result.
    // cf https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md
    // and https://github.com/gatk-workflows/gatk4-data-processing/blob/8ffa26ff4580df4ac3a5aa9e272a4ff6bab44ba2/processing-for-variant-discovery-gatk4.b37.wgs.inputs.json#L29
    CN = params.sequencingCenter ? "CN:${params.sequencingCenter}\\t" : ""
    readGroup = "@RG\\tID:${runId}\\t${CN}PU:${runId}\\tSM:${sampleId}\\tLB:${sampleId}\\tPL:illumina"
    // adjust mismatch penalty for tumor samples
    status = statusMap[sampleId]
    extra = status == 1 ? "-B 3" : ""
    convertToFastq = hasExtension(inputFile1, "bam") ? "gatk --java-options -Xmx${task.memory.toGiga()}g SamToFastq --INPUT=${inputFile1} --FASTQ=/dev/stdout --INTERLEAVE=true --NON_PF=true | \\" : ""
    input = hasExtension(inputFile1, "bam") ? "-p /dev/stdin - 2> >(tee ${inputFile1}.bwa.stderr.log >&2)" : "${inputFile1} ${inputFile2}"
    """
        ${convertToFastq}
        bwa mem ${params.bwaOptions} -R \"${readGroup}\" ${extra} -t ${task.cpus} -M ${fasta} \
        ${input} | \
        samtools sort --threads ${task.cpus} -m 2G - > ${sampleName}_${runId}.bam
        samtools --version &> v_samtools.txt 2>&1 || true
    """
}

bamMappedCh = bamMappedCh.dump(tag:'Mapped BAM')
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
singleBamCh = singleBamCh.dump(tag:'Single BAM')

// STEP 1': MERGING BAM FROM MULTIPLE LANES

process MergeBamMapped {
    label 'samtools'
    label 'cpus8'

    tag {sampleId + "-" + sampleName}

    input:
        set sampleId, sampleName, runId, file(bam) from multipleBamCh

    output:
        set sampleId, sampleName, file("${sampleName}.bam") into mergedBamCh, mergedBamUCh
        file 'v_samtools.txt' into samtoolsMergeBamMappedVersionCh

    script:
    """
    samtools merge --threads ${task.cpus} ${sampleName}.bam ${bam}
    samtools --version &> v_samtools.txt 2>&1 || true
    """
}

mergedBamCh = mergedBamCh.dump(tag:'Merged BAM')

mergedBamCh = mergedBamCh.mix(singleBamCh)

mergedBamCh = mergedBamCh.dump(tag:'BAMs for MD')

(mergedBamCh, mergedBamToIndexCh, mergedBamUCh) = mergedBamCh.into(3)

process IndexBamFile {
    label 'samtools'
    label 'cpus8'

    tag {sampleId + "-" + sampleName}

    input:
        set sampleId, sampleName, file(bam) from mergedBamToIndexCh

    output:
        set sampleId, sampleName, file(bam), file("*.bai") into indexedBamCh
        file 'v_samtools.txt' into samtoolsIndexBamFileVersionCh

    when: !params.knownIndels

    script:
    """
    samtools index ${bam}
    mv ${bam}.bai ${bam.baseName}.bai
    samtools --version &> v_samtools.txt 2>&1 || true
    """
}


/*
================================================================================
                                  FILTERING
================================================================================
*/

// STEP 2: BWAMEM UNIQ FILTER
// Mapping Quality Filter
process BwaMemUniq {
    label 'samtools'
    label 'cpus2'

    tag {sampleId + "-" + sampleName}

    publishDir "${params.outputDir}/Reports/${sampleName}/Uniq", mode: params.publishDirMode

    input:
        set sampleId, sampleName, file(bam) from mergedBamCh

    output:
        set sampleId, sampleName, file("${sampleName}.bam") into memUbamCh
        file("${sampleName}.mapping.stats") into mapUReport
        file 'v_samtools.txt' into samtoolsBwaMemUniqVersionCh

    when: !('uniq' in skipFilterSNV)

    script:

    """
    #removed unmapped also with -F 4
    samtools view  -@ ${task.cpus} -h ${params.samtoolsUniqOptions} ${bam} | grep -v \"XA:Z\" | samtools view  -@ ${task.cpus} -bS > ${sampleName}.temp.bam 2> ${sampleName}.temp.txt 
    samtools sort -@ ${task.cpus} -o ${sampleName}.bam ${sampleName}.temp.bam
    samtools index ${sampleName}.bam
    samtools index ${bam}

    UniqueHits=\$(samtools idxstats ${sampleName}.bam |  awk '{ UNIQ_HIT+=\$3 } END { print UNIQ_HIT }')
    samtools idxstats ${bam} |  awk -v Unique_hits="\$UniqueHits" '{
    Total_reads+=\$3+\$4; Mapped_reads+=\$3; Unmapped+=\$4 } END {
          printf("Total_reads\\t%d\\nMapped_reads\\t%d\\nUnique_hits\\t%d\\nMulti_hits\\t%d\\nUnmapped\\t%d\\n.uniq(%%)\\t%.2f \\n", \
          Total_reads, Mapped_reads, Unique_hits, (Mapped_reads - Unique_hits), Unmapped, (Unique_hits*100/Total_reads))
    }' > ${sampleName}.mapping.stats 
    # clean
    rm ./${sampleName}.temp.* ./*.bam.bai 
    samtools --version &> v_samtools.txt 2>&1 || true
    """
}

if ('uniq' in skipFilterSNV) {
 	memUbamCh = mergedBamUCh
}

// STEP 2: MARKING DUPLICATES FILTER
process MarkDuplicates {
    label 'sambamba'
    label 'cpus16'
    label 'memoryMax'

    tag {sampleId + "-" + sampleName}

    publishDir params.outputDir, mode: params.publishDirMode,
        saveAs: {
            if (it == "${sampleName}.bam.metrics" && (('markduplicates' in skipFilterSNV) || ('markduplicates' in skipFilterSV))) null
            else if (it == "${sampleName}.bam.metrics") "Reports/${sampleName}/MarkDuplicates/${it}"
            else "Preprocessing/${sampleName}/DuplicateMarked/${it}"
        }

    input:
        set sampleId, sampleName, file("${sampleName}.bam") from memUbamCh

    output:
        set sampleId, sampleName, file("${sampleName}.md.bam"), file("${sampleName}.md.bam.bai") into duplicateMarkedBamsCh, duplicateMarkedBamsMQCh
        file ("${sampleName}.bam.metrics") into markDuplicatesReportCh

    when: (params.knownIndels && (!('markduplicates' in skipFilterSNV) || !('markduplicates' in skipFilterSV)))

    script:
    """

    sambamba markdup --remove-duplicates --nthreads ${task.cpus} --tmpdir . ${sampleName}.bam ${sampleName}.md.bam 
    sambamba flagstat --nthreads ${task.cpus} ${sampleName}.md.bam > ${sampleName}.bam.metrics

    """
}

if (('markduplicates' in skipFilterSNV) || ('markduplicates' in skipFilterSV)) markDuplicatesReportCh.close()

// STEP 2: MAPQ FILTER
// Mapping Quality Filter
// TODO: Do we have to use mapQReportCh for multiqc ?
process MapQ {
    label 'samtools'
    label 'cpus2'

    tag {sampleId + "-" + sampleName}

    publishDir "${params.outputDir}/Reports/${sampleName}/MapQ", mode: params.publishDirMode
   // publishDir "${params.outputDir}/Reports/${sampleName}/MapQ", pattern: '*.{bam,bam.bai}', mode: 'copy', overwrite: true


    input:
        set sampleId, sampleName, file(bam), file(bai) from duplicateMarkedBamsMQCh

    output:
        set sampleId, sampleName, file("${sampleName}.recal.bam"), file("${sampleName}.recal.bam.bai") into mapQbamCh
        file("${bam.baseName}.${params.mapQual}.mapping.stats") into mapQReportCh
        file 'v_samtools.txt' into samtoolsMapQVersionCh

    when: !('mapq' in skipFilterSNV)

    script:

    """
    samtools view -@ ${task.cpus} -q ${params.mapQual} -b ${bam} > ${sampleName}.recal.bam
    samtools index ${sampleName}.recal.bam 
    samtools idxstats ${sampleName}.recal.bam |  awk -v id_sample="${sampleName}" -v map_qual="${params.mapQual}" '{
    mapped+=\$3; unmapped+=\$4 } END {
          printf("SAMPLE\\t%s\\nNB\\t%d\\nNB_MAPPED\\t%d\\n.q%d(%%)\\t%.2f \\n", id_sample, mapped+unmapped, mapped, map_qual, (mapped*100/(mapped+unmapped)))
    }' > ${bam.baseName}.${params.mapQual}.mapping.stats
    samtools --version &> v_samtools.txt 2>&1 || true
    """
}

duplicateMarkedBamsCh = duplicateMarkedBamsCh.dump(tag:'MD BAM')
markDuplicatesReportCh = markDuplicatesReportCh.dump(tag:'MD Report')

if ('mapq' in skipFilterSNV) {
 	mapQbamCh = duplicateMarkedBamsCh
}

(bamMDCh, bamMDToJoinCh) = mapQbamCh.into(2) // duplicateMarked + MapQ


/*
================================================================================
                                  RECALIBRATING
================================================================================
*/

bamBaseRecalibratorCh = bamMDCh.combine(intBaseRecalibratorCh)
bamBaseRecalibratorCh = bamBaseRecalibratorCh.dump(tag:'BAM FOR BASERECALIBRATOR')

// STEP 3: CREATING RECALIBRATION TABLES
process BaseRecalibrator {
    label 'gatk'
    label 'cpus1'

    tag {sampleId + "-" + sampleName + "-" + intervalBed.baseName}

    input:
        set sampleId, sampleName, file(bam), file(bai), file(intervalBed) from bamBaseRecalibratorCh
        file(dbsnp) from dbsnpCh
        file(dbsnpIndex) from dbsnpIndexCh
        file(fasta) from fastaCh
        file(dict) from dictCh
        file(fastaFai) from fastaFaiCh
        file(knownIndels) from knownIndelsCh
        file(knownIndelsIndex) from knownIndelsIndexCh

    output:
        set sampleId, sampleName, file("${prefix}${sampleName}.recal.table") into tableGatherBQSRReportsCh
        set sampleId, sampleName into recalTableTSVnoIntCh

    when: params.knownIndels

    script:
    dbsnpOptions = params.dbsnp ? "--known-sites ${dbsnp}" : ""
    knownOptions = params.knownIndels ? knownIndels.collect{"--known-sites ${it}"}.join(' ') : ""
    prefix = params.noIntervals ? "" : "${intervalBed.baseName}_"
    intervalsOptions = params.noIntervals ? "" : "-L ${intervalBed}"
    // TODO: --use-original-qualities ???
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        BaseRecalibrator \
        -I ${bam} \
        -O ${prefix}${sampleName}.recal.table \
        --tmp-dir ${params.baseRecalibratorOpts} \
        -R ${fasta} \
        ${intervalsOptions} \
        ${dbsnpOptions} \
        ${knownOptions} \
        --verbosity INFO
    """
}

if (!params.noIntervals) tableGatherBQSRReportsCh = tableGatherBQSRReportsCh.groupTuple(by:[0, 1])

tableGatherBQSRReportsCh = tableGatherBQSRReportsCh.dump(tag:'BQSR REPORTS')

if (params.noIntervals) {
    (tableGatherBQSRReportsCh, tableGatherBQSRReportsNoIntCh) = tableGatherBQSRReportsCh.into(2)
    recalTableCh = tableGatherBQSRReportsNoIntCh
} else recalTableTSVnoIntCh.close()

// STEP 3.5: MERGING RECALIBRATION TABLES
process GatherBQSRReports {
    label 'gatk'
    label 'memorySingleCPU2Task'
    label 'cpus2'

    tag {sampleId + "-" + sampleName}

    publishDir "${params.outputDir}/Preprocessing/${sampleName}/DuplicateMarked", mode: params.publishDirMode, overwrite: false

    input:
        set sampleId, sampleName, file(recal) from tableGatherBQSRReportsCh

    output:
        set sampleId, sampleName, file("${sampleName}.recal.table") into recalTableCh
        set sampleId, sampleName into recalTableTSVCh

    when: !(params.noIntervals)

    script:
    input = recal.collect{"-I ${it}"}.join(' ')
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        GatherBQSRReports \
        ${input} \
        -O ${sampleName}.recal.table \
    """
}

recalTableCh = recalTableCh.dump(tag:'RECAL TABLE')

(recalTableTSVCh, recalTableSampleTSVCh) = recalTableTSVCh.mix(recalTableTSVnoIntCh).into(2)

// Create TSV files to restart from this step
recalTableTSVCh.map { sampleId, sampleName ->
    status = statusMap[sampleId]
    gender = genderMap[sampleId]
    bam = "${params.outputDir}/Preprocessing/${sampleName}/DuplicateMarked/${sampleName}.md.bam"
    bai = "${params.outputDir}/Preprocessing/${sampleName}/DuplicateMarked/${sampleName}.md.bai"
    recalTable = "${params.outputDir}/Preprocessing/${sampleName}/DuplicateMarked/${sampleName}.recal.table"
    "${sampleId}\t${gender}\t${status}\t${sampleName}\t${bam}\t${bai}\t${recalTable}\n"
}.collectFile(
    name: 'duplicateMarked.tsv', sort: true, storeDir: "${params.outputDir}/Preprocessing/TSV"
)

recalTableSampleTSVCh
    .collectFile(storeDir: "${params.outputDir}/Preprocessing/TSV/") {
        sampleId, sampleName ->
        status = statusMap[sampleId]
        gender = genderMap[sampleId]
        bam = "${params.outputDir}/Preprocessing/${sampleName}/DuplicateMarked/${sampleName}.md.bam"
        bai = "${params.outputDir}/Preprocessing/${sampleName}/DuplicateMarked/${sampleName}.md.bai"
        recalTable = "${params.outputDir}/Preprocessing/${sampleName}/DuplicateMarked/${sampleName}.recal.table"
        ["duplicateMarked_${sampleName}.tsv", "${sampleId}\t${gender}\t${status}\t${sampleName}\t${bam}\t${bai}\t${recalTable}\n"]
}

bamApplyBQSRCh = step in 'recalibrate' ? samplePlanCh : bamMDToJoinCh.join(recalTableCh, by:[0,1])

bamApplyBQSRCh = bamApplyBQSRCh.dump(tag:'BAM + BAI + RECAL TABLE')
// [DUMP: recal.table] ['normal', 'normal', normal.md.bam, normal.md.bai, normal.recal.table]

bamApplyBQSRCh = bamApplyBQSRCh.combine(intApplyBQSRCh)

bamApplyBQSRCh = bamApplyBQSRCh.dump(tag:'BAM + BAI + RECAL TABLE + INT')
// [DUMP: BAM + BAI + RECAL TABLE + INT] ['normal', 'normal', normal.md.bam, normal.md.bai, normal.recal.table, 1_1-200000.bed]

// STEP 4: RECALIBRATING
process ApplyBQSR {
    label 'gatk'
    label 'memorySingleCPU2Task'
    label 'cpus2'

    tag {sampleId + "-" + sampleName + "-" + intervalBed.baseName}

    input:
        set sampleId, sampleName, file(bam), file(bai), file(recalibrationReport), file(intervalBed) from bamApplyBQSRCh
        file(dict) from dictCh
        file(fasta) from fastaCh
        file(fastaFai) from fastaFaiCh

    output:
        set sampleId, sampleName, file("${prefix}${sampleName}.recal.bam") into bamMergeBamRecalCh
        file("v_gatk.txt") into gatkVersionCh

    script:
    prefix = params.noIntervals ? "" : "${intervalBed.baseName}_"
    intervalsOptions = params.noIntervals ? "" : "-L ${intervalBed}"
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        ApplyBQSR \
        -R ${fasta} \
        --input ${bam} \
        --output ${prefix}${sampleName}.recal.bam \
        ${intervalsOptions} \
        --bqsr-recal-file ${recalibrationReport}
    gatk ApplyBQSR --help &> v_gatk.txt 2>&1 || true
    """
}

bamMergeBamRecalCh = bamMergeBamRecalCh.groupTuple(by:[0, 1])
(bamMergeBamRecalCh, bamMergeBamRecalNoIntCh) = bamMergeBamRecalCh.into(2)

// EP 4.5: MERGING THE RECALIBRATED BAM FILES
process MergeBamRecal {
    label 'samtools'
    label 'cpus8'

    tag {sampleId + "-" + sampleName}

    publishDir "${params.outputDir}/Preprocessing/${sampleName}/Recalibrated", mode: params.publishDirMode

    input:
        set sampleId, sampleName, file(bam) from bamMergeBamRecalCh

    output:
        set sampleId, sampleName, file("${sampleName}.recal.bam"), file("${sampleName}.recal.bam.bai") into bamRecalCh
        set sampleId, sampleName, file("${sampleName}.recal.bam") into bamRecalQCCh
        set sampleId, sampleName into bamRecalTSVCh
        file 'v_samtools.txt' into samtoolsMergeBamRecalVersionCh

    when: !(params.noIntervals)

    script:
    """
    samtools merge --threads ${task.cpus} ${sampleName}.recal.bam ${bam}
    samtools index ${sampleName}.recal.bam
    samtools --version &> v_samtools.txt 2>&1 || true
    """
}

// STEP 4.5': INDEXING THE RECALIBRATED BAM FILES
process IndexBamRecal {
    label 'samtools'
    label 'cpus8'

    tag {sampleId + "-" + sampleName}

    publishDir "${params.outputDir}/Preprocessing/${sampleName}/Recalibrated", mode: params.publishDirMode

    input:
        set sampleId, sampleName, file("${sampleName}.recal.bam") from bamMergeBamRecalNoIntCh

    output:
        set sampleId, sampleName, file("${sampleName}.recal.bam"), file("${sampleName}.recal.bam.bai") into bamRecalNoIntCh
        set sampleId, sampleName, file("${sampleName}.recal.bam") into bamRecalQCnoIntCh
        set sampleId, sampleName into bamRecalTSVnoIntCh
        file 'v_samtools.txt' into samtoolsIndexBamRecalVersionCh


    when: params.noIntervals

    script:
    """
    samtools index ${sampleName}.recal.bam
    samtools --version &> v_samtools.txt 2>&1 || true
    """
}

bamRecalCh = bamRecalCh.mix(bamRecalNoIntCh)
bamRecalQCCh = bamRecalQCCh.mix(bamRecalQCnoIntCh)
bamRecalTSVCh = bamRecalTSVCh.mix(bamRecalTSVnoIntCh)

(bamRecalBamQCCh, bamRecalSamToolsStatsCh) = bamRecalQCCh.into(2)
(bamRecalTSVCh, bamRecalSampleTSVCh) = bamRecalTSVCh.into(2)

// Creating a TSV file to restart from this step
bamRecalTSVCh.map { sampleId, sampleName ->
    gender = genderMap[sampleId]
    status = statusMap[sampleId]
    bam = "${params.outputDir}/Preprocessing/${sampleName}/Recalibrated/${sampleName}.recal.bam"
    bai = "${params.outputDir}/Preprocessing/${sampleName}/Recalibrated/${sampleName}.recal.bam.bai"
    "${sampleId}\t${gender}\t${status}\t${sampleName}\t${bam}\t${bai}\n"
}.collectFile(
    name: 'recalibrated.tsv', sort: true, storeDir: "${params.outputDir}/Preprocessing/TSV"
)

bamRecalSampleTSVCh
    .collectFile(storeDir: "${params.outputDir}/Preprocessing/TSV") {
        sampleId, sampleName ->
        status = statusMap[sampleId]
        gender = genderMap[sampleId]
        bam = "${params.outputDir}/Preprocessing/${sampleName}/Recalibrated/${sampleName}.recal.bam"
        bai = "${params.outputDir}/Preprocessing/${sampleName}/Recalibrated/${sampleName}.recal.bam.bai"
        ["recalibrated_${sampleName}.tsv", "${sampleId}\t${gender}\t${status}\t${sampleName}\t${bam}\t${bai}\n"]
}

// When no knownIndels for mapping, Channel bamRecalCh is indexedBamCh
bamRecalCh = (params.knownIndels && step == 'mapping') ? bamRecalCh : indexedBamCh


/*
================================================================================
                                  QUALITY CHECK
================================================================================
*/

process SamtoolsStats {
    label 'samtools'
    label 'cpus2'

    tag {sampleId + "-" + sampleName}

    publishDir "${params.outputDir}/Reports/${sampleName}/SamToolsStats", mode: params.publishDirMode

    input:
        set sampleId, sampleName, file(bam) from bamRecalSamToolsStatsCh

    output:
        file ("${bam}.samtools.stats.out") into samtoolsStatsReportCh
        file 'v_samtools.txt' into samtoolsStatsVersionCh

    when: !('samtoolsStats' in skipQC)

    script:
    """
    samtools stats ${bam} > ${bam}.samtools.stats.out
    samtools --version &> v_samtools.txt 2>&1 || true
    """
}

samtoolsStatsReportCh = samtoolsStatsReportCh.dump(tag:'SAMTools')

bamBamQCCh = bamMappedBamQCCh.mix(bamRecalBamQCCh) // Mapreads + MapQ + MarkDuplicates + ApplyBQSR

process BamQC {
    label 'qualimap'
    label 'memoryMax'
    label 'cpus16'

    tag {sampleId + "-" + sampleName}

    publishDir "${params.outputDir}/Reports/${sampleName}/bamQC", mode: params.publishDirMode

    input:
        set sampleId, sampleName, file(bam) from bamBamQCCh
        file(targetBED) from targetBEDCh

    output:
        file("${bam.baseName}") into bamQCReportCh
        file 'v_qualimap.txt' into qualimapVersionCh

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
    qualimap --version &> v_qualimap.txt 2>&1 || true
    """
}

bamQCReportCh = bamQCReportCh.dump(tag:'BamQC')


/*
================================================================================
                            VARIANT CALLING
================================================================================
*/


// When starting with variant calling, Channel bamRecalCh is inputSampleCh
bamRecalCh = step in 'variantcalling' ? inputSampleCh : bamRecalCh
bamRecalCh = bamRecalCh.dump(tag:'BAM')

// Here we have a recalibrated bam set
// The TSV file is formatted like: "sampleId status sampleName bamFile baiFile"
// Manta will be run in Germline mode, or in Tumor mode depending on status
// HaplotypeCaller will be run for Normal and Tumor samples

(bamMantaSingleCh, bamAscatCh, bamRecalAllCh, bamRecalAllTempCh) = bamRecalCh.into(4)
//(bamAscatCh, bamRecalAllCh) = bamRecalAllCh.into(2)

// separate BAM by status for somatic variant calling
bamRecalAllCh.branch{
    normalCh: statusMap[it[0]] == 0
    tumorCh: statusMap[it[0]] == 1
}.set { bamRecalAllForks }
(bamRecalNormalCh, bamRecalTumorCh) = [bamRecalAllForks.normalCh, bamRecalAllForks.tumorCh]
// Crossing Normal and Tumor to get a T/N pair for Somatic Variant Calling
// Remapping channel to remove common key sampleId
pairBamCh = bamRecalNormalCh.combine(bamRecalTumorCh).filter{ pairMap.containsKey([it[0], it[4]]) }

pairBamCh = pairBamCh.dump(tag:'BAM Somatic Pair')

// Manta,  Mutect2
(pairBamMantaCh, pairBamCalculateContaminationCh, pairBamCh) = pairBamCh.into(3)

intervalPairBamCh = pairBamCh.combine(bedIntervalsCh)

// intervals for Mutect2 calls and pileups for Mutect2 filtering
(pairBamMutect2Ch, pairBamPileupSummariesCh) = intervalPairBamCh.into(2)


/*
================================================================================
                            SNV VARIANT CALLING
================================================================================
*/

// To speed Variant Callers up we are chopping the reference into smaller pieces
// Do variant calling by this intervals, and re-merge the VCFs

bamHaplotypeCallerCh = bamRecalAllTempCh.combine(intHaplotypeCallerCh)

// STEP GATK HAPLOTYPECALLER.1

process HaplotypeCaller {
    label 'gatk'
    label 'memorySingleCPUTaskSq'
    label 'cpus2'

    tag {sampleName + "-" + intervalBed.baseName}

    input:
        set sampleId, sampleName, file(bam), file(bai), file(intervalBed) from bamHaplotypeCallerCh
        file(dbsnp) from dbsnpCh
        file(dbsnpIndex) from dbsnpIndexCh
        file(dict) from dictCh
        file(fasta) from fastaCh
        file(fastaFai) from fastaFaiCh

    output:
        set val("HaplotypeCallerGVCF"), sampleId, sampleName, file("${intervalBed.baseName}_${sampleName}.g.vcf") into gvcfHaplotypeCallerCh
        set sampleId, sampleName, file(intervalBed), file("${intervalBed.baseName}_${sampleName}.g.vcf") into gvcfGenotypeGVCFsCh

    when: 'haplotypecaller' in tools

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g -Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
        HaplotypeCaller \
        -R ${fasta} \
        -I ${bam} \
        -L ${intervalBed} \
        -D ${dbsnp} \
        -O ${intervalBed.baseName}_${sampleName}.g.vcf \
        -ERC GVCF
    """
}

gvcfHaplotypeCallerCh = params.noGVCF ? gvcfHaplotypeCallerCh.close() :  gvcfHaplotypeCallerCh.groupTuple(by:[0, 1, 2]).dump(tag:'GVCF HaplotypeCaller')

// STEP GATK HAPLOTYPECALLER.2

process GenotypeGVCFs {
    label 'gatk'
    tag {sampleName + "-" + intervalBed.baseName}

    input:
        set sampleId, sampleName, file(intervalBed), file(gvcf) from gvcfGenotypeGVCFsCh
        file(dbsnp) from dbsnpCh
        file(dbsnpIndex) from dbsnpIndexCh
        file(dict) from dictCh
        file(fasta) from fastaCh
        file(fastaFai) from fastaFaiCh

    output:
    set val("HaplotypeCaller"), sampleId, sampleName, file("${intervalBed.baseName}_${sampleName}.vcf") into vcfGenotypeGVCFsCh

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
        -O ${intervalBed.baseName}_${sampleName}.vcf
    """
}
vcfGenotypeGVCFsCh = vcfGenotypeGVCFsCh.groupTuple(by:[0, 1, 2])

// STEP GATK MUTECT2.1 - RAW CALLS

process Mutect2 {
    tag {sampleNameTumor + "_vs_" + sampleNameNormal + "-" + intervalBed.baseName}
    label 'gatk'
    label 'cpus_1'

    input:
    set sampleIdNormal, sampleNameNormal, file(bamNormal), file(baiNormal), sampleIdTumor, sampleNameTumor, file(bamTumor), file(baiTumor), file(intervalBed) from pairBamMutect2Ch
    file(dict) from dictCh
    file(fasta) from fastaCh
    file(fastaFai) from fastaFaiCh
    file(germlineResource) from germlineResourceCh
    file(germlineResourceIndex) from germlineResourceIndexCh
    file(intervals) from intervalsCh
    file(ponIndex) from Channel.value(params.ponIndex ? file(params.ponIndex) : ponIndexBuiltCh)


    output:
    set val("Mutect2"),
            pairName,
            val("${sampleNameTumor}_vs_${sampleNameNormal}"),
            file("${intervalBed.baseName}_${sampleNameTumor}_vs_${sampleNameNormal}.vcf") into mutect2OutputCh
    set pairName,
            sampleNameTumor,
            sampleNameNormal,
            file("${intervalBed.baseName}_${sampleNameTumor}_vs_${sampleNameNormal}.vcf.stats") optional true into mutect2StatsCh

    when: 'mutect2' in tools

    script:
    pairName = pairMap[[sampleIdNormal, sampleIdTumor]]
    // please make a panel-of-normals, using at least 40 samples
    // https://gatkforums.broadinstitute.org/gatk/discussion/11136/how-to-call-somatic-mutations-using-gatk4-mutect2
    PON = params.pon ? "--panel-of-normals ${pon}" : ""
    """
    # Get raw calls
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
      Mutect2 \
      -R ${fasta}\
      -I ${bamTumor}  -tumor ${sampleIdTumor} \
      -I ${bamNormal} -normal ${sampleIdNormal} \
      -L ${intervalBed} \
      --germline-resource ${germlineResource} \
      ${PON} \
      -O ${intervalBed.baseName}_${sampleNameTumor}_vs_${sampleNameNormal}.vcf
    """
}

mutect2OutputCh = mutect2OutputCh.groupTuple(by:[0,1,2])
(mutect2OutputCh, mutect2OutForStats) = mutect2OutputCh.into(2)

(mutect2StatsCh, intervalStatsFilesCh) = mutect2StatsCh.into(2)
mutect2StatsCh = mutect2StatsCh.groupTuple(by:[0,1,2])

// STEP GATK MUTECT2.2 - MERGING STATS

process MergeMutect2Stats {
    tag {sampleNameTumor + "_vs_" + sampleNameNormal}
    label 'gatk'

    publishDir "${params.outputDir}/VariantCalling/${sampleNameTumor}_vs_${sampleNameNormal}/Mutect2", mode: params.publishDirMode

    input:
    set caller, pairName, sampleNameTumor_vs_sampleNameNormal, file(vcfFiles) from mutect2OutForStats // corresponding small VCF chunks
    set pairName, sampleNameTumor, sampleNameNormal, file(statsFiles) from mutect2StatsCh               // the actual stats files
    file(dict) from dictCh
    file(fasta) from fastaCh
    file(fastaFai) from fastaFaiCh
    file(germlineResource) from germlineResourceCh
    file(germlineResourceIndex) from germlineResourceIndexCh
    file(intervals) from intervalsCh

    output:
    file("${sampleNameTumor_vs_sampleNameNormal}.vcf.gz.stats") into mergedStatsFileCh

    when: 'mutect2' in tools

    script:
    stats = statsFiles.collect{ "-stats ${it} " }.join(' ')
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        MergeMutectStats \
        ${stats} \
        -O ${sampleNameTumor}_vs_${sampleNameNormal}.vcf.gz.stats
    """
}

// we are merging the VCFs that are called separatelly for different intervals
// so we can have a single sorted VCF containing all the calls for a given caller

// STEP MERGING VCF - GATK HAPLOTYPECALLER & GATK MUTECT2 (UNFILTERED)

vcfConcatenateVCFsCh = mutect2OutputCh.mix(vcfGenotypeGVCFsCh, gvcfHaplotypeCallerCh)
vcfConcatenateVCFsCh = vcfConcatenateVCFsCh.dump(tag:'VCF to merge')

process ConcatVCF {
    label 'bcftools'
    label 'cpus8'

    tag {variantCaller + "-" + sampleName}

    publishDir "${params.outputDir}/VariantCalling/${sampleName}/${"$variantCaller"}", mode: params.publishDirMode

    input:
    set variantCaller, sampleId, sampleName, file(vcFiles) from vcfConcatenateVCFsCh
    file(fastaFai) from fastaFaiCh
    file(targetBED) from targetBEDCh

    output:
    // we have this funny *_* pattern to avoid copying the raw calls to publishdir
    set variantCaller, sampleId, sampleName, file("*_*.vcf.gz"), file("*_*.vcf.gz.tbi") into vcfConcatenatedCh
    file("v_bcftools.txt") into bcftoolsVersionCh

    when: ('haplotypecaller' in tools || 'mutect2' in tools)

    script:
    if (variantCaller == 'HaplotypeCallerGVCF')
        outputFile = "HaplotypeCaller_${sampleName}.g.vcf"
    else if (variantCaller == "Mutect2")
        outputFile = "unfiltered_${variantCaller}_${sampleName}.vcf"
    else
        outputFile = "${variantCaller}_${sampleName}.vcf"
    options = params.targetBED ? "-t ${targetBED}" : ""
    """
    apConcatenateVCFs.sh -i ${fastaFai} -c ${task.cpus} -o ${outputFile} ${options}
    bcftools --version &> v_bcftools.txt 2>&1 || true
    """
}

(vcfConcatenatedCh, vcfConcatenatedForFilterCh) = vcfConcatenatedCh.into(2)
vcfConcatenatedCh = vcfConcatenatedCh.dump(tag:'VCF')

// STEP GATK MUTECT2.3 - GENERATING PILEUP SUMMARIES

process PileupSummariesForMutect2 {
    tag {sampleNameTumor + "_vs_" + sampleNameNormal + "_" + intervalBed.baseName }
    label 'gatk'
    label 'cpus_1'

    input:
    set sampleIdNormal, sampleNameNormal, file(bamNormal), file(baiNormal), sampleIdTumor, sampleNameTumor, file(bamTumor), file(baiTumor), file(intervalBed) from pairBamPileupSummariesCh
    set sampleId, sampleNameNormal, sampleNameTumor, file(statsFile) from intervalStatsFilesCh
    file(germlineResource) from germlineResourceCh
    file(germlineResourceIndex) from germlineResourceIndexCh

    output:
    set pairName,
            sampleNameTumor,
            file("${intervalBed.baseName}_${sampleNameTumor}_pileupsummaries.table") into pileupSummariesCh

    when: 'mutect2' in tools

    script:
    pairName = pairMap[[sampleIdNormal, sampleIdTumor]]
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        GetPileupSummaries \
        -I ${bamTumor} \
        -V ${germlineResource} \
        -L ${intervalBed} \
        -O ${intervalBed.baseName}_${sampleNameTumor}_pileupsummaries.table
    """
}

pileupSummariesCh = pileupSummariesCh.groupTuple(by:[0,1])

// STEP GATK MUTECT2.4 - MERGING PILEUP SUMMARIES

process MergePileupSummaries {
    label 'gatk'
    label 'cpus_1'

    tag {sampleId + "_" + sampleNameTumor}

    publishDir "${params.outputDir}/VariantCalling/${sampleNameTumor}/Mutect2", mode: params.publishDirMode

    input:
    set pairName, sampleNameTumor, file(pileupSums) from pileupSummariesCh
    file(dict) from dictCh

    output:
    file("${sampleNameTumor}_pileupsummaries.table.tsv") into mergedPileupFileCh

    when: 'mutect2' in tools
    script:
    allPileups = pileupSums.collect{ "-I ${it} " }.join(' ')
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        GatherPileupSummaries \
        --sequence-dictionary ${dict} \
        ${allPileups} \
        -O ${sampleNameTumor}_pileupsummaries.table.tsv
    """
}

// STEP GATK MUTECT2.5 - CALCULATING CONTAMINATION

process CalculateContamination {
    label 'gatk'
    label 'cpus_1'

    tag {sampleNameTumor + "_vs_" + sampleNameNormal}

    publishDir "${params.outputDir}/VariantCalling/${sampleNameTumor}/Mutect2", mode: params.publishDirMode

    input:
    set sampleIdNormal, sampleNameNormal, file(bamNormal), file(baiNormal), sampleIdTumor, sampleNameTumor, file(bamTumor), file(baiTumor) from pairBamCalculateContaminationCh
    file("${sampleNameTumor}_pileupsummaries.table") from mergedPileupFileCh

    output:
    file("${sampleNameTumor}_contamination.table") into contaminationTableCh

    when: 'mutect2' in tools

    script:
    """
    # calculate contamination
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        CalculateContamination \
        -I ${sampleNameTumor}_pileupsummaries.table \
        -O ${sampleNameTumor}_contamination.table
    """
}

// STEP GATK MUTECT2.6 - FILTERING CALLS

process FilterMutect2Calls {
    label 'gatk'
    label 'medCpu'
    label 'medMem'

    tag {sampleNameTN}

    publishDir "${params.outputDir}/VariantCalling/${sampleNameTN}/${"$variantCaller"}", mode: params.publishDirMode

    input:
    set variantCaller, sampleId, sampleNameTN, file(unfiltered), file(unfilteredIndex) from vcfConcatenatedForFilterCh
    file("${sampleNameTN}.vcf.gz.stats") from mergedStatsFileCh
    file("${sampleNameTN}_contamination.table") from contaminationTableCh
    file(dict) from dictCh
    file(fasta) from fastaCh
    file(fastaFai) from fastaFaiCh
    file(germlineResource) from germlineResourceCh
    file(germlineResourceIndex) from germlineResourceIndexCh
    file(intervals) from intervalsCh

    output:
    set val("Mutect2"), sampleId, sampleNameTN,
            file("filtered_${variantCaller}_${sampleNameTN}.vcf.gz"),
            file("filtered_${variantCaller}_${sampleNameTN}.vcf.gz.tbi"),
            file("filtered_${variantCaller}_${sampleNameTN}.vcf.gz.filteringStats.tsv") into filteredMutect2OutputCh

    when: 'mutect2' in tools

    script:
    """
    # do the actual filtering
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        FilterMutectCalls \
        -V ${unfiltered} \
        --contamination-table ${sampleNameTN}_contamination.table \
        --stats ${sampleNameTN}.vcf.gz.stats \
        -R ${fasta} \
        -O filtered_${variantCaller}_${sampleNameTN}.vcf.gz
    """
}


/*
================================================================================
                            SV VARIANT CALLING
================================================================================
*/

// STEP MANTA.1 - SINGLE MODE

process MantaSingle {
    label 'manta'
    label 'cpusMax'
    label 'memoryMax'

    tag {sampleName}

    publishDir "${params.outputDir}/VariantCalling/${sampleName}/Manta", mode: params.publishDirMode

    input:
        set sampleId, sampleName, file(bam), file(bai) from bamMantaSingleCh
        file(fasta) from fastaCh
        file(fastaFai) from fastaFaiCh
        file(targetBED) from targetBEDCh

    output:
        set val("Manta"), sampleId, sampleName, file("*.vcf.gz"), file("*.vcf.gz.tbi") into vcfMantaSingleCh
        file 'v_manta.txt' into mantaSingleVersionCh

    when: 'manta' in tools

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
        Manta_${sampleName}.candidateSmallIndels.vcf.gz
    mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
        Manta_${sampleName}.candidateSmallIndels.vcf.gz.tbi
    mv Manta/results/variants/candidateSV.vcf.gz \
        Manta_${sampleName}.candidateSV.vcf.gz
    mv Manta/results/variants/candidateSV.vcf.gz.tbi \
        Manta_${sampleName}.candidateSV.vcf.gz.tbi
    mv Manta/results/variants/${vcftype}SV.vcf.gz \
        Manta_${sampleName}.${vcftype}SV.vcf.gz
    mv Manta/results/variants/${vcftype}SV.vcf.gz.tbi \
        Manta_${sampleName}.${vcftype}SV.vcf.gz.tbi
    configManta.py --version &> v_manta.txt 2>&1 || true
    """
}

vcfMantaSingleCh = vcfMantaSingleCh.dump(tag:'Single Manta')

// STEP MANTA.2 - SOMATIC PAIR

process Manta {
    label 'manta'
    label 'cpusMax'
    label 'memoryMax'

    tag {sampleNameTumor + "_vs_" + sampleNameNormal}

    publishDir "${params.outputDir}/VariantCalling/${sampleNameTumor}_vs_${sampleNameNormal}/Manta", mode: params.publishDirMode

    input:
        set sampleIdNormal, sampleNameNormal, file(bamNormal), file(baiNormal), sampleIdTumor, sampleNameTumor, file(bamTumor), file(baiTumor) from pairBamMantaCh
        file(fasta) from fastaCh
        file(fastaFai) from fastaFaiCh
        file(targetBED) from targetBEDCh

    output:
        set val("Manta"), pairName, val("${sampleNameTumor}_vs_${sampleNameNormal}"), file("*.vcf.gz"), file("*.vcf.gz.tbi") into vcfMantaCh
        file 'v_manta.txt' into mantaVersionCh

    when: 'manta' in tools

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
        Manta_${sampleNameTumor}_vs_${sampleNameNormal}.candidateSmallIndels.vcf.gz
    mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
        Manta_${sampleNameTumor}_vs_${sampleNameNormal}.candidateSmallIndels.vcf.gz.tbi
    mv Manta/results/variants/candidateSV.vcf.gz \
        Manta_${sampleNameTumor}_vs_${sampleNameNormal}.candidateSV.vcf.gz
    mv Manta/results/variants/candidateSV.vcf.gz.tbi \
        Manta_${sampleNameTumor}_vs_${sampleNameNormal}.candidateSV.vcf.gz.tbi
    mv Manta/results/variants/diploidSV.vcf.gz \
        Manta_${sampleNameTumor}_vs_${sampleNameNormal}.diploidSV.vcf.gz
    mv Manta/results/variants/diploidSV.vcf.gz.tbi \
        Manta_${sampleNameTumor}_vs_${sampleNameNormal}.diploidSV.vcf.gz.tbi
    mv Manta/results/variants/somaticSV.vcf.gz \
        Manta_${sampleNameTumor}_vs_${sampleNameNormal}.somaticSV.vcf.gz
    mv Manta/results/variants/somaticSV.vcf.gz.tbi \
        Manta_${sampleNameTumor}_vs_${sampleNameNormal}.somaticSV.vcf.gz.tbi
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
process AlleleCounter {
    label 'canceritAllelecount'
    label 'memorySingleCPU2Task'

    tag {sampleName}

    input:
        set sampleId, sampleName, file(bam), file(bai) from bamAscatCh
        file(acLoci) from acLociCh
        file(dict) from dictCh
        file(fasta) from fastaCh
        file(fastaFai) from fastaFaiCh

    output:
        set sampleId, sampleName, file("${sampleName}.alleleCount") into alleleCounterOutCh
        file("v_allelecount.txt") into alleleCountsVersionCh

    when: 'ascat' in tools

    script:
    """
    alleleCounter \
        -l ${acLoci} \
        -r ${fasta} \
        -b ${bam} \
        -o ${sampleName}.alleleCount;
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
    [sampleIdNormal, sampleNameNormal, sampleIdTumor, sampleNameTumor, alleleCountOutNormal, alleleCountOutTumor]
}
// STEP ASCAT.2 - CONVERTALLELECOUNTS

// R script from Malin Larssons bitbucket repo:
// https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
process ConvertAlleleCounts {
    label 'ascat'
    label 'memorySingleCPU2Task'

    tag {sampleNameTumor + "_vs_" + sampleNameNormal}

    publishDir "${params.outputDir}/VariantCalling/${sampleNameTumor}_vs_${sampleNameNormal}/ASCAT", mode: params.publishDirMode

    input:
        set sampleIdNormal, sampleNameNormal, sampleIdTumor, sampleNameTumor, file(alleleCountNormal), file(alleleCountTumor) from alleleCounterOutCh

    output:
        set sampleIdNormal, sampleNameNormal, sampleNameTumor, file("${sampleNameNormal}.BAF"), file("${sampleNameNormal}.LogR"), file("${sampleNameTumor}.BAF"), file("${sampleNameTumor}.LogR") into convertAlleleCountsOutCh
        file("v_ascat.txt") into convertAlleleCountsVersionCh

    when: 'ascat' in tools

    script:
    pairName = pairMap[[sampleIdNormal, sampleIdTumor]]
    gender = genderMap[sampleIdNormal]
    """
    Rscript ${workflow.projectDir}/bin/apConvertAlleleCounts.r ${sampleNameTumor} ${alleleCountTumor} ${sampleNameNormal} ${alleleCountNormal} ${gender}
    R -e "packageVersion('ASCAT')" > v_ascat.txt
    """
}

// STEP ASCAT.3 - ASCAT

// R scripts from Malin Larssons bitbucket repo:
// https://bitbucket.org/malinlarsson/somatic_wgs_pipeline
process Ascat {
    label 'ascat'
    label 'memorySingleCPU2Task'

    tag {sampleNameTumor + "_vs_" + sampleNameNormal}

    publishDir "${params.outputDir}/VariantCalling/${sampleNameTumor}_vs_${sampleNameNormal}/ASCAT", mode: params.publishDirMode

    input:
        set sampleIdNormal, sampleNameNormal, sampleNameTumor, file(bafNormal), file(logrNormal), file(bafTumor), file(logrTumor) from convertAlleleCountsOutCh
        file(acLociGC) from acLociGCCh

    output:
        set val("ASCAT"), sampleIdNormal, sampleNameNormal, sampleNameTumor, file("${sampleNameTumor}.*.{png,txt}") into ascatOutCh
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
        --tumorname ${sampleNameTumor} \
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
            [variantcaller, sampleName, vcf]
    },
    vcfConcatenatedCh.map{
        variantcaller, sampleId, sampleName, vcf, tbi ->
            [variantcaller, sampleName, vcf]
    },
    vcfMantaSingleCh.map {
        variantcaller, sampleId, sampleName, vcf, tbi ->
        [variantcaller, sampleName, vcf[2]]
    },
    vcfMantaDiploidSVCh.map {
        variantcaller, sampleId, sampleName, vcf, tbi ->
        [variantcaller, sampleName, vcf[2]]
    },
    vcfMantaSomaticSVCh.map {
        variantcaller, sampleId, sampleName, vcf, tbi ->
        [variantcaller, sampleName, vcf[3]]
    })

if (step == 'annotate') {
    vcfToAnnotateCh = Channel.create()
    vcfNoAnnotateCh = Channel.create()

    if (inputPath == []) {
        // By default, annotates all available vcfs that it can find in the VariantCalling directory
        // Excluding vcfs from and g.vcf from HaplotypeCaller
        // Basically it's: results/VariantCalling/*/{HaplotypeCaller,Manta,Mutect2}/*.vcf.gz
        // Without *SmallIndels.vcf.gz from Manta
        // The small snippet `vcf.minus(vcf.fileName)[-2]` catches sampleName
        // This field is used to output final annotated VCFs in the correct directory
        Channel.empty().mix(
            Channel.fromPath("${params.outputDir}/VariantCalling/*/HaplotypeCaller/*.vcf.gz")
              .flatten().map{vcf -> ['HaplotypeCaller', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
            Channel.fromPath("${params.outputDir}/VariantCalling/*/Manta/*[!candidate]SV.vcf.gz")
              .flatten().map{vcf -> ['Manta', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
            Channel.fromPath("${params.outputDir}/VariantCalling/*/Mutect2/*.vcf.gz")
              .flatten().map{vcf -> ['Mutect2', vcf.minus(vcf.fileName)[-2].toString(), vcf]}
        ).choice(vcfToAnnotateCh, vcfNoAnnotateCh) {
        annotateTools == [] || (annotateTools != [] && it[0] in annotateTools) ? 0 : 1
      }
    } else if (annotateTools == []) {
    // Annotate user-submitted VCFs
    // If user-submitted, assume that the sampleName should be assumed automatically
      vcfToAnnotateCh = Channel.fromPath(inputPath)
        .map{vcf -> ['userspecified', vcf.minus(vcf.fileName)[-2].toString(), vcf]}
    } else exit 1, "specify only tools or files to annotate, not both"

    vcfNoAnnotateCh.close()
    vcfAnnotationCh = vcfAnnotationCh.mix(vcfToAnnotateCh)
}
// as now have the list of VCFs to annotate, the first step is to annotate with allele frequencies, if there are any

// STEP SNPEFF

process Snpeff {
    tag {"${sampleName} - ${variantCaller} - ${vcf}"}
    label 'snpeff'

    publishDir params.outputDir, mode: params.publishDirMode, saveAs: {
        if (it == "${reducedVCF}_snpEff.ann.vcf") null
        else "Reports/${sampleName}/snpEff/${it}"
    }

    input:
        set variantCaller, sampleName, file(vcf) from vcfAnnotationCh
        file(dataDir) from snpEffCacheCh
        val snpeffDb from snpeffDbCh

    output:
        set file("${reducedVCF}_snpEff.txt"), file("${reducedVCF}_snpEff.html"), file("${reducedVCF}_snpEff.csv") into snpeffReportCh
        set variantCaller, sampleName, file("${reducedVCF}_snpEff.ann.vcf") into snpeffVCFCh
        file 'v_snpeff.txt' into snpeffVersionCh

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

process CompressVCFsnpEff {
    tag {"${sampleName} - ${vcf}"}
    label 'tabix'

    publishDir "${params.outputDir}/Annotation/${sampleName}/snpEff", mode: params.publishDirMode

    input:
        set variantCaller, sampleName, file(vcf) from snpeffVCFCh

    output:
        set variantCaller, sampleName, file("*.vcf.gz"), file("*.vcf.gz.tbi") into (compressVCFsnpEffOutCh)

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
process GetSoftwareVersions {
    label 'python'

    publishDir path:"${params.outputDir}/pipeline_info", mode: params.publishDirMode

    input:
        file 'v_ascat.txt' from ascatVersionCh.mix(convertAlleleCountsVersionCh).first().ifEmpty('')
        file 'v_allelecount.txt' from alleleCountsVersionCh.first().ifEmpty('')
        file 'v_bcftools.txt' from bcftoolsVersionCh.first().ifEmpty('')
        file 'v_bwa.txt' from bwaVersionCh.ifEmpty('')
        file 'v_fastqc.txt' from fastqcVersionCh.mix(fastqcBamVersionCh).first().ifEmpty('')
        file 'v_gatk.txt' from gatkVersionCh.first().ifEmpty('')
        file 'v_manta.txt' from mantaVersionCh.mix(mantaSingleVersionCh).first().ifEmpty('')
        file 'v_qualimap.txt' from qualimapVersionCh.first().ifEmpty('')
        file 'v_samtools.txt' from samtoolsBwaMemUniqVersionCh.mix(samtoolsIndexBamFileVersionCh).mix(samtoolsIndexBamRecalVersionCh).mix(samtoolsMapQVersionCh).mix(samtoolsMapReadsVersionCh).mix(samtoolsMergeBamMappedVersionCh).mix(samtoolsMergeBamRecalVersionCh).mix(samtoolsStatsVersionCh).first().ifEmpty('')
        file 'v_snpeff.txt' from snpeffVersionCh.first().ifEmpty('')

    output:
        file 'software_versions_mqc.yaml' into yamlSoftwareVersionCh

    when: !('versions' in skipQC)

    script:
    """
    echo "${workflow.manifest.version}" &> v_pipeline.txt 2>&1 || true
    echo "${workflow.nextflow.version}" &> v_nextflow.txt 2>&1 || true

    apScrapeSoftwareVersions.py &> software_versions_mqc.yaml
    """
}

yamlSoftwareVersionCh = yamlSoftwareVersionCh.dump(tag:'SOFTWARE VERSIONS')

process MultiQC {
    label 'multiqc'
    publishDir "${params.outputDir}/Reports/MultiQC", mode: params.publishDirMode

    input:
        file multiqcConfig from Channel.value(params.multiqcConfig ? file(params.multiqcConfig) : "")
        file workflow_summary from workflowSummaryCh.collectFile(name: "workflow_summary_mqc.yaml")
        file (versions) from yamlSoftwareVersionCh
        file ('bamQC/*') from bamQCReportCh.collect().ifEmpty([])
        file ('FastQC/*') from fastQCReportCh.collect().ifEmpty([])
        file ('MarkDuplicates/*') from markDuplicatesReportCh.collect().ifEmpty([])
        file ('SamToolsStats/*') from samtoolsStatsReportCh.collect().ifEmpty([])
        file ('snpEff/*') from snpeffReportCh.collect().ifEmpty([])

    output:
        file "*multiqc_report.html" into multiQCOutCh
        file "*_data"

    when: !('multiqc' in skipQC)

    script:
    rtitle = customRunName ? "--title \"$customRunName\"" : ''
    rfilename = customRunName ? "--filename " + customRunName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f ${rtitle} ${rfilename} --config ${multiqcConfig} .
    """
}

multiQCOutCh.dump(tag:'MultiQC')

workflow.onComplete {

    def reportFields = params + [
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
