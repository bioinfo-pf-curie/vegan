/*
 * BQSR step
 */

include { baseRecalibrator } from '../../common/process/gatk/baseRecalibrator'
include { gatherBQSRReports } from '../../common/process/gatk/gatherBQSRReports'
include { applyBQSR } from '../../common/process/gatk/applyBqsr'
include { samtoolsIndex } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsMerge } from '../../common/process/samtools/samtoolsMerge'

workflow bqsrFlow {

  take:
    bam
    intervals
    dbsnp
    dbsnpIndex
    fasta
    fastaFai
    knownIndels
    knownIndelsIndex
    dict

  main:
    chVersions = Channel.empty()
    chBamIntervals = params.noIntervals ? bam.map{ meta,bam,bai -> [meta,bam,bai,[]]} : bam.combine(intervals)

    baseRecalibrator(
      chBamIntervals,
      dbsnp,
      dbsnpIndex,
      fasta,
      fastaFai,
      knownIndels,
      knownIndelsIndex,
      dict
    )

    gatherBQSRReports(
      baseRecalibrator.out.table.groupTuple()
    )

    chBaseRecalTable = params.noIntervals ? bam.join(baseRecalibrator.out.table).map{ meta,bam,bai,table -> [meta,bam,bai,table,[]]} :
                       params.targetBed ? bam.join(baseRecalibrator.out.table).combine(intervals) : bam.join(gatherBQSRReports.out.table).combine(intervals)

    applyBQSR(
      chBaseRecalTable,
      fasta,
      fastaFai,
      dict
    )

    samtoolsMerge(
      applyBQSR.out.bam.groupTuple()
    )

    chBamBQSR = params.noIntervals || params.targetBed ? applyBQSR.out.bam : samtoolsMerge.out.bam

    samtoolsIndex(
      chBamBQSR
    )

  emit:
  table = params.noIntervals ? baseRecalibrator.out.table : gatherBQSRReports.out.table
  bam = chBamBQSR.join(samtoolsIndex.out.bai)
  versions = applyBQSR.out.versions
}
