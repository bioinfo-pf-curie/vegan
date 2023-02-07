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

    baseRecalibrator(
      bam.combine(intervals),
      dbsnp,
      dbsnpIndex,
      fasta,
      fastaFai,
      knownIndels,
      knownIndelsIndex,
      dict
    )

    if (!params.noIntervals && !params.targetBed){
      gatherBQSRReports(
        baseRecalibrator.out.table.groupTuple()
      )
      chBaseRecalTable = gatherBQSRReports.out.table
    }else{
      chBaseRecalTable = baseRecalibrator.out.table
    }

    applyBQSR(
      bam.join(chBaseRecalTable).combine(intervals),
      fasta,
      fastaFai,
      dict
    )

    if (!params.noIntervals && !params.targetBed){
      samtoolsMerge(
        applyBQSR.out.bam.groupTuple()
      )
      chBamRecal = samtoolsMerge.out.bam
    }else{
      chBamRecal = applyBQSR.out.bam
    }

    samtoolsIndex(
      chBamRecal
    )

  emit:
  table = baseRecalibrator.out.table
  bam = chBamRecal.join(samtoolsIndex.out.bai)
  versions = applyBQSR.out.versions
}
