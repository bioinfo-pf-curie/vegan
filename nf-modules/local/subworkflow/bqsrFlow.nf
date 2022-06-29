/*
 * BQSR step
 */

include { baseRecalibrator } from '../../common/process/gatk/baseRecalibrator'
include { applyBQSR } from '../../common/process/gatk/applyBqsr'
include { samtoolsIndex } from '../../common/process/samtools/samtoolsIndex'

workflow bqsrFlow {

  take:
    bam
    bed
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
      bam,
      bed,
      dbsnp,
      dbsnpIndex,
      fasta,
      fastaFai,
      knownIndels,
      knownIndelsIndex,
      dict
    )

    applyBQSR(
      baseRecalibrator.out.table,
      bed,
      fasta,
      fastaFai,
      dict
    )

    samtoolsIndex(
      applyBQSR.out.bam
    )

  emit:
    bqsrTable = baseRecalibrator.out.table
    bqsrBam = applyBQSR.out.bam.join(samtoolsIndex.out.bai)
    versions = applyBQSR.out.versions
}
