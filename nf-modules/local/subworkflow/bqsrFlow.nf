/*
 * BQSR step
 */

include { baseRecalibrator } from '../../local/process/gatk/baseRecalibrator'
include { applyBQSR } from '../../local/process/gatk/applyBqsr'
include { samtoolsIndex } from '../../common/process/samtoolsIndex'
//include { indexBamRecal } from '../../local/process/indexBamRecal'

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
