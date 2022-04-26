/*
 * BQSR step
 */

include { baseRecalibrator } from '../../local/process/baseRecalibrator'
include { applyBQSR } from '../../local/process/applyBqsr'
include { indexBamRecal } from '../../local/process/indexBamRecal'

workflow bqsrFlow {

  take:
    bamFiltered
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
      bamFiltered,
      bed.collect(),
      dbsnp.collect(),
      dbsnpIndex.collect(),
      fasta.collect(),
      fastaFai.collect(),
      knownIndels.collect(),
      knownIndelsIndex.collect(),
      dict.collect()
      )

    applyBQSR(
      baseRecalibrator.out.table,
      bed.collect(),
      fasta.collect(),
      fastaFai.collect(),
      dict.collect(),
      )

    if(params.noIntervals){indexBamRecal(
      applyBQSR.out.bqsrBam
      )
    }

  emit:
    bqsrTable = baseRecalibrator.out.table
    bqsrBam = applyBQSR.out.bqsrBam
}
