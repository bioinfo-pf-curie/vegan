/*
 * QC Worflow for BAM files
 */

include { collectInsertSizeMetrics } from '../../common/process/picard/collectInsertSizeMetrics'
include { mosdepth } from '../../common/process/mosdepth/mosdepth'
include { collectWgsMetrics } from '../../common/process/gatk/collectWgsMetrics'

workflow bamQcFlow {

  take:
  bam
  bed
  fasta
  fai
  dict

  main:
  chVersions = Channel.empty()

  // Insert size for paired-end data
  collectInsertSizeMetrics(
    bam,
    fasta,
    fai
  )
  chVersions = chVersions.mix(collectInsertSizeMetrics.out.versions)

  // Sequencing depth
  mosdepth(
    bam,
    bed,
    fasta
  )
  chVersions = chVersions.mix(mosdepth.out.versions)
  chMosdepthLog = params.targetBed ? mosdepth.out.regionsTxt.mix(mosdepth.out.regionsBed) : mosdepth.out.globalTxt.mix(mosdepth.out.regionsBed)

  // WGS metrics
  collectWgsMetrics(
    bam,
    bed,
    fasta,
    fai,
    dict
  )
  chVersions = chVersions.mix(collectWgsMetrics.out.versions)

  emit:
  fragSize = collectInsertSizeMetrics.out.results
  depth = chMosdepthLog
  wgsMetrics = collectWgsMetrics.out.metrics
  versions = chVersions
}
