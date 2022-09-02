/*
 * MSI
 */

include { msisensorproScan } from '../../common/process/msisensorpro/msisensorproScan.nf'
include { msisensorproMsi } from '../../common/process/msisensorpro/msisensorproMsi.nf'

workflow msiFlow {

  take:
  bam
  fasta
  bed

  main:
  chVersions = Channel.empty()

  msisensorproScan(
    fasta
  )

  msisensorproMsi(
    bam,
    fasta.collect(),
    msisensorproScan.out.list.collect(),
    bed
  )
  chVersions = chVersions.mix(msisensorproMsi.out.versions)

  emit:
  logs = msisensorproMsi.out.outputReport
  versions = chVersions
}
