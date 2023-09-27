/*
 * MSI
 */

include { prepareBaselineConfig } from '../../common/process/msisensorpro/prepareBaselineConfig.nf'
include { msisensorproScan } from '../../common/process/msisensorpro/msisensorproScan.nf'
include { msisensorproMsi } from '../../common/process/msisensorpro/msisensorproMsi.nf'
include { msisensorproBaseline } from '../../common/process/msisensorpro/msisensorproBaseline.nf'
include { msisensorproPro } from '../../common/process/msisensorpro/msisensorproPro.nf'

workflow msiFlow {

  take:
  pairedBam
  tumorOnlyBam
  baselineConfig
  fasta
  bed

  main:
  chVersions = Channel.empty()

  msisensorproScan(
    fasta
  )

  /* Tumor only BAM files */
  msisensorproBaseline(
    msisensorproScan.out.list.collect(),
    baselineConfig
  )
  chVersions = chVersions.mix(msisensorproBaseline.out.versions)

  msisensorproPro(
    tumorOnlyBam,
    msisensorproBaseline.out.baseline,
    bed
  )
  chVersions = chVersions.mix(msisensorproPro.out.versions)

  /* Paired BAM files */

  msisensorproMsi(
    pairedBam,
    fasta.collect(),
    msisensorproScan.out.list.collect(),
    bed
  )
  chVersions = chVersions.mix(msisensorproMsi.out.versions)

  emit:
  report = msisensorproMsi.out.outputReport.mix(msisensorproPro.out.outputReport)
  versions = chVersions
}
