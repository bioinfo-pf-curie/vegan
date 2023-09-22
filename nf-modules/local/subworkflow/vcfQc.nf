/*
 * QC Worflow for VCF files
 */

include { collectVCFmetrics } from '../../local/process/collectVCFmetrics'
include { computeTransition } from '../../local/process/computeTransition'

workflow vcfQcFlow {

  take:
  vcfs //meta, vcf_raw, tbi_raw, vcf_filt, tbi_filt

  main:
  chVersions = Channel.empty()

  collectVCFmetrics(
    vcfs
  )

  computeTransition(
    vcfs.map{ it -> [it[0], it[3], it[4]]}
  )

  emit:
  transition = computeTransition.out.metrics
  mqc = collectVCFmetrics.out.mqc
  versions = chVersions
}
