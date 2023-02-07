/*
 * QC Worflow for VCF files
 */

include { collectVCFmetrics } from '../../local/process/collectVCFmetrics'
include { computeTransition } from '../../local/process/computeTransition'

workflow vcfQcFlow {

  take:
  vcf //[meta, vcf, tbi]

  main:
  chVersions = Channel.empty()

  collectVCFmetrics(
    vcf.map{ it -> [it[0], it[1], it[2], [], [], []] }
  )

  computeTransition(
    vcf
  )

  emit:
  transition = computeTransition.out.metrics
  mqc = collectVCFmetrics.out.mqc
  versions = chVersions
}
