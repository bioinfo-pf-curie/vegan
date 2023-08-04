/*
 * QC Worflow for VCF files
 */

include { collectVCFmetrics } from '../../local/process/collectVCFmetrics'
include { computeTransition } from '../../local/process/computeTransition'

workflow vcfQcFlow {

  take:
  allVcf

  main:
  chVersions = Channel.empty()

  collectVCFmetrics(
    allVcf
  )

  computeTransition(
    allVcf.map{ it -> [it[0], it[3], it[4]]}
  )

  emit:
  transition = computeTransition.out.metrics
  mqc = collectVCFmetrics.out.mqc
  versions = chVersions
}
