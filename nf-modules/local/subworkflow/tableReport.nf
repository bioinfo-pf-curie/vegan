/*
 * VCF to table report sub-workflow
 */

include { snpSiftExtractFieldsPatch } from '../../common/process/snpSift/snpSiftExtractFieldsPatch'

workflow tableReportFlow {

  take:
  vcf

  main:
  chVersions = Channel.empty()

  snpSiftExtractFieldsPatch(
    vcf
  )

  emit:
  tsv = snpSiftExtractFieldsPatch.out.tsv
  versions = chVersions
}
