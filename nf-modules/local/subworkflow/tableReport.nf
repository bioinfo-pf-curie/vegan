/*
 * VCF to table report sub-workflow
 */

include { snpSiftExtractFields } from '../../common/process/snpSift/snpSiftExtractFields'

workflow tableReportFlow {

  take:
  vcf

  main:
  chVersions = Channel.empty()

  snpSiftExtractFields(
    vcf
  )

  emit:
  tsv = snpSiftExtractFields.out.tsv
  versions = chVersions
}
