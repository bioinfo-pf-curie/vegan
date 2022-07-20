/*
 * VCF annotation sub-workflow
 */

include { snpEff } from '../../common/process/snpEff/snpEff'
include { tabix } from '../../common/process/tabix/tabix'

workflow vcfAnnotFlow {

  take:
  vcf
  db
  cache

  main:
  chVersions = Channel.empty()

  snpEff(
    vcf,
    db,
    cache
  )
  chVersions = chVersions.mix(snpEff.out.versions)

  tabix(
    snpEff.out.vcf
  )
  chVersions = chVersions.mix(tabix.out.versions)

  emit:
  vcf = tabix.out.vcf
  logs = snpEff.out.report
  versions = chVersions
}
