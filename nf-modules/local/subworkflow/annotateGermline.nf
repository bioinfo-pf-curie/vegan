/*
 * VCF germline annotation sub-workflow
 */

include { snpEff } from '../../common/process/snpEff/snpEff'
include { snpSiftAnnotate as snpSiftCosmic } from '../../common/process/snpSift/snpSiftAnnotate'
include { snpSiftAnnotate as snpSiftIcgc } from '../../common/process/snpSift/snpSiftAnnotate'
include { snpSiftAnnotate as snpSiftCancerHotspots } from '../../common/process/snpSift/snpSiftAnnotate'
include { snpSiftAnnotate as snpSiftGnomAD } from '../../common/process/snpSift/snpSiftAnnotate'
include { snpSiftDbnsfp } from '../../common/process/snpSift/snpSiftDbnsfp'

workflow annotateGermlineFlow {

  take:
  vcf
  snpEffDb
  snpEffCache
  gnomAd
  gnomAdIndex
  dbnsfp
  dbnsfpIndex

  main:
  annotDb = params.annotDb ? params.annotDb.split(',').collect{it.trim().toLowerCase()} : []
  chVersions = Channel.empty()

  /*
   * Main annotation with snpEff
   */

  snpEff(
    vcf,
    snpEffDb,
    snpEffCache
  )
  chVersions = chVersions.mix(snpEff.out.versions)
  chAnnotVcf = snpEff.out.vcf

  /*
   * GnomAD annotations
   */

  snpSiftGnomAD(
    chAnnotVcf,
    gnomAd.combine(gnomAdIndex).collect()
  )
  chVersions = chVersions.mix(snpSiftGnomAD.out.versions)
  chAnnotVcf = 'gnomad' in annotDb ? snpSiftGnomAD.out.vcf : chAnnotVcf

  /*
   * SnpSift dbNSFP
   */

  snpSiftDbnsfp(
    chAnnotVcf,
    dbnsfp.combine(dbnsfpIndex).collect()
  )
  chVersions = chVersions.mix(snpSiftDbnsfp.out.versions)
  chAnnotVcf = 'dbnsfp' in annotDb ? snpSiftDbnsfp.out.vcf : chAnnotVcf

  emit:
  vcf = chAnnotVcf
  snpEffReport = snpEff.out.report
  versions = chVersions
}
