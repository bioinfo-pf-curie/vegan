/*
 * VCF annotation sub-workflow
 */

include { snpEff } from '../../common/process/snpEff/snpEff'
include { snpSiftAnnotate as snpSiftCosmic } from '../../common/process/snpSift/snpSiftAnnotate'
include { snpSiftAnnotate as snpSiftIcgc } from '../../common/process/snpSift/snpSiftAnnotate'
include { snpSiftAnnotate as snpSiftCancerHotspot } from '../../common/process/snpSift/snpSiftAnnotate'
include { snpSiftAnnotate as snpSiftGnomAD } from '../../common/process/snpSift/snpSiftAnnotate'
include { snpSiftDbnsfp } from '../../common/process/snpSift/snpSiftDbnsfp'

workflow annotateSomaticFlow {

  take:
  vcf
  snpEffDb
  snpEffCache
  cosmic
  icgc
  cancerHotspot
  gnomAd

  main:
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
   * COSMIC annotations
   */

  snpSiftCosmic(
    chAnnotVcf,
    cosmic
  )
  chVersions = chVersions.mix(snpSiftCosmic.out.versions)
  chAnnotVcf = params.annotDb.contains('cosmic') ? snpSiftCosmic.out.vcf : chAnnotVcf

  /*
   * ICGC annotations
   */
  
  snpSiftIcgc(
    chAnnotVcf,
    icgc
  )
  chVersions = chVersions.mix(snpSiftIcgc.out.versions)
  chAnnotVcf = params.annotDb.contains('icgc') ? snpSiftIcgc.out.vcf : chAnnotVcf

  /*
   * CancerHotspot annotations
   */

  snpSiftCancerHotspot(
    chAnnotVcf,
    cancerHotspot
  )
  chVersions = chVersions.mix(snpSiftCancerHotspot.out.versions)
  chAnnotVcf = params.annotDb.contains('cancerhotspot') ? snpSiftCancerHotspot.out.vcf : chAnnotVcf

  /*
   * GnomAD annotations
   */

  snpSiftGnomAD(
    chAnnotVcf,
    gnomAD
  )
  chVersions = chVersions.mix(snpSiftGnomAD.out.versions)
  chAnnotVcf = params.annotDb.contains('gnomad') ? snpSiftGnomAD.out.vcf : chAnnotVcf

  /*
   * SnpSift dbNSFP
   */

  snpSiftDbnsfp(
    chAnnotVcf
  )
  chVersions = chVersions.mix(snpSiftDbnsfp.out.versions)
  chAnnotVcf = params.annotDb.contains('dbnsfp') ? snpSiftDbnsfp.out.vcf : chAnnotVcf

  emit:
  vcf = chAnnotVcf
  versions = chVersions
}
