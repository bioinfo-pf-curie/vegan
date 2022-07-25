/*
 * VCF annotation sub-workflow
 */

include { snpEff } from '../../common/process/snpEff/snpEff'
include { snpSiftAnnotate as snpSiftCosmic } from '../../common/process/snpSift/snpSiftAnnotate'
include { snpSiftAnnotate as snpSiftIcgc } from '../../common/process/snpSift/snpSiftAnnotate'
include { snpSiftAnnotate as snpSiftCancerHotspots } from '../../common/process/snpSift/snpSiftAnnotate'
include { snpSiftAnnotate as snpSiftGnomAD } from '../../common/process/snpSift/snpSiftAnnotate'
include { snpSiftDbnsfp } from '../../common/process/snpSift/snpSiftDbnsfp'

workflow annotateSomaticFlow {

  take:
  vcf
  snpEffDb
  snpEffCache
  cosmic
  cosmicIndex
  icgc
  icgcIndex
  cancerHotspots
  cancerHotspotsIndex
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
   * COSMIC annotations
   */

  snpSiftCosmic(
    chAnnotVcf,
    cosmic.combine(cosmicIndex)
  )
  chVersions = chVersions.mix(snpSiftCosmic.out.versions)
  chAnnotVcf = 'cosmic' in annotDb ? snpSiftCosmic.out.vcf : chAnnotVcf

  /*
   * ICGC annotations
   */
  
  icgc.view()
  icgcIndex.view()

  snpSiftIcgc(
    chAnnotVcf,
    icgc.combine(icgcIndex)
  )
  chVersions = chVersions.mix(snpSiftIcgc.out.versions)
  chAnnotVcf = 'icgc' in annotDb ? snpSiftIcgc.out.vcf : chAnnotVcf

  /*
   * CancerHotspots annotations
   */

  snpSiftCancerHotspots(
    chAnnotVcf,
    cancerHotspots.combine(cancerHotspotsIndex)
  )
  chVersions = chVersions.mix(snpSiftCancerHotspots.out.versions)
  chAnnotVcf = 'cancerhotspots' in annotDb ? snpSiftCancerHotspots.out.vcf : chAnnotVcf

  /*
   * GnomAD annotations
   */

  snpSiftGnomAD(
    chAnnotVcf,
    gnomAd.combine(gnomAdIndex)
  )
  chVersions = chVersions.mix(snpSiftGnomAD.out.versions)
  chAnnotVcf = 'gnomad' in annotDb ? snpSiftGnomAD.out.vcf : chAnnotVcf

  /*
   * SnpSift dbNSFP
   */

  snpSiftDbnsfp(
    chAnnotVcf,
    dbnsfp.combine(dbnsfpIndex)
  )
  chVersions = chVersions.mix(snpSiftDbnsfp.out.versions)
  chAnnotVcf = 'dbnsfp' in annotDb ? snpSiftDbnsfp.out.vcf : chAnnotVcf

  emit:
  vcf = chAnnotVcf
  versions = chVersions
}
