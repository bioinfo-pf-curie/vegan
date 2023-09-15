/*
 * Filtering subworkflow
 */

include { bcftoolsNorm } from '../../common/process/bcftools/bcftoolsNorm'
include { snpSiftAnnotate as snpSiftGnomAD } from '../../common/process/snpSift/snpSiftAnnotate'
// include { filter as filterDp } from '../../common/process/bcftools/bcftoolsFilter'
// include { filter as filterFreq } from '../../common/process/bcftools/bcftoolsFilter'
// include { filter as filterMaf } from '../../common/process/bcftools/bcftoolsFilter'
include { bcftoolsFilter } from '../../common/process/bcftools/bcftoolsFilter'
// include { filter as keepPassOnly } from '../../common/process/bcftools/bcftoolsFilter'

workflow filterSomaticFlow {

  take:
  vcf // [meta, vcf, tbi]
  fasta 
  gnomAd
  gnomAdIndex

  main:
  annotDb = params.annotDb ? params.annotDb.split(',').collect{it.trim().toLowerCase()} : []
  chVersions = Channel.empty()

  // Normalization
  bcftoolsNorm(
    vcf,
    fasta
  )
  chVersions = chVersions.mix(bcftoolsNorm.out.versions)
  chFilterVcf = bcftoolsNorm.out.vcf

  // Annotation with gnomAD
  snpSiftGnomAD(
    chFilterVcf,
    gnomAd.combine(gnomAdIndex).collect()
  )
  chVersions = chVersions.mix(snpSiftGnomAD.out.versions)
  chFilterVcf = (params.annotDb && params.annotDb.contains('gnomad')) ? snpSiftGnomAD.out.vcf : chFilterVcf
  
  // Filtering on DP,FREQ MAF via bcftools
  // filterDp(
  //   chFilterVcf
  // )
  // chFilterVcf = params.filterDp != 0 ? filterDp.out.vcf : chFilterVcf

  // filterFreq(
  //   chFilterVcf
  // )
  // chFilterVcf = params.filterFreq != 0 ? filterFreq.out.vcf : chFilterVcf

  // filterMaf(
  //   chFilterVcf
  // )
  // chVersions = chVersions.mix(filterMaf.out.versions)
  // chFilterVcf = (params.filterMaf &&  params.annotDb && params.annotDb.contains('gnomad')) ? filterMaf.out.vcf : chFilterVcf

  bcftoolsFilter(
    chFilterVcf
  )
  chVersions = chVersions.mix(bcftoolsFilter.out.versions)
  chFilterVcf = bcftoolsFilter.out.vcf


  // // Keep only PASS variants
  // keepPassOnly(
  //   chFilterVcf
  // )
  // chFilterVcf = keepPassOnly.out.vcf

  emit:
  versions = chVersions
  vcfRaw = bcftoolsNorm.out.vcf
  vcfPass = chFilterVcf
  vcfFiltered = chFilterVcf
}
