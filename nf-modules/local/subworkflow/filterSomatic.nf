/*
 * Apply basic filtering on somatic variants
 */

include { bcftoolsNorm } from '../../common/process/bcftools/bcftoolsNorm'
include { snpSiftAnnotate as snpSiftGnomAD } from '../../common/process/snpSift/snpSiftAnnotate'
include { bcftoolsFilter } from '../../common/process/bcftools/bcftoolsFilter'

workflow filterSomaticFlow {

  take:
  vcf // [meta, vcf, tbi]
  fasta 
  gnomAd
  gnomAdIndex

  main:
  annotDb = params.annotDb ? params.annotDb.split(',').collect{it.trim().toLowerCase()} : []
  chVersions = Channel.empty()

  // Normalize VCFs
  bcftoolsNorm(
    vcf,
    fasta
  )
  chVersions = chVersions.mix(bcftoolsNorm.out.versions)

  // Annotate with gnomAD
  snpSiftGnomAD(
    bcftoolsNorm.out.vcf,
    gnomAd.combine(gnomAdIndex).collect()
  )
  chVersions = chVersions.mix(snpSiftGnomAD.out.versions)
  chVcf2Filter = (params.annotDb && params.annotDb.contains('gnomad')) ? snpSiftGnomAD.out.vcf : bcftoolsNorm.out.vcf
  
  // Filter VCF with PASS/DP/AF
  bcftoolsFilter(
    chVcf2Filter
  )
  chVersions = chVersions.mix(bcftoolsFilter.out.versions)

  emit:
  versions = chVersions
  vcfRaw = bcftoolsNorm.out.vcf
  vcfFiltered = bcftoolsFilter.out.vcf
}
