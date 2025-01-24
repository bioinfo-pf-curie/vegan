/*
 * Apply basic filtering on somatic variants
 */

include { bcftoolsNorm } from '../../common/process/bcftools/bcftoolsNorm'
include { snpSiftAnnotate as snpSiftGnomAD } from '../../common/process/snpSift/snpSiftAnnotate'
include { bcftoolsFilter } from '../../common/process/bcftools/bcftoolsFilter'
include { bcftoolsFilter as rescueVariants } from '../../common/process/bcftools/bcftoolsFilter'
include { bcftoolsConcat } from '../../common/process/bcftools/bcftoolsConcat'

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

  // Rescue non-PASS variants of interest
  rescueVariants(
    chVcf2Filter
  )

  // Concat filtered vcf and rescue ones
  bcftoolsConcat(
    bcftoolsNorm.out.vcf.join(rescueVariants.out.vcf).map{it->[it[0],[it[1][0],it[2][0]],[it[1][1],it[2][1]]]}
  )

  vcfFiltered = params.rescueVariants ? bcftoolsConcat.out.vcf : bcftoolsFilter.out.vcf

  emit:
  versions = chVersions
  vcfRaw = bcftoolsNorm.out.vcf
  vcfFiltered = vcfFiltered
}
