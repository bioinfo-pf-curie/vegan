/*
 * Filtering subworkflow
 */

include { bcftoolsNorm } from '../../common/process/bcftools/bcftoolsNorm'
include { filterVcf as keepPassOnly } from '../../local/process/filterVcf'
include { snpSiftAnnotate as snpSiftGnomAD } from '../../common/process/snpSift/snpSiftAnnotate'
include { filterDpFreqMaf } from '../../common/process/bcftools/filterDpFreqMaf'

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

  // Keep only PASS variants
  keepPassOnly(
    bcftoolsNorm.out.vcf
  )
  chFilterVcf = keepPassOnly.out.vcf.map{it -> [it[0], [it[1],it[2]]] }

  // Annotation with gnomAD
  snpSiftGnomAD(
    chFilterVcf,
    gnomAd.combine(gnomAdIndex).collect()
  )
  chVersions = chVersions.mix(snpSiftGnomAD.out.versions)

  chFilterVcf = ('gnomad' in annotDb) ? snpSiftGnomAD.out.vcf : chFilterVcf
  chFilterVcf = params.filter ? chFilterVcf : chFilterVcf.map{it -> [it[0],it[1][0],it[1][1]]}

  // Filtering on DP,FREQ MAF via bcftools
  filterDpFreqMaf(
    chFilterVcf
  )
  chVersions = chVersions.mix(filterDpFreqMaf.out.versions)
  chFilterVcf = params.filter ? filterDpFreqMaf.out.vcf : chFilterVcf

  emit:
  versions = chVersions
  vcfRaw = bcftoolsNorm.out.vcf
  vcfPass = keepPassOnly.out.vcf
  vcfFiltered = chFilterVcf
}