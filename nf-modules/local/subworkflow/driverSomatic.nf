/*
 * Apply decision algorithm to select drivers mutation
 */

include { oncoDriver } from '../../local/process/oncoDriver'
include { tabix } from '../../common/process/tabix/tabix'

workflow driverSomaticFlow {

  take:
  snv // [meta, vcf, tbi]
  cnv // [meta, cnv]
  drivers
  decision
  canonicals

  main:
  chVersions = Channel.empty()

  tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase()} : []
  if('oncodriver' in tools){
    if('facets' in tools && 'snpEff' in tools){
      chVariants = snv
        .combine(cnv, by:0)
    }else if ('snpEff' in tools){
      chVariants = snv
        .map{meta, vcf, tbi -> [meta, vcf, tbi, []]}
    }else if ('facets' in tools){
      chVariants = cnv
        .map{meta, seg -> [meta, [], [], seg]}
    }
  }

  // Detect driver mutations
  oncoDriver(
    chVariants,
    drivers,
    decision,
    canonicals
  )
  chVersions = chVersions.mix(oncoDriver.out.versions)

  tabix(
    oncoDriver.out.vcf
  )
  chVersions = chVersions.mix(tabix.out.versions)

  emit:
  versions = chVersions
  vcf = oncoDriver.out.vcf.join(tabix.out.tbi)
}