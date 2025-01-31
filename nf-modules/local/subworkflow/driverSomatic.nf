/*
 * Apply decision algorithm to select drivers mutation
 */

include { oncoDriver } from '../../local/process/oncoDriver'
include { tabix } from '../../common/process/tabix/tabix'

workflow driverSomaticFlow {

  take:
  variants // [meta, vcf, tbi, cnv]
  drivers
  decision
  canonicals

  main:
  chVersions = Channel.empty()

  // Detect driver mutations
  oncoDriver(
    variants,
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