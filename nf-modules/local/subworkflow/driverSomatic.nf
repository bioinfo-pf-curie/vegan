/*
 * Apply decision algorithm to select drivers mutation
 */

include { oncoDriver } from '../../local/process/oncoDriver'
include { annotateCNV } from '../../local/process/annotateCNV'
include { tabix } from '../../common/process/tabix/tabix'

workflow driverSomaticFlow {

  take:
  variants // [meta, vcf, tbi, cnv]
  gtf
  drivers
  decision
  canonicals

  main:
  chVersions = Channel.empty()

  annotateCNV(
    variants,
    gtf
  )
  
  chVcfCnv = variants.map{it -> [it[0],it[1],it[2]]}
    .combine(annotateCNV.out.cnv, by:0)
  
  // Detect driver mutations
  oncoDriver(
    chVcfCnv,
    drivers,
    decision,
    canonicals
  )
  chVersions = chVersions.mix(oncoDriver.out.versions)

  tabix(
    oncoDriver.out.vcf
  )

  emit:
  versions = chVersions
  vcf = oncoDriver.out.vcf.join(tabix.out.tbi)
}