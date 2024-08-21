/*
 * Apply decision algorithm to select drivers mutation
 */

include { oncoDrivers } from '../../local/process/oncoDrivers'
include { annotateCNV } from '../../local/process/annotateCNV'

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
  oncoDrivers(
    chVcfCnv,
    drivers,
    decision,
    canonicals
  )
  chVersions = chVersions.mix(oncoDrivers.out.versions)

  emit:
  versions = chVersions
  vcf = oncoDrivers.out.vcf
}