/*
 * Apply decision algorithm to select drivers mutation
 */

include { oncoDrivers } from '../../local/process/oncoDrivers'
include { annotateCNV } from '../../local/process/annotateCNV'

workflow driverSomaticFlow {

  take:
  vcf // [meta, vcf, tbi]
  cnv
  gtf
  drivers
  decision
  canonicals

  main:
  chVersions = Channel.empty()

  annotateCNV(
    cnv,
    gtf
  )
  
  // Detect driver mutations
  oncoDrivers(
    vcf,
    annotateCNV.out.cnv,
    drivers,
    decision,
    canonicals
  )
  chVersions = chVersions.mix(oncoDrivers.out.versions)

  emit:
  versions = chVersions
  vcf = oncoDrivers.out.vcf
}