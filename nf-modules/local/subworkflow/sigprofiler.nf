/*
 * Worflow for mutational analysis
 */

include { sigprofilerAssignment } from '../../common/process/sigprofiler/sigprofilerAssignment'

workflow sigprofilerFlow {

  take:
  vcf
  tsb

  main:
  chVersions = Channel.empty()

  chAllVcf = vcf
  .map{ meta, vcf, tbi -> vcf }
  .collect()
  .map{vcf -> 
    def meta = [ id: "sigprofiler" ]
    [meta, vcf]
  }.view()

  sigprofilerAssignment(
    chAllVcf,
    tsb
  )

  emit:
  tsv = sigprofilerAssignment.out.sbs
  versions = chVersions
}