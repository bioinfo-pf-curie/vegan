/*
 * Worflow for mutational analysis
 */

//include { preprocessTsv } from '../../local/process/preprocessTsv'
include { sigprofilerAssignment } from '../../common/process/sigprofiler/sigprofilerAssignment'

workflow sigprofilerFlow {

  take:
  vcf
  tsb

  main:
  chVersions = Channel.empty()

  //preprocessTsv(
  //  tsv
  //)

  sigprofilerAssignment(
    vcf,
    tsb
  )

  emit:
  tsv = sigprofilerAssignment.out.sbs
  versions = chVersions
}