/*
 * Worflow for mutational analysis
 */

include { preprocessTsv } from '../../local/process/preprocessTsv'
include { sigprofilerAssignment } from '../../common/process/sigprofiler/sigprofilerAssignment'

workflow sigprofilerFlow {

  take:
  tsv //basic tableReport
  tsb

  main:
  chVersions = Channel.empty()

  preprocessTsv(
    tsv
  )

  sigprofilerAssignment(
    preprocessTsv.out.table,
    tsb
  )

  emit:
  tsv = sigprofilerAssignment.out.SBS
  versions = chVersions
}