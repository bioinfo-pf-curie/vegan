/*
* STEP ASCAT.2 - CONVERTALLELECOUNTS
*/

process convertAlleleCounts {
  label 'ascat'
  label 'lowMem'
  label 'minCpu'

  tag "${meta.tumor_id}_vs_${meta.normal_id}"

  input:
  tuple val(meta), path(alleleCountTumor), path(alleleCountNormal)

  output:
  tuple val(meta),
        path("${meta.normal_id}.BAF"), path("${meta.normal_id}.LogR"),
        path("${meta.tumor_id}.BAF"), path("${meta.tumor_id}.LogR"),emit: counts
  path("versions.txt"),emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  """
  Rscript ${workflow.projectDir}/bin/apConvertAlleleCounts.r ${meta.tumor_id} ${alleleCountTumor} ${meta.normal_id} ${alleleCountNormal} ${meta.sex}
  R -e "packageVersion('ASCAT')" > versions.txt
  """
}
