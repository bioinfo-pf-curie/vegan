/*
* ASCAT : infers tumour purity and ploidy and calculates allele-specific copy number propaths.
*/

process ascat {
  label 'ascat'
  label 'lowMem'
  label 'minCpu'

  tag "${meta.id}"

  input:
  tuple val(meta),
        path(bafNormal), path(logrNormal),
        path(bafTumor), path(logrTumor)
  path(acLociGC)

  output:
  tuple val(meta), path("*.{png,txt}"), emit: results
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ""
  """
  runAscat.r \
    --tumorbaf ${bafTumor} \
    --tumorlogr ${logrTumor} \
    --normalbaf ${bafNormal} \
    --normallogr ${logrNormal} \
    --tumorname ${meta.tumor_id} \
    --basedir ${workflow.projectDir} \
    --gcfile ${acLociGC} \
    --gender ${meta.sex} \
    ${args}

  R --slave -e "message('Facets ', packageVersion('ASCAT'))" > versions.txt
  """
}
