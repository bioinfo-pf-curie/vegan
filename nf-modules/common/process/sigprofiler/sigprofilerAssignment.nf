/*
 * Mutational signatures
 */

process sigprofilerAssignment {
  tag "${meta.id}"
  label 'lowCpu'
  label 'medMem'
  label 'sigprofiler'

  input:
  tuple val(meta), path(vcf)
  path(tsb)

  output:
  tuple val(meta), path("Assignment_Solution/"), emit: sbs
  path "versions.txt", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir input
  gunzip -c ${vcf} > input/${prefix}.vcf
  sigProfilerAssignment.py ${args} --volume ./${tsb}
  """
}