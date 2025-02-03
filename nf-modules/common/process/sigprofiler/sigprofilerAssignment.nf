/*
 * Mutational signatures
 */

process sigprofilerAssignment {
  tag "${meta.id}"
  label 'lowCpu'
  label 'medMem'
  label 'sigprofiler'

  input:
  tuple val(meta), path(vcfs)
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
  for f in *vcf.gz; do
    pfix=\$(basename \${f} .gz)
    gunzip -c \${f} > input/\${pfix}
  done

  sigProfilerAssignment.py \
    ${args} \
    --volume ./${tsb}
  """
}