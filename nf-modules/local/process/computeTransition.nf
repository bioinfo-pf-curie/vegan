/*
 * Compute Transition/Tansversion ratio
 */

process computeTransition {
  label 'minCpu'
  label 'lowMem'
  label 'transition'
  tag "${meta.id}"

  input:
  tuple val(meta), path(vcf), path(index)

  output:
  path("*table.tsv"), emit: metrics

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  parseTransition.py -i ${vcf} ${args} -o ${prefix}.transi.tsv
  transitionTable.r ${prefix}.transi.tsv ${prefix}.table.tsv
  """
}
