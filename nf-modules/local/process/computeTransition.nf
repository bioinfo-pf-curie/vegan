/*
 * Compute Transition/Tansversion ratio
 */

process computeTransition {
  label 'minCpu'
  label 'lowMem'
  label 'transition'
  tag "${prefix}"

  input:
  tuple val(meta), path(vcf), path(index)

  output:
  path("*table.tsv"), emit: metrics

  when:
  task.ext.when == null || task.ext.when

  script:
  prefix = task.ext.prefix ?: "${meta.id}"
  """
  parseTransition.py -i ${vcf} --sample ${prefix} -o ${prefix}_${meta.status}.transi.tsv
  transitionTable.r ${prefix}_${meta.status}.transi.tsv ${prefix}_${meta.status}.table.tsv
  """
}
