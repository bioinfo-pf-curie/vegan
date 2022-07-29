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
  """
  id=\$(bcftools query -l ${vcf} | tail -n1)
  apParseTransition.py -i ${vcf} --sample \$id -o ${prefix}.transi.tsv
  apTransition.R ${prefix}.transi.tsv ${prefix}.table.tsv
  echo "tafsta"
  """
}
