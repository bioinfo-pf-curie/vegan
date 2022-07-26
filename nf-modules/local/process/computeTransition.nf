/*
 * Compute Transition/Tansversion ratio
 */

process computeTransition {
  label 'minCpu'
  label 'lowMem'
  label 'transition'
  tag "${fileID}"

  input:
  tuple val(meta), path(vcf), path(index)

  output:
  path("*table.tsv"), emit: metrics

  when:
  task.ext.when == null || task.ext.when

  script:
  //def prefix = task.ext.prefix ?: "${meta.id}"
  fileID = "${meta.status}" == "pair" ? "${meta.tumor_id}_vs_${meta.normal_id}" : "${meta.status}" == "tumor" ? "${meta.tumor_id}" : "${meta.normal_id}"
  """
  id=\$(bcftools query -l ${vcf} | tail -n1)
  apParseTransition.py -i ${vcf} --sample \$id -o ${fileID}.transi.tsv
  apTransition.R ${fileID}.transi.tsv ${fileID}.table.tsv
  """
}
