/*
 * getFragmentSize:
 * External parameters :
 */

process getFragmentSize {
  tag "${meta.id}"
  label 'picard'
  label 'minCpu'
  label 'medMem'

  input:
  tuple val(meta), path(bamFiltered), path(baiFiltered)

  output:
  path("*_insert_size_{hist.pdf,metrics.txt}"), emit: metrics
  path("versions.txt"), emit: versions

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  echo \$(picard CollectInsertSizeMetrics --version 2>&1) &> versions.txt
  picard CollectInsertSizeMetrics \
      I=${bamFiltered} \
      O=${prefix}_insert_size_metrics.txt \
      H=${prefix}_insert_size_hist.pdf \
      M=0.5
  """
}
