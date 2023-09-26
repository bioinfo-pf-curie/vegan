process collectVCFmetrics {
  label 'minCpu'
  label 'minMem'
  label 'unix'
  tag "${meta.id}"

  input:
  tuple val(meta), path(vcf), path(index), file(vcfFilt), path(indexFilt)

  output:
  tuple val(meta), path("*.stats"), emit: stats

  script:
  def prefix = "${meta.id}"
  """
  getCallingMetrics.sh -i ${vcf} \
                       -f ${vcfFilt} \
                       -n ${prefix} > ${prefix}_callingMetrics.stats
  """
}
