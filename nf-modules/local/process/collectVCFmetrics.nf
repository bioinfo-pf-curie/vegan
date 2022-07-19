process collectVCFmetrics {
  label 'minCpu'
  label 'minMem'
  label 'onlyLinux'

  input:
  tuple val(meta), path(vcf)

  output:
  path("*.mqc"), emit: mqc

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  getCallingMetrics.sh -i ${vcf[0]} \
                       -n ${meta.id} > ${prefix}_callingMetrics.mqc
  """
}
