process collectVCFmetrics {
  label 'minCpu'
  label 'minMem'
  label 'onlyLinux'

  input:
  tuple val(meta), val(variantCaller), path(vcf), path(tbi)

  output:
  path("*.mqc"), emit: hcCallingMetricsMqc

  """
  getCallingMetrics.sh -i ${vcf} \
                       -n ${meta.id} > ${meta.id}_${variantCaller}_callingMetrics.mqc
  """
}
