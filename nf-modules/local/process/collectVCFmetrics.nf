process collectVCFmetrics {
  label 'minCpu'
  label 'minMem'
  label 'unix'
  tag "${meta.id}"

  input:
  tuple val(meta), path(vcf), path(index), file(vcfFilt), path(indexFilt), file(conta)

  output:
  path("*.mqc"), emit: mqc

  script:
  def prefix = "${meta.id}"
  def filtOpt = vcfFilt.name != [] ? "-f $vcfFilt" : ""
  def contaOpt = conta.name != [] && !params.skipMutectContamination ? "-c $conta" : ""
  """
  getCallingMetrics.sh -i ${vcf} \
                        ${filtOpt} \
                        ${contaOpt} \
                       -n ${prefix} > ${prefix}_callingMetrics.mqc
  """
}
