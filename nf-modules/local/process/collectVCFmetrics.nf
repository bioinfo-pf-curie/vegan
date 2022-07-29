process collectVCFmetrics {
  label 'minCpu'
  label 'minMem'
  label 'onlyLinux'

  tag "${meta.id}"

  input:
  tuple val(meta), path(vcf), path(index), file(vcfFilt), path(indexFilt), file(conta)

  output:
  path("*.mqc"), emit: mqc

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def filtOpt = vcfFilt.name != [] ? "-f $vcfFilt" : ""
  def contaOpt = conta.name != [] && !params.skipMutectContamination ? "-c $conta" : ""
  """
  getCallingMetrics.sh -i ${vcf} \
                        $filtOpt \
                        $contaOpt \
                       -n ${prefix} > ${prefix}_callingMetrics.mqc
  """
}


// process collectVCFmetrics {
//   label 'minCpu'
//   label 'minMem'
//   label 'onlyLinux'
//   tag "${fileID}"
//
//   input:
//   tuple val(meta), path(vcf)
//   tuple val(meta), path(unfilteredVcf)
//
//   output:
//   path("*.mqc"), emit: mqc
//
//   script:
//   //def prefix = task.ext.prefix ?: "${meta.id}"
//   def args = task.ext.args ?: ''
//   fileID = "${meta.status}" == "pair" ? "${meta.tumor_id}_vs_${meta.normal_id}" : "${meta.status}" == "tumor" ? "${meta.tumor_id}" : "${meta.normal_id}"
//   """
//   getCallingMetrics.sh -i ${unfilteredVcf} \
//                        -f ${vcf[0]} \
//                        ${args} \
//                        -n ${fileID} > ${fileID}_callingMetrics.mqc
//   """
// }
