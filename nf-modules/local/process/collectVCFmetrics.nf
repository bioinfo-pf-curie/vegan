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
