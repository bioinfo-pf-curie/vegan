/*
 * Msisensor-pro msi -  evaluate MSI using paired tumor-normal sequencing data
 */

process msisensorproMsi {
  tag "$meta.id"
  label 'minCup'
  label 'lowMem'
  label 'msisensorpro'

  input:
  tuple val(meta), path(tumor), path(tumorIndex), path(normal), path(normalIndex)
  path (fasta)
  path (msisensorScan)

  output:
  tuple val(meta), path("${prefix}")         , emit: outputReport
  tuple val(meta), path("${prefix}_dis")     , emit: outputDis
  tuple val(meta), path("${prefix}_germline"), emit: outputGermline
  tuple val(meta), path("${prefix}_somatic") , emit: outputSomatic
  path "versions.txt"                        , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args   ?: ''
  def prefix   = task.ext.prefix ?: "${meta.id}"
  def fasta = fasta ? "-g ${fasta}" : ""
  //def intervals = intervals ? " -e ${intervals} " : ""
  """
  msisensor-pro \\
      msi \\
      -d ${msisensorScan} \\
      -n ${normal} \\
      -t ${tumor} \\
      ${fasta} \\
      -o $prefix \\
      -b ${task.cpus} \\
      $args

  echo "missensor-pro "\$(msisensor-pro 2>&1 | sed -nE 's/Version:\\sv([0-9]\\.[0-9])/\\1/ p') > versions.txt
  """
}
