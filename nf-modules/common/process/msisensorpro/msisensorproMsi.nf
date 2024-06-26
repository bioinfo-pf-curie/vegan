/*
 * Msisensor-pro msi -  evaluate MSI using paired tumor-normal sequencing data
 */

process msisensorproMsi {
  tag "$meta.id"
  label 'minCpu'
  label 'lowMem'
  label 'msisensorpro'

  input:
  tuple val(meta), path(tumor), path(tumorIndex), path(normal), path(normalIndex)
  path(fasta)
  path(msisensorScan)
  path(bed)

  output:
  tuple val(meta), path("*msi.txt") , emit: outputReport
  tuple val(meta), path("*dis")              , emit: outputDis
  tuple val(meta), path("*germline")         , emit: outputGermline
  tuple val(meta), path("*somatic")          , emit: outputSomatic
  path "versions.txt"                        , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args   ?: ''
  def prefix   = task.ext.prefix ?: "${meta.id}"
  def fasta = fasta ? "-g ${fasta}" : ""
  def bed = bed ? " -e ${bed} " : ""
  """
  msisensor-pro \\
      msi \\
      -d ${msisensorScan} \\
      -n ${normal} \\
      -t ${tumor} \\
      ${fasta} \\
      ${bed} \\
      -o ${prefix}_msi.txt \\
      -b ${task.cpus} \\
      $args

  echo "missensor-pro "\$(msisensor-pro 2>&1 | sed -nE 's/Version:\\sv([0-9]\\.[0-9])/\\1/ p') > versions.txt
  """
}
