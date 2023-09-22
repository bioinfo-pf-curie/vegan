/*
 * Msisensorpro Pro - evaluate MSI using single (tumor) sample sequencing data
 */

process msisensorproPro {
  label 'minCpu'
  label 'lowMem'
  label 'msisensorpro'

  input:
  tuple val(meta), path(bam), path(bai)
  path(baseline)
  path(bed)

  output:
  tuple val(meta), path("${meta.id}_msi.txt") , emit: outputReport
  tuple val(meta), path("${meta.id}_dis")     , emit: outputDis
  tuple val(meta), path("${meta.id}_all")     , emit: outputAll
  tuple val(meta), path("${meta.id}_unstable"), emit: outputUnstable
  path "versions.txt" , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix   = task.ext.prefix ?: "${meta.id}"
  def bed = bed ? " -e ${bed} " : ""
  """
  msisensor-pro \\
      pro \\
      ${args} \\
      ${bed} \\
      -d ${baseline} \\
      -t ${bam} \\
      -o ${prefix}_msi.txt

  echo "msisensor-pro "\$(msisensor-pro 2>&1 | sed -nE 's/Version:\\sv([0-9]\\.[0-9])/\\1/ p') > versions.txt
  """
}
