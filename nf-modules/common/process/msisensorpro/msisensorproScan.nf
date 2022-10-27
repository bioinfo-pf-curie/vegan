/*
 * Msisensorpro Scan - scan the reference genome to get microsatellites information
 */

process msisensorproScan {
  label 'minCpu'
  label 'lowMem'
  label 'msisensorpro'

  input:
  path(fasta)

  output:
  path("*.list"), emit: list
  path "versions.txt" , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  """
  msisensor-pro \\
      scan \\
      -d $fasta \\
      -o msisensor_scan.list \\
      $args

  echo "msisensor-pro "\$(msisensor-pro 2>&1 | sed -nE 's/Version:\\sv([0-9]\\.[0-9])/\\1/ p') > versions.txt
  """
}
