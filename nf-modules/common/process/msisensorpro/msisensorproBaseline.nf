/*
 * Msisensorpro Baseline - build baseline for tumor only detection
 */

process msisensorproBaseline {
  label 'minCpu'
  label 'lowMem'
  label 'msisensorpro'

  input:
  path(list)
  path(config)

  output:
  path("msisensor_baseline/*list_baseline"), emit: baseline
  path "versions.txt" , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  """
  msisensor-pro \\
      baseline \\
      ${args} \\
      -d ${list} \\
      -i ${config} \\
      -o  msisensor_baseline

  echo "msisensor-pro "\$(msisensor-pro 2>&1 | sed -nE 's/Version:\\sv([0-9]\\.[0-9])/\\1/ p') > versions.txt
  """
}
