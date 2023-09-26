/*
 * Msisensorpro Baseline - prepare a config file from a list of BAMs
 */

process prepareBaselineConfig {
  label 'minCpu'
  label 'minMem'
  label 'msisensorpro'

  input:
  path(bams) 
  path(bai)

  output:
  path("msi_baseline.txt"), emit: config

  when:
  task.ext.when == null || task.ext.when

  script:
  """
  for f in ${bams}
  do
    echo -e "\${f}\t\${PWD}/\${f}" >> msi_baseline.txt 
  done
  """
}
