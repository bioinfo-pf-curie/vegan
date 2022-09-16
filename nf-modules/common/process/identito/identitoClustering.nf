/*
 * Identito Monitoring - create the clustering
 */

process identitoClustering {
  label 'lowCpu'
  label 'lowMem'
  label 'identito'

  input:
  path(matrix)

  output:
  path("*.{tsv,csv,png}"), optional: true, emit: results
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  """
  echo \$(R --version | awk 'NR==1{print \$1,\$3}') > versions.txt
  apComputeClust.R --input ${matrix} --dist ejaccard
  """
}
