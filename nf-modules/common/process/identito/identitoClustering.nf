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
  path("*.{tsv,csv,png,txt}"), emit: results
  //path("versions.txt"), emit: versions

  script:
  """
  echo \$(R --version | awk 'NR==1{print \$1,\$3}') > versions.txt
  apComputeClust.R --input ${matrix} --dist ejaccard
  """
}
