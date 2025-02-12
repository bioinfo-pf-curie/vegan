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
  def prefix = task.ext.prefix ? "--prefix ${task.ext.prefix}" : ""
  """
  echo \$(R --version | awk 'NR==1{print \$1,\$3}') > versions.txt
  computeClust.R --input ${matrix} --dist ejaccard ${prefix}
  """
}
