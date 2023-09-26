/*
 * Identito Monitoring - combine per sample polym information
 */

process identitoCombine {
  label 'lowCpu'
  label 'lowMem'
  label 'identito'

  input:
  path(matrix)
  val(customRunName)

  output:
  path("*.tsv"), emit: tsv
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def outfilename = customRunName ? customRunName + "_identito_polym.tsv" : "identito_polym.tsv"
  """
  echo \$(R --version | awk 'NR==1{print \$1,\$3}') > versions.txt
  (head -1 "${matrix[0]}"; tail -n +2 -q *matrix.tsv) > ${outfilename}
  """
}
