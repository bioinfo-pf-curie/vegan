/*
 * Process to convert table into input files for sigprofiler
 */

process preprocessTsv {
  tag "${meta.id}"
  label 'lowCpu'
  label 'medMem'
  label 'sigprofiler'

  input:
  tuple val(meta), path(tsv)

  output:
  tuple val(meta), path("*sigprofiler.vcf"), emit: table
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  echo \$(R --version | awk 'NR==1{print \$1,\$3}') > versions.txt
  vcfToFilesForSignatures.R ./ ./ ${meta.id}
  """
}