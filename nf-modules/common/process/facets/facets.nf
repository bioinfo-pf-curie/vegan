/*
 * Compute Cellular Fraction and Copy Numbers from Tumor Sequencing
 */

process facets{
  label 'facets'
  label 'minCpu'
  label 'highMem'

  tag "${prefix}"

  when:
  task.ext.when == null || task.ext.when

  input:
  tuple val(meta), path(snppileupCounts)

  output:
  path("*.{txt,pdf}"), emit: facetsResultsCh
  path("versions.txt"), emit: versions

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  genome=\$(basename ${params.genome} "_base")

  facets.r -i ${snppileupCounts} \\
	   --name ${prefix} \\
	   --assembly \$genome \\
           ${args}
  R -e "packageVersion('facets')" > versions.txt
  """
 }
