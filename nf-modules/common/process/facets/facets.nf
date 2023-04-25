/*
 * Compute Cellular Fraction and Copy Numbers from Tumor Sequencing
 */

process facets{
  label 'facets'
  label 'minCpu'
  label 'highMem'
  tag "${prefix}"

  input:
  tuple val(meta), path(snppileupCounts)

  output:
  path("*.{txt,pdf}"), emit: results
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  prefix = task.ext.prefix ?: "${meta.id}"
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
