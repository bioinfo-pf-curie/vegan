/*
 * Samtools - Stats
 */

process samtoolsStats {
  tag "${meta.id}"
  label 'samtools'
  label 'minCpu'
  label 'lowMem'

  input:
  tuple val(meta), path(bam)
  path(bed)

  output:
  tuple val(meta), path("*stats"), emit: stats
  path("versions.txt") , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  targetOpts = bed ? "-t ${bed}" : ""
  """
  echo \$(samtools --version | head -1) > versions.txt
  samtools stats ${args} ${targetOpts} ${bam} > ${prefix}.stats
  """
}
