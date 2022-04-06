/*
 * BAM Filtering based on samtools
 */

process samtoolsFilter {
  tag "${meta.id}"
  label 'samtools'
  label 'lowCpu'
  label 'lowMem'

  input:
  tuple val(meta), path(bam)

  output:
  tuple val(meta), path("*filtered.bam"), emit: bam
  path("versions.txt"), emit: versions

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  echo \$(samtools --version | head -1) > versions.txt
  samtools view \\
    ${args} \\
    -b ${bam} > ${prefix}_filtered.bam
  """
}
