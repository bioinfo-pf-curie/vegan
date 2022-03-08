/*
 * samtools merge - Merge BAM files
 */

process samtoolsMerge{
  label 'samtools'
  label 'highCpu'
  label 'lowMem'

  tag "${sampleId}"

  input:
  tuple val(prefix), path(bams)

  output:
  tuple val(prefix), file("*_merged.bam"), emit: bam
  path("versions.txt"), emit: versions

  script:
  def args = task.ext.args ?: ''
  inputs = bams.collect{"${it}"}.join(' ')
  """
  samtools merge --threads ${task.cpus} ${args} ${prefix}_merged.bam ${inputs}
  echo \$(samtools --version | head -1) > versions.txt
  """
}

