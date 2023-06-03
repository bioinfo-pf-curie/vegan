/*
 * Samtools - Sort
 */

process samtoolsSort {
  tag "${meta.id}"
  label 'samtools'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path (bam)

  output:
  tuple val(meta), path("*.bam"), optional: true, emit: bam
  tuple val(meta), path("*.cram"), optional: true, emit: cram
  path("versions.txt") , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${bam.baseName}_sorted"
  def extension = bam.getExtension()
  """
  echo \$(samtools --version | head -1 ) > versions.txt
  samtools sort \\
    ${args} \\
    -@  ${task.cpus}  \\
    -o ${prefix}.${extension}  \\
    ${bam}
  """
}
