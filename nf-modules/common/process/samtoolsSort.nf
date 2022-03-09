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
  tuple val(meta), path ("*_sorted.bam"), emit: bam
  path("versions.txt") , emit: versions

  script:
  def args = task.ext.args ?: ''
  """
  echo \$(samtools --version | head -1 ) > versions.txt
  samtools sort \\
    ${args} \\
    -@  ${task.cpus}  \\
    -o ${bam.baseName}_sorted.bam  \\
    ${bam}
  """
}
    