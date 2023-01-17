/*
 * Samtools - View
 */

process samtoolsView {
  tag "${meta.id}"
  label 'samtools'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path (sam)

  output:
  tuple val(meta), path ("*.bam"), emit: bam
  path("*.log"), emit: logs
  path("versions.txt") , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${sam.baseName}"
  """
  echo \$(samtools --version | head -1 ) > versions.txt
  samtools view -bS\\
    -@  ${task.cpus}  \\
    -o ${prefix}.bam  \\
    ${sam}

  getBWAstats.sh -i ${prefix}.bam -p ${task.cpus} > ${prefix}_bwa.log
  """
}
