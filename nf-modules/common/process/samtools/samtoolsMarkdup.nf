/*
 * Samtools - MarkDup
 */

process samtoolsMarkdup {
  tag "${meta.id}"
  label 'samtools'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path(bam), path(bai)
  path fasta
  path fai

  output:
  tuple val(meta), path("*.bam"),  emit: bam, optional: true
  tuple val(meta), path("*.cram"), emit: cram, optional: true
  tuple val(meta), path("*.sam"),  emit: sam,  optional: true
  path "versions.txt",             emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}_markdup"
  def reference = fasta ? "--reference ${fasta}" : ""
  def extension = bam.getExtension()
  """
  samtools \\
      markdup \\
      $args \\
      ${reference} \\
      --output-fmt ${extension} \\
      -@ $task.cpus \\
      -T $prefix \\
      ${bam} \\
      ${prefix}.${extension}

  echo \$(samtools --version | head -1) > versions.txt
  """
}