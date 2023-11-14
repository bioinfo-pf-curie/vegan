/*
 * samtools merge - Merge BAM/CRAM files
 */

process samtoolsMerge{
  tag "${meta.id}"
  label 'samtools'
  label 'highCpu'
  label 'lowMem'

  input:
  tuple val(meta), path(bams)
  path(fasta)

  output:
  tuple val(meta), path("*.bam"), optional: true,  emit: bam
  tuple val(meta), path("*.cram"), optional: true, emit: cram
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}_merged"
  def inputs = bams.collect{"${it}"}.join(' ')
  def reference = fasta ? "--reference ${fasta}" : ""
  def extension = bams instanceof List ? bams[0].getExtension() : bams.getExtension()
  """
  samtools merge \
    --threads ${task.cpus} \
    ${args} \
    ${reference} \
    ${prefix}.${extension} \
    ${inputs}
  echo \$(samtools --version | head -1) > versions.txt
  """
}
