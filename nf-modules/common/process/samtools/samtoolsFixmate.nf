/*
 * Samtools - Fixmate
 */

process samtoolsFixmate {
  tag "${meta.id}"
  label 'samtools'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path(bam), path(bai)

  output:
  tuple val(meta), path ("*_fixmate.sam"), optional: true, emit: sam
  tuple val(meta), path ("*_fixmate.bam"), optional: true, emit: bam
  tuple val(meta), path ("*_fixmate.cram"), optional: true, emit: cram

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def extension = args.contains("--output-fmt sam") ? "sam" :
                  args.contains("--output-fmt bam") ? "bam" :
                  args.contains("--output-fmt cram") ? "cram" :
                  bam.getExtension()
  """
  samtools fixmate \\
    ${args} \\
    -@  ${task.cpus}  \\
    ${bam} \\
    ${prefix}_fixmate.${extension}
  """
}