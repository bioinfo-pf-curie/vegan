/*
 * Samtools - Index BAM/CRAM file
 */

process samtoolsIndex {
  tag "${meta.id}"
  label 'samtools'
  label 'minCpu'
  label 'lowMem'
 
  input:
  tuple val(meta), path (bam)

  output:
  tuple val(meta), path("*bam.bai")  , optional: true, emit: bai
  tuple val(meta), path("*cram.crai"), optional: true, emit: crai
  path("versions.txt")               , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  """
  echo \$(samtools --version | head -1) > versions.txt
  samtools index ${bam}
  """
}
    