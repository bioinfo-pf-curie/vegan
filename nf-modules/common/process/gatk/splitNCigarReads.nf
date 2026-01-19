/*
 * GATK SplitNCigarReads
 */

process splitNCigarReads {
  tag "${meta.id}"
  label 'gatk'
  label 'minCpu'
  label 'medMem'

  input:
  tuple val(meta), path(bam), path(bai)
  path(fasta)
  path(fai)
  path(dict)

  output:
  tuple val(meta), path('*splitNCigar.bam'), emit: bam
  tuple val(meta), path('*splitNCigar.bai'), emit: bai
  path('versions.txt'), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" SplitNCigarReads \\
        --input ${bam} \\
        --output ${prefix}_splitNCigar.bam \\
        --reference ${fasta} \\
        --tmp-dir . \\
        $args

    echo "GATK "\$(gatk --version 2>&1 | grep \\(GATK\\) | sed 's/^.*(GATK) v//') > versions.txt
  """
}