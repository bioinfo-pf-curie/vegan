/*
 * GATK: GENERATING PILEUP SUMMARIES
 */

process getPileupSummaries {
  label 'gatk'
  label 'minCpu'
  label 'extraMem'
  tag "${meta.id}"

  input:
  tuple val(meta), path(bam), path(bai), path(intervals)
  path(pileupSum)
  path(pileupSumIndex)
  path(fasta)
  path(fai)
  path(dict)

  output:
  tuple val(meta), path("*_pileups.table"), emit: table
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def intervalCmd = intervals ? "-L ${intervals}" : "-L ${pileupSum}"
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g" GetPileupSummaries \\
    --input ${bam} \\
    --variant ${pileupSum} \\
    --output ${prefix}_pileups.table \\
    --reference ${fasta} \\
    --sequence-dictionary ${dict} \\
    --tmp-dir . \\
    ${intervalCmd} \\
    ${args}

  echo "GATK "\$(gatk --version 2>&1 | grep \\(GATK\\) | sed 's/^.*(GATK) v//') > versions.txt
  """
}
