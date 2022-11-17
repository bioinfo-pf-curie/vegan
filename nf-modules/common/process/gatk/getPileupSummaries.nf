/*
 * GATK: GENERATING PILEUP SUMMARIES
 */

process getPileupSummaries {
  label 'gatk'
  label 'minCpu'
  label 'extraMem'
  tag "${meta.id}"

  input:
  tuple val(meta), path(bam), path(bai)
  path(intervalBed)
  path(pileupSum)
  path(pileupSumIndex)
  path(targetBed)

  output:
  tuple val(meta), path("*_pileups.table"), emit: table
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  //pairName = pairMap[[sampleIdNormal, sampleIdTumor]]
  def intervalCmd= params.targetBed ? "-L ${targetBed}" : "-L ${pileupSum}"
  //params.noIntervals ? params.targetBed ? "-L ${targetBed}" : "-L ${pileupSum}" : "-L ${intervalBed}"
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    GetPileupSummaries \
    --input ${bam} \
    --variant ${pileupSum} \
    --output ${prefix}_pileups.table \
    --tmp-dir . \
    ${intervalCmd} \
    ${args}

  echo "GATK "\$(gatk --version 2>&1 | grep \\(GATK\\) | sed 's/^.*(GATK) v//') > versions.txt
  """
}
