/*
 * GATK: GENERATING PILEUP SUMMARIES
 */

process getPileupSummaries {
  label 'gatk'
  label 'minCpu'
  label 'extraMem'

  tag "${meta.tumor_id}_vs_${meta.normal_id}"

  input:
  tuple val(meta), path(bamTumor),path(baiTumor), path(bamNormal), path(baiNormal)
  path(intervalBed)
  path(pileupSum)
  path(pileupSumIndex)
  path(targetBed)

  output:
  tuple val(meta), path("*_pileupsummaries.table"), emit: pileupSummaries
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  //pairName = pairMap[[sampleIdNormal, sampleIdTumor]]
  //intervalOpts = params.noIntervals ? params.targetBed ? "-L ${targetBed}" : "-L ${pileupSum}" : "-L ${intervalBed}"

  def args = task.ext.args ?: ''
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    GetPileupSummaries \
    -I ${bamTumor} \
    -V ${pileupSum} \
    ${args} \
    -O ${meta.tumor_id}_vs_${meta.normal_id}_pileupsummaries.table
  echo "GATK "\$(gatk --version 2>&1 | grep \\(GATK\\) | sed 's/^.*(GATK) v//') > versions.txt
  """
}
