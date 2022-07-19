/*
 * GATK: GENERATING PILEUP SUMMARIES
 */

process getPileupSummaries {
  label 'gatk'
  label 'minCpu'
  label 'extraMem'

  tag "${meta.tumor_id}_vs_${meta.normal_id}_${intervalBed.baseName}"

  input:
  tuple val(meta), path(bamTumor),path(baiTumor), path(bamNormal), path(baiNormal)
  path(intervalBed)
  path(germlineResource)
  path(germlineResourceIndex)
  path(targetBed)

  output:
  tuple val(meta), path("*_pileupsummaries.table"), emit: pileupSummaries

  when:
  task.ext.when == null || task.ext.when

  script:
  //pairName = pairMap[[sampleIdNormal, sampleIdTumor]]
  //intervalOpts = params.noIntervals ? params.targetBed ? "-L ${targetBed}" : "-L ${germlineResource}" : "-L ${intervalBed}"

  def args = task.ext.args ?: ''
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    GetPileupSummaries \
    -I ${bamTumor} \
    -V ${germlineResource} \
    ${args} \
    -O ${intervalBed.baseName}_${meta.tumor_id}_pileupsummaries.table
  """
}
