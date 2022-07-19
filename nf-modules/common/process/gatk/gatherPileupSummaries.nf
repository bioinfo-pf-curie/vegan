// STEP GATK MUTECT2.4 - GatherPileupSummaries

process gatherPileupSummaries {
  label 'gatk'
  label 'minCpu'
  label 'medMem'

  tag "${meta.id}_${meta.tumor_id}"

  input:
  tuple val(meta), path(pileupSums)
  path(dict)

  output:
  tuple val(meta), path("${meta.tumor_id}_pileupsummaries.table.tsv"), emit: mergedPileupFileCh

  script:
  allPileups = pileupSums.collect{ "-I ${it} " }.join(' ')
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    GatherPileupSummaries \
    --sequence-dictionary ${dict} \
    ${allPileups} \
    -O ${meta.tumor_id}_pileupsummaries.table.tsv
  """
}
