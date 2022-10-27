// STEP GATK MUTECT2.4 - GatherPileupSummaries

process gatherPileupSummaries {
  label 'gatk'
  label 'minCpu'
  label 'medMem'

  tag "${meta.tumor_id}_vs_${meta.normal_id}"

  input:
  tuple val(meta), path(pileupSums)
  path(dict)

  output:
  tuple val(meta), path("${meta.tumor_id}_vs_${meta.normal_id}_pileupsummaries.table.tsv"), emit: mergedPileupFileCh
  path("versions.txt"), emit: versions

  script:
  allPileups = pileupSums.collect{ "-I ${it} " }.join(' ')
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    GatherPileupSummaries \
    --sequence-dictionary ${dict} \
    ${allPileups} \
    -O ${meta.tumor_id}_vs_${meta.normal_id}_pileupsummaries.table.tsv

  echo "GATK "\$(gatk --version 2>&1 | grep \\(GATK\\) | sed 's/^.*(GATK) v//') > versions.txt
  """
}
