// STEP GATK MUTECT2.4 - GatherPileupSummaries

process gatherPileupSummaries {
  label 'gatk'
  label 'minCpu'
  label 'medMem'

  tag "${meta.id}"

  input:
  tuple val(meta), path(pileupSums)
  path(dict)

  output:
  tuple val(meta), path("*_pileupsummaries.table"), emit: table
  path("versions.txt"), emit: versions

  script:
  def inputs = pileupSums.collect{ "-I ${it} " }.join(' ')
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    GatherPileupSummaries \
    --sequence-dictionary ${dict} \
    ${inputs} \
    -O ${prefix}_pileupsummaries.table \
    ${args}

  echo "GATK "\$(gatk --version 2>&1 | grep \\(GATK\\) | sed 's/^.*(GATK) v//') > versions.txt
  """
}
