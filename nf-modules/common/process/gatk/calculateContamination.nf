// STEP GATK MUTECT2.5 - CALCULATING CONTAMINATION

process calculateContamination {
  label 'gatk'
  label 'minCpu'
  label 'lowMem'

  tag "${meta.id}"

  input:
  tuple val(meta), path(pileup), path(matched)

  output:
  tuple val(meta), path("*contamination.table.tsv"), emit: contamination
  tuple val(meta), path("*segmentation.table.tsv"), emit: segmentation
  path("versions.txt"), emit: versions

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def matchedCmd = matched ? "--matched-normal $matched" : ''
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    CalculateContamination \
    --input ${pileup} \
    --output ${prefix}_contamination.table.tsv \
    --tumor-segmentation ${prefix}_segmentation.table.tsv \
    --tmp-dir . \
    ${matchedCmd} \
    ${args}

  echo "GATK "\$(gatk --version 2>&1 | grep \\(GATK\\) | sed 's/^.*(GATK) v//') > versions.txt
  """
}
