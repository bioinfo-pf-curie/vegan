// STEP GATK MUTECT2.5 - CALCULATING CONTAMINATION

process calculateContamination {
  label 'gatk'
  label 'minCpu'
  label 'lowMem'

  tag "${meta.tumor_id}_vs_${meta.normal_id}"

  input:
  tuple val(meta), path(bamTumor), path(baiTumor),path(bamNormal), path(baiNormal), path(mergedPileup)

  output:
  tuple val(meta), path("${meta.tumor_id}_vs_${meta.normal_id}_contamination.table.tsv"), emit: contaminationTable

  script:
  """
  # calculate contamination
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    CalculateContamination \
    -I ${meta.tumor_id}_vs_${meta.normal_id}_pileupsummaries.table.tsv \
    -O ${meta.tumor_id}_vs_${meta.normal_id}_contamination.table.tsv
  """
}
