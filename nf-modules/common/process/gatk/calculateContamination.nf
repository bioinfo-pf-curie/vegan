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
  path("versions.txt"), emit: versions

  script:
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    CalculateContamination \
    -I ${meta.tumor_id}_vs_${meta.normal_id}_pileupsummaries.table.tsv \
    -O ${meta.tumor_id}_vs_${meta.normal_id}_contamination.table.tsv

  echo "GATK "\$(gatk --version 2>&1 | grep \\(GATK\\) | sed 's/^.*(GATK) v//') > versions.txt
  """
}
