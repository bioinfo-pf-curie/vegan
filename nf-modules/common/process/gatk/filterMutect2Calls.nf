/*
 * GATK: Filter Mutect Calls
 */

process filterMutect2Calls {
  label 'gatk'
  label 'medCpu'
  label 'medMem'

  tag "${meta.tumor_id}_vs_${meta.normal_id}"

  input:
  tuple val(meta), path(unfiltered), path(unfilteredIndex),
        path(statspath),
        path(contaminationTable)
  path(dict)
  path(fasta)
  path(fastaFai)
  path(germlineResource)
  path(germlineResourceIndex)
  path(intervals)

  output:
  tuple val(meta), val("Mutect2"), path("*_filtered_pass.vcf.gz"), path("*_filtered_pass.vcf.gz.tbi"), path("*filteringStats.tsv"), emit: filteredMutect2OutputCh
  path("*.mqc") ,emit: mutect2CallingMetricsMqcCh

  script:
  def args = task.ext.args ?: ''
  def args2 = task.ext.args2 ?: ''
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    FilterMutectCalls \
    -V ${unfiltered} \
    ${args} \
    --stats ${meta.tumor_id}_vs_${meta.normal_id}.vcf.gz.stats \
    -R ${fasta} \
    -O ${meta.tumor_id}_vs_${meta.normal_id}_Mutect2_filtered.vcf.gz

  awk '\$0~"^#" || \$7 == "PASS"{print}' <(bgzip -dc ${meta.tumor_id}_vs_${meta.normal_id}_Mutect2_filtered.vcf.gz) | bgzip > ${meta.tumor_id}_vs_${meta.normal_id}_Mutect2_filtered_pass.vcf.gz
  tabix ${meta.tumor_id}_vs_${meta.normal_id}_Mutect2_filtered_pass.vcf.gz

  getCallingMetrics.sh -i ${unfiltered} \
                       -f ${meta.tumor_id}_vs_${meta.normal_id}_Mutect2_filtered_pass.vcf.gz \
                       ${args2} \
                       -n ${meta.tumor_id}_vs_${meta.normal_id} > ${meta.tumor_id}_vs_${meta.normal_id}_Mutect2_callingMetrics.mqc
  """
}
