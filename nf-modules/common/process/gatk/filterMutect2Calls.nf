/*
 * GATK: Filter Mutect Calls
 */

process filterMutect2Calls {
  label 'gatk'
  label 'medCpu'
  label 'medMem'

  tag "${meta.tumor_id}_vs_${meta.normal_id}"

  input:
  tuple val(meta), path(vcf), path(stats), path(contaminationTable)
  path(dict)
  path(fasta)
  path(fastaFai)
  path(germlineResource)
  path(germlineResourceIndex)
  path(intervals)

  output:
  tuple val(meta), path("*_filtered*.vcf.gz*"), emit: vcf
  tuple val(meta), path("*filteringStats.tsv"), emit: stats
  path("versions.txt"), emit: versions

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.tumor_id}_vs_${meta.normal_id}"
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    FilterMutectCalls \
    -V ${vcf[0]} \
    ${args} \
    --stats ${stats} \
    -R ${fasta} \
    -O ${prefix}_Mutect2_filtered.vcf.gz

  awk '\$0~"^#" || \$7 == "PASS"{print}' <(bgzip -dc ${meta.tumor_id}_vs_${meta.normal_id}_Mutect2_filtered.vcf.gz) | bgzip > ${meta.tumor_id}_vs_${meta.normal_id}_Mutect2_filtered_pass.vcf.gz
  tabix ${meta.tumor_id}_vs_${meta.normal_id}_Mutect2_filtered_pass.vcf.gz

  gatk FilterMutectCalls --version 2>&1 | sed 's/^.*(GATK) v/GATK /; s/Version: //; s/ .*\$//' | tail -3 > versions.txt
  echo "tot"
  """
}
