/*
 * GATK: Filter Mutect Calls
 */

process filterMutect2Calls {
  label 'gatk'
  label 'medCpu'
  label 'medMem'

  tag "${meta.id}"

  input:
  tuple val(meta), path(vcf), path(index), path(stats), path(orientation), path(contamination), path(segmentation)
  path(dict)
  path(fasta)
  path(fastaFai)

  output:
  tuple val(meta), path("*_filtered.vcf.gz"), path("*_filtered.vcf.gz.tbi"), emit:vcf
  tuple val(meta), path("*filteringStats.tsv"), emit: stats
  path("versions.txt"), emit: versions

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def orientationCmd  = orientation  ? orientation.collect{"--orientation-bias-artifact-priors $it"}.join(' ') : ''
  def segmentationCmd = segmentation ? segmentation.collect{"--tumor-segmentation $it"}.join(' ') : ''
  def contaminationCmd = contamination ? " --contamination-table ${contamination} " : ''
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    FilterMutectCalls \
    --variant ${vcf} \
    --stats ${stats} \
    --reference ${fasta} \
    --output ${prefix}_Mutect2_filtered.vcf.gz \
    ${orientationCmd} \
    ${segmentationCmd} \
    ${contaminationCmd} \
    --tmp-dir . \
    ${args}

  awk '\$0~"^#" || \$7 == "PASS"{print}' <(bgzip -dc ${prefix}_Mutect2_filtered.vcf.gz) | bgzip > ${prefix}_Mutect2_filtered_pass.vcf.gz
  tabix ${prefix}_Mutect2_filtered_pass.vcf.gz

  echo "GATK "\$(gatk --version 2>&1 | grep \\(GATK\\) | sed 's/^.*(GATK) v//') > versions.txt
  """
}
