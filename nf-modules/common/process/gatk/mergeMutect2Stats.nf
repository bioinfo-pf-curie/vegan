/*
 * merge stat files from mutect2 somatic variant calling
 */

process mergeMutect2Stats {
  tag "${meta.tumor_id}_vs_${meta.normal_id}"
  label 'gatk'
  label 'minCpu'
  label 'medMem'

  input:
  tuple val(meta), path(statsFiles)

  output:
  tuple val(meta), path("${meta.tumor_id}_vs_${meta.normal_id}.vcf.gz.stats") , emit: mergedStatsFile
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def stats = statsFiles.collect{"-stats ${it} " }.join(" ")
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \
    MergeMutectStats \
    ${stats} \
    -O ${meta.tumor_id}_vs_${meta.normal_id}.vcf.gz.stats

  echo "GATK "\$(gatk --version 2>&1 | grep \\(GATK\\) | sed 's/^.*(GATK) v//') > versions.txt
  """
}
