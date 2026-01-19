/*
 * merge stat files from mutect2 somatic variant calling
 */

process mergeMutect2Stats {
  tag "${meta.id}"
  label 'gatk'
  label 'minCpu'
  label 'medMem'

  input:
  tuple val(meta), path(stats)

  output:
  tuple val(meta), path("*.vcf.gz.stats") , emit: stats
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def statsFiles = stats.collect{"-stats ${it} " }.join(" ")
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \\
    MergeMutectStats \\
    ${statsFiles} \\
    ${args} \\
    -O ${prefix}.vcf.gz.stats

  echo "GATK "\$(gatk --version 2>&1 | grep \\(GATK\\) | sed 's/^.*(GATK) v//') > versions.txt
  """
}
