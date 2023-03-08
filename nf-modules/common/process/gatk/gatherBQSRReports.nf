/*
 * BQSR - gather BQSR reports
 */

process gatherBQSRReports {
  tag "${meta.id}"
  label 'gatk'
  label 'minCpu'
  label 'medMem'

  input:
  tuple val(meta), path(table)


  output:
  tuple val(meta), path("*.recal.table"), emit: table
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  input = table.collect{"-I ${it}"}.join(' ')
  """
  gatk --java-options -Xmx${task.memory.toGiga()}g \
      GatherBQSRReports \
      ${input} \
      -O ${prefix}.recal.table
  echo "GATK "\$(gatk --version 2>&1 | grep \\(GATK\\) | sed 's/^.*(GATK) v//') > versions.txt
  """
}
