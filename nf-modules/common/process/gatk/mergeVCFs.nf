/*
 * GATK merge vcfs file
 */

process mergeVCFs {
  tag "$meta.id"
  label 'gatk'
  label 'minCpu'
  label 'medMem'

  input:
  tuple val(meta), path(vcf), path(tbi)
  path  dict

  output:
  tuple val(meta), path('*.vcf.gz'), path("*.tbi"), emit: vcf
  path  "versions.txt"             , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def inputList = vcf.collect{ "--INPUT $it"}.join(' ')
  def referenceCmd = dict ? "--SEQUENCE_DICTIONARY $dict" : ""

  """
  gatk --java-options -Xmx${task.memory.toGiga()}g MergeVcfs \\
      $inputList \\
      --OUTPUT ${prefix}.vcf.gz \\
      $referenceCmd \\
      --TMP_DIR . \\
      $args
  echo "GATK "\$(gatk --version 2>&1 | grep \\(GATK\\) | sed 's/^.*(GATK) v//') > versions.txt
  """
}
