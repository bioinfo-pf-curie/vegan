/*
 * Germine variant calling
 */

process genotypeGVCFs {
  tag "${meta.id}"
  label 'gatk'

  input:
  tuple val(meta), path(gvcf), path(index), path(intervals)
  path(dbsnp)
  path(dbsnpIndex)
  path(fasta)
  path(fastaFai)
  path(dict)

  output:
  tuple val(meta), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcf
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def intervalCmd = intervals ? "-L ${intervals}" : ""
  def dbsnpCmd = dbsnp ? "--D ${dbsnp}" : ""
  """
  gatk --java-options -Xmx${task.memory.toGiga()}g \
    GenotypeGVCFs \
    -R ${fasta} \
    ${args} \
    ${intervalCmd} \
    ${dbsnpCmd} \
    -V ${gvcf} \
    -O ${prefix}.vcf.gz

  echo "GATK "\$(gatk --version 2>&1 | grep \\(GATK\\) | sed 's/^.*(GATK) v//') > versions.txt
  """
}
