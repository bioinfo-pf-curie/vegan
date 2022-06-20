/*
 * Germine variant calling
 */

process genotypeGVCFs {
  tag "${meta.id}"
  label 'gatk'

  input:
  tuple val(meta), path(gvcf)
  path(bed)
  path(dbsnp)
  path(dbsnpIndex)
  path(fasta)
  path(fastaFai)
  path(dict)

  output:
  tuple val(meta), path("${prefix}.vcf"),emit: hcVcf

  script:
  def args = task.ext.args ?: ''
  def args2 = task.ext.args2 ?: ''
  prefix = task.ext.prefix ?: "${meta.id}"
  """
  gatk --java-options -Xmx${task.memory.toGiga()}g \
    IndexFeatureFile -I ${gvcf}

  gatk --java-options -Xmx${task.memory.toGiga()}g \
    GenotypeGVCFs \
    -R ${fasta} \
    ${args} \
    ${args2} \
    -V ${gvcf} \
    -O ${prefix}.vcf

  """
}
