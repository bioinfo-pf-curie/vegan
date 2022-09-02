/*
 * Germine variant calling
 */

process genotypeGVCFs {
  tag "${prefix}"
  label 'gatk'

  input:
  tuple val(meta), path(gvcf), path(index)
  path(bed)
  path(dbsnp)
  path(dbsnpIndex)
  path(fasta)
  path(fastaFai)
  path(dict)

  output:
  tuple val(meta), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"), emit: vcf
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

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
    -O ${prefix}.vcf.gz

  echo "GATK "\$(gatk --version 2>&1 | grep \\(GATK\\) | sed 's/^.*(GATK) v//') > versions.txt
  """
}
