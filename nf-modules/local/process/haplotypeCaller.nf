/*
 * Germine variant calling with HaplotypeCaller
 */

process haplotypeCaller {
  tag "${meta.id}"
  label 'gatk'
  label 'medMemSq'
  label 'lowCpu'

  input:
  tuple val(meta), path(bam), path(bai)
  path(bed)
  path(dbsnp)
  path(dbsnpIndex)
  path(fasta)
  path(fastaFai)
  path(dict)

  output:
  tuple val(meta), path("*.g.vcf"), emit:hcGvcf

  script:
  def args = task.ext.args ?: ''
  def args2 = task.ext.args2 ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g -Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
    HaplotypeCaller \
    -R ${fasta} \
    -I ${bam} \
    ${args} \
    ${args2} \
    -O ${prefix}.g.vcf \
    -ERC GVCF
  """
}
