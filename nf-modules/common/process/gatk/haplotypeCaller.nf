/*
 * Germine variant calling with HaplotypeCaller
 */

process haplotypeCaller {
  tag "${prefix}"
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
  tuple val(meta), path("*.g.vcf"), path("*.g.vcf.idx"), emit: gvcf
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def args2 = task.ext.args2 ?: ''
  prefix = task.ext.prefix ?: "${meta.id}"

  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g -Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
    HaplotypeCaller \
    -R ${fasta} \
    -I ${bam} \
    ${args} \
    ${args2} \
    -O ${prefix}.g.vcf \
    -ERC GVCF

  echo "GATK "\$(gatk --version 2>&1 | grep \\(GATK\\) | sed 's/^.*(GATK) v//') > versions.txt
  """
}
