/*
 * Germine variant calling with HaplotypeCaller
 */

process haplotypeCaller {
  tag "${meta.id}"
  label 'gatk'
  label 'medMemSq'
  label 'lowCpu'

  input:
  tuple val(meta), path(bam), path(bai), path(intervals)
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
  def prefix = task.ext.prefix ?: "${meta.id}"
  def intervalCmd = intervals ? "-L ${intervals}" : ""
  def dbsnpCmd = dbsnp ? "--D ${dbsnp}" : ""
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g -Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
    HaplotypeCaller \
    -R ${fasta} \
    -I ${bam} \
    ${args} \
    ${dbsnpCmd} \
    ${intervalCmd} \
    -O ${prefix}.g.vcf \
    -ERC GVCF

  echo "GATK "\$(gatk --version 2>&1 | grep \\(GATK\\) | sed 's/^.*(GATK) v//') > versions.txt
  """
}