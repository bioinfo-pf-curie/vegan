/*
 * Germine variant calling with HaplotypeCaller
 */

process haplotypeCaller {
  tag "${meta.id}"
  label 'gatk'
  label 'highMem'
  label 'medCpu'

  input:
  tuple val(meta), path(bam), path(bai), path(intervals)
  path(dbsnp)
  path(dbsnpIndex)
  path(fasta)
  path(fastaFai)
  path(dict)

  output:
  tuple val(meta), path("*.vcf.gz"), emit: vcf
  tuple val(meta), path("*.tbi"), optional: true, emit: tbi
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
    --reference ${fasta} \
    --input ${bam} \
    --output ${prefix}.vcf.gz \
    ${dbsnpCmd} \
    ${intervalCmd} \
    --tmp-dir . \
    ${args}

  echo "GATK "\$(gatk --version 2>&1 | grep \\(GATK\\) | sed 's/^.*(GATK) v//') > versions.txt
  """
}
