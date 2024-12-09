/*
 * RNA-seq variant filtering
 */

process variantFiltration {
  tag "${meta.id}"
  label 'gatk'
  label 'highMem'
  label 'medCpu'

  input:
  tuple val(meta), path(vcf), path(tbi)
  path(fasta)
  path(fastaFai)
  path(dict)

  output:
  tuple val(meta), path("*.vcf.gz"), path("*.tbi"), emit: vcf
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def availMem = (task.memory.mega*0.8).intValue()
  """
  gatk --java-options -Xmx${availMem}M \
    VariantFiltration \
    --reference ${fasta} \
    --variant ${vcf} \
    --output ${prefix}_filtered.vcf.gz \
    --tmp-dir . \
    ${args}

  echo "GATK "\$(gatk --version 2>&1 | grep \\(GATK\\) | sed 's/^.*(GATK) v//') > versions.txt
  """
}
