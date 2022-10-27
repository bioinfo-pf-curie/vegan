/*
 * VCF annotation with snpSnift
 */

process snpSiftExtractFields {
  label 'snpsift'
  label 'lowMem'
  label 'lowCpu'

  tag "${meta.id}"

  input:
  tuple val(meta), path(vcf)

  output:
  tuple val(meta), path("*.tsv"), emit: tsv
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def args2 = task.ext.args2 ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  SnpSift -Xmx${task.memory.toGiga()}g \
    extractFields \
    ${args} \
    ${vcf[0]} \
    ${args2} \
    > ${prefix}.tsv

  echo "snpSift "\$(SnpSift 2>&1 | awk '\$0~"SnpSift version"{print \$3}') > versions.txt
  echo "tot"
  """
}
