/*
 * VCF annotation with snpSnift
 */

process snpSiftAnnotate {
  label 'snpsift'
  label 'lowMem'
  label 'lowCpu'

  tag "${meta.id}"

  input:
  tuple val(meta), path(vcf)
  tuple path(db), path(tbi)

  output:
  tuple val(meta), path("*.vcf{.gz,.gz.tbi}"), emit: vcf
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}.ann"
  """
  SnpSift -Xmx${task.memory.toGiga()}g \\
    annotate \\
    ${args} \\
    ${db[0]} \\
    ${vcf[0]} \\
    > ${prefix}.vcf

  bgzip < ${prefix}.vcf > ${prefix}.vcf.gz
  tabix ${prefix}.vcf.gz

  echo "snpSift "\$(SnpSift 2>&1 | awk '\$0~"SnpSift version"{print \$3}') > versions.txt
  echo "tabix "\$(tabix 2>&1 | awk '\$1~"Version"{print \$2}') >> versions.txt
  """
}

