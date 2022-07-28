/*
 * VCF annotation with dbNSFP
 */

process snpSiftDbnsfp {
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
  prefix = task.ext.prefix ?: "${meta.id}.ann"
  """
  SnpSift -Xmx${task.memory.toGiga()}g \\
    dbnsfp \\
    ${args} \\
    -db ${db[0]} \\
    ${vcf[0]} \\
    > ${prefix}.vcf

  bgzip < ${prefix}.vcf > ${prefix}.vcf.gz
  tabix ${prefix}.vcf.gz

  echo "snpSift "\$(SnpSift 2>&1 | awk '\$0~"SnpSift version"{print \$3}') > versions.txt
  """
}
