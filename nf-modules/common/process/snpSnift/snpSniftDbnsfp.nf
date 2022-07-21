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
  val(db)

  output:
  tuple val(meta), path("*.ann.vcf.gz*"), emit: vcf
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  snpSift -Xmx${task.memory.toGiga()}g \\
    dbnsfp \\
    ${args} \\
    -db ${db} \\
    ${vcf[0]} \\
    > ${prefix}.ann.vcf

  bgzip < ${prefix}.ann.vcf > ${prefix}.ann.vcf.gz
  tabix ${prefix}.ann.vcf.gz

  echo \$(snpSift -version | cut -d" " -f1,2) > versions.txt
  """
}

