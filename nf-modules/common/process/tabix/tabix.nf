/*
 * Compress and index vcf file with tabix
 */

process tabix {
  label 'tabix'
  label 'lowMem'
  label 'lowCpu'

  tag "${meta.id}"

  input:
  tuple val(meta), path(vcf)

  output:
  tuple val(meta), path("*.vcf.gz*"), emit: vcf
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  prefix = task.ext.prefix ?: "${meta.id}"
  compress = vcf.getExtension() == "vcf" ? "bgzip < ${vcf} > ${vcf}.gz" : ""
  """
  $compress
  tabix ${vcf}.gz
  echo "tabix "\$(tabix 2>&1 | awk '\$1~"Version"{print \$2}') > versions.txt
  """
}
