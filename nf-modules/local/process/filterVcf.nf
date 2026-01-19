process filterVcf {
  label 'minCpu'
  label 'minMem'
  label 'tabix'
  tag "${meta.id}"

  input:
  tuple val(meta), path(vcf), path(tbi)

  output:
  tuple val(meta), path("*_pass.vcf.gz"), path("*_pass.vcf.gz.tbi"), emit: vcf

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  awk '\$0~"^#" || \$7 == "PASS"{print}' <(bgzip -dc ${vcf}) | bgzip > ${prefix}_pass.vcf.gz
  tabix ${prefix}_pass.vcf.gz
  """
}