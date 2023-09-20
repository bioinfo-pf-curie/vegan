/*
 * Filter on DP, AF and MAF
 */

process bcftoolsFilter {
  label 'bcftools'
  label 'medCpu'
  label 'medMem'
  tag "${meta.id}"

  input:
  tuple val(meta), path(vcf)

  output:
  tuple val(meta), path("*.vcf{.gz,.gz.tbi}"), emit: vcf
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  def id = 0
  
  """
  if [ ${meta.tumor_id} == \$(zgrep "#CHROM" ${vcf[0]} | cut -f10) ] ; then id=0; else id=1; fi
  bcftools view \
    -Oz -i '${args}' \
    ${vcf[0]} \
    --threads ${task.cpus} \
    -o ${prefix}_filtered.vcf.gz
  tabix ${prefix}_filtered.vcf.gz

  echo \$(bcftools --version | head -1) > versions.txt
  echo "tabix "\$(tabix 2>&1 | awk '\$1~"Version"{print \$2}') >> versions.txt
  """
}