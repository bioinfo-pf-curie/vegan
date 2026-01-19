/*
 * Rescue missed variants and merge with the PASS vcf
 */

process rescueVariants {
  label 'bcftools'
  label 'medCpu'
  label 'medMem'
  tag "${meta.id}"

  input:
  tuple val(meta), path(vcfRaw), path(tbi)
  tuple val(meta), path(vcfFiltered), path(tbi)

  output:
  tuple val(meta), path("*.vcf.gz"), path("*.tbi"), emit: vcf
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''

  """
  gunzip -c ${vcfFiltered} > ${prefix}_merge.vcf
  if [ ${meta.tumor_id} == \$(zgrep "#CHROM" ${vcfFiltered[0]} | cut -f10) ] ; then id=0; else id=1; fi
  bcftools view ${vcfRaw} -i "${args}" -H -Ov >> ${prefix}_merge.vcf
  bcftools sort ${prefix}_merge.vcf -Oz -o ${prefix}_merge.vcf.gz
  tabix ${prefix}_merge.vcf.gz

  echo \$(bcftools --version | head -1) > versions.txt
  echo "tabix "\$(tabix 2>&1 | awk '\$1~"Version"{print \$2}') >> versions.txt
  """
}