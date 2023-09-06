/*
 * Filter on DP, AF and MAF
 */

process filterDpFreqMaf {
  label 'bcftools'
  label 'medCpu'
  label 'medMem'
  tag "${meta.id}"

  input:
  tuple val(meta), path(vcf)

  output:
  tuple val(meta), path("*_filtered.vcf.gz"), path("*_filtered.vcf.gz.tbi"), emit: vcf
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  def args2 = task.ext.args2 ?: ''
  def id = 0
  
  """
  if [ \$(zgrep "tumor_sample" ${vcf[0]} | cut -d "=" -f2) = \$(zgrep "#CHROM" ${vcf[0]} | cut -f10) ] ; then id=0; else id=1; fi
  echo \${id}
  echo \$(zgrep "##tumor_sample" ${vcf[0]} | cut -d '=' -f2)
  echo \$(zgrep "#CHROM" D262E02_T_Mutect2_tagged_norm_pass.vcf.gz | cut -f10)
  bcftools filter \
    -Oz -i '${args} ${args2}' \
    ${vcf[0]} \
    --threads ${task.cpus} \
    -o ${prefix}_filtered.vcf.gz
  tabix ${prefix}_filtered.vcf.gz

  echo \$(bcftools --version | head -1) > versions.txt
  echo "tabix "\$(tabix 2>&1 | awk '\$1~"Version"{print \$2}') >> versions.txt
  """
}