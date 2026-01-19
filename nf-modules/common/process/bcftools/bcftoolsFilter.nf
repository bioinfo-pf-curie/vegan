/*
 * Filter VCF file
 */

process bcftoolsFilter {
  label 'bcftools'
  label 'medCpu'
  label 'medMem'
  tag "${meta.id}"

  input:
  tuple val(meta), path(vcf), path(tbi)

  output:
  tuple val(meta), path("*.vcf.gz"), path("*.tbi"), emit: vcf
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}_filtered"
  def args = task.ext.args ?: ''

  """
  id=0
  if [ "${meta.status}" == 'pair' ]; then
    if [ "${meta.tumor_id}" == \$(zgrep "#CHROM" ${vcf[0]} | cut -f10) ]
    then  
      id=0 
    else 
      id=1 
    fi 
  fi

  bcftools view \
    -Oz -i "${args}" \
    ${vcf[0]} \
    --threads ${task.cpus} \
    -o ${prefix}.vcf.gz
  tabix ${prefix}.vcf.gz

  echo \$(bcftools --version | head -1) > versions.txt
  echo "tabix "\$(tabix 2>&1 | awk '\$1~"Version"{print \$2}') >> versions.txt
  """
}