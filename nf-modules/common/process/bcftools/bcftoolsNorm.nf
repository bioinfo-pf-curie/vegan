/*
 * Normalize VCF with bcftools
 */

process bcftoolsNorm {
  label 'bcftools'
  label 'medCpu'
  label 'medMem'
  tag "${meta.id}"

  input:
  tuple val(meta), path(vcf), path(index)
  path(fasta)

  output:
  tuple val(meta), path("*_norm.vcf.gz"),path("*_norm.vcf.gz.tbi"), emit: vcf
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  bcftools norm -Oz -m -both -f ${fasta} --threads ${task.cpus} ${vcf} -o ${prefix}_norm.vcf.gz
  tabix ${prefix}_norm.vcf.gz

  echo \$(bcftools --version | head -1) > versions.txt
  echo "tabix "\$(tabix 2>&1 | awk '\$1~"Version"{print \$2}') >> versions.txt
  """
}
