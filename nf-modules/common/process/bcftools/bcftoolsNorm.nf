/*
 * Normalize VCF with bcftools
 */

process bcftoolsNorm {
  label 'bcftools'
  label 'medCpu'
  label 'medMem'
  tag "${fileID}"

  input:
  tuple val(meta), path(vcf)
  path(fasta)

  output:
  tuple val(meta), path("*_norm.vcf.gz"),path("*_norm.vcf.gz.tbi"), emit: vcf
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  //prefix = task.ext.prefix ?: "${meta.id}"

  fileID = "${meta.status}" == "pair" ? "${meta.tumor_id}_vs_${meta.normal_id}" : "${meta.status}" == "tumor" ? "${meta.tumor_id}" : "${meta.normal_id}"

  """
  bcftools norm -Oz -m -both -f ${fasta} --threads ${task.cpus} ${vcf[0]} -o ${fileID}_norm.vcf.gz
  tabix ${fileID}_norm.vcf.gz

  echo \$(bcftools --version | head -1) > versions.txt
  echo "tabix "\$(tabix 2>&1 | awk '\$1~"Version"{print \$2}') >> versions.txt
  """
}
