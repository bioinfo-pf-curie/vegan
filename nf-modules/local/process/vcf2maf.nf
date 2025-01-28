/*
 * Export VCF to MAF
 */

process vcf2maf {
  tag "${meta.id}"
  label 'lowCpu'
  label 'medMem'
  label 'vcf2maf'

  input:
  tuple val(meta), path(vcf), path(tbi)
  path(fasta)

  output:
  tuple val(meta), path("*.maf"), emit: maf
//   path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  gunzip -c ${vcf} > ${meta.id}.vcf
  vcf2maf.pl --input-vcf ${meta.id}.vcf --output-maf ${prefix}.maf --ref-fasta ${fasta} --inhibit-vep --vcf-tumor-id ${meta.tumor_id} --ncbi-build GRCh38
  """
}
