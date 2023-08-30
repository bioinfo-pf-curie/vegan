/*
 * VCF annotation with snpSnift
 */

process snpSiftExtractFieldsPatch {
  label 'snpsift'
  label 'lowMem'
  label 'lowCpu'

  tag "${meta.id}"

  input:
  tuple val(meta), path(vcf)

  output:
  tuple val(meta), path("*.tsv"), emit: tsv
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def args2 = task.ext.args2 ?: ''
  def args3 = task.ext.args3 ?: ''
  def gnomadHeader = task.ext.gnomadHeader ?: ''
  def gnomadExtract = task.ext.gnomadExtract ?: ''
  def db = task.ext.db ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def var = ''
  def j = ''
  """
  SAMPLES=\$(zgrep CHROM ${vcf[0]}  | sed -n 's/^.*FORMAT//p')
  echo \$SAMPLES
  var="" ; for j in \$(echo AD DP AF); do var=\${var}" "\$(echo \${SAMPLES} | sed "s/ /.\${j}\\t/g" | sed "s/\$/.\${j}/g"); done
  echo "Chrom Pos Ref Alt Gene Effect Transcript Exon HGVS_C HGVS_P ${gnomadHeader} ${db} \${var}" | sed "s/ /\\t/g"  > ${prefix}.tsv
  echo \$SAMPLES
  echo \$var
  SnpSift -Xmx${task.memory.toGiga()}g \
    extractFields \
    ${args} \
    ${vcf[0]} \
    ${args2} ${gnomadExtract} ${db} ${args3} | sed 's/;/\\t/g'  | grep -v CHROM | head >> ${prefix}.tsv

  echo "snpSift "\$(SnpSift 2>&1 | awk '\$0~"SnpSift version"{print \$3}') > versions.txt
  echo "tot"
  """
}