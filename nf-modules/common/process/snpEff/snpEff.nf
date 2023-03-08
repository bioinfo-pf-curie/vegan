/*
 * VCF annotation with snpEff
 */

process snpEff {
  label 'snpeff'
  label 'lowMem'
  label 'lowCpu'
  tag "${meta.id}"

  input:
  tuple val(meta), path(vcf), path(index)
  val(db)
  path(cache)

  output:
  tuple val(meta), path("*.ann.vcf.gz*"), emit: vcf
  path("*.csv")                         , emit: report
  path("*.html")                        , emit: summary_html
  path("*.genes.txt")                   , emit: genes_txt
  path("versions.txt")                  , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def cacheCmd = cache ? "-dataDir \${PWD}/${cache}" : ""
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  snpEff -Xmx${task.memory.toGiga()}g \\
    ${db} \\
    ${args} \\
    -csvStats ${prefix}.csv \\
    $cacheCmd \\
    ${vcf[0]} \\
    > ${prefix}.ann.vcf

  bgzip < ${prefix}.ann.vcf > ${prefix}.ann.vcf.gz
  tabix ${prefix}.ann.vcf.gz

  echo \$(snpEff -version | cut -d" " -f1,2) > versions.txt
  echo "tabix "\$(tabix 2>&1 | awk '\$1~"Version"{print \$2}') >> versions.txt
  """
}
