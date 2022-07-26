/*
 * VCF annotation with snpEff
 */

process snpEff {
  label 'snpeff'
  label 'lowMem'
  label 'lowCpu'

  tag "${fileID}"

  input:
  tuple val(meta), path(vcf)
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
  fileID = "${meta.status}" == "pair" ? "${meta.tumor_id}_vs_${meta.normal_id}" : "${meta.status}" == "tumor" ? "${meta.tumor_id}" : "${meta.normal_id}"
  """
  snpEff -Xmx${task.memory.toGiga()}g \\
    ${db} \\
    ${args} \\
    -csvStats ${fileID}.csv \\
    $cacheCmd \\
    ${vcf[0]} \\
    > ${fileID}.ann.vcf

  bgzip < ${fileID}.ann.vcf > ${fileID}.ann.vcf.gz
  tabix ${fileID}.ann.vcf.gz

  echo \$(snpEff -version | cut -d" " -f1,2) > versions.txt
  echo "tabix "\$(tabix 2>&1 | awk '\$1~"Version"{print \$2}') >> versions.txt
  """
}
