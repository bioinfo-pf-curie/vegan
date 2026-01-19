/*
 * Facets snp-pileup -  Compute SNP pileup at reference positions in one or more input bam files. Input for facets
 */

process facetsPileup {
  label 'facets'
  label 'minCpu'
  label 'medMem'

  tag "${meta.tumor_id}_vs_${meta.normal_id}"

  input:
  tuple val(meta), path(bamTumor), path(baiTumor), path(bamNormal), path(baiNormal)
  path(polym)

  output:
  tuple val(meta), path("${meta.tumor_id}_vs_${meta.normal_id}.csv.gz") ,emit: snppileup

  when:
  task.ext.when == null || task.ext.when

  script:
  """
  snp-pileup \\
    -A -d 1000000 \\
    --gzip \\
    --min-map-quality ${params.mapQual} \\
    --min-base-quality ${params.baseQual} \\
     ${polym} ${meta.tumor_id}_vs_${meta.normal_id}.csv.gz ${bamNormal} ${bamTumor}
  """
}
