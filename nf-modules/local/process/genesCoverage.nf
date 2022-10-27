/*
 * genesCoverage:
 * External parameters :
 */

process genesCoverage {
  tag "${meta.id}"
  label 'mosdepth'
  label 'minCpu'
  label 'medMem'

  input:
  tuple val(meta), path(bamFiltered), path(baiFiltered)
  path(exonBed)

  output:
  path("*.mqc"), emit: geneCovMqc
  path("*.pdf"), emit: geneCovOutput

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  mosdepth -n -t ${task.cpus} --by ${exonBed} ${meta.id}.genecov ${bamFiltered}
  geneCov.r --cov ${meta.id}.genecov.regions.bed.gz --oprefix ${meta.id}_covdensity
  """
}
