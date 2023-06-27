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
  tuple val(meta), path(bam), path(bai)
  path(exonBed)
  path(fasta)

  output:
  path("*.mqc"), emit: geneCovMqc
  path("*.pdf"), emit: geneCovOutput

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  def fastaOpts = fasta ? "--fasta ${fasta}" : ''
  """
  mosdepth \
    -n 
    -t ${task.cpus} \
    --by ${exonBed} \
    ${fastaOpts}
    ${meta.id}.genecov \
    ${bam}
  geneCov.r --cov ${meta.id}.genecov.regions.bed.gz --oprefix ${meta.id}_covdensity
  """
}
