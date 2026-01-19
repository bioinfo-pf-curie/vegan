/*
 * Mostdepth - calculate sequencing depth
 */

process mosdepth {
  tag "${meta.id}"
  label 'mosdepth'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path(bamFiltered), path(baiFiltered)
  path(bed)
  path(fasta)

  output:
  tuple val(meta), path('*.summary.txt')          , emit: summary
  tuple val(meta), path('*.global.dist.txt')      , emit: globalTxt
  tuple val(meta), path('*.region.dist.txt')      , optional:true, emit: regionsTxt
  tuple val(meta), path('*.per-base.bed*')        , optional:true, emit: perBaseBed
  tuple val(meta), path('*.regions.bed*')         , optional:true, emit: regionsBed
  tuple val(meta), path('*.quantized.bed*')       , optional:true, emit: quantizedBed
  tuple val(meta), path('*.thresholds.bed*')      , optional:true, emit: thresholdsBed

  path("*.*.txt"), emit: metrics
  path("*{.bed.gz,.bed.gz.csi}"), emit: bedcov
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  bedOpts = bed ? "--by ${bed}" : ''
  fastaOpts = fasta ? "--fasta ${fasta}" : ''
  """
  mosdepth --version &> versions.txt 2>&1 || true
  mosdepth -t ${task.cpus} \
    ${args} \
    ${bedOpts} \
    ${fastaOpts} \
    ${prefix} \
    ${bamFiltered}
  """
}
