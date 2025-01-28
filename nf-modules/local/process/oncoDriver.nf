process oncoDriver {
  label 'minCpu'
  label 'minMem'
  label 'oncodriver'
  tag "${meta.id}"

  input:
  tuple val(meta), path(vcf), path(tbi), path(cnv)
  path(genes)
  path(config)
  path(gtf)

  output:
  tuple val(meta), path("*_oncoDriver.vcf.gz"), emit: vcf
  tuple val(meta), path("*.log"), emit: log
  path ("versions.txt"), emit: versions
  tuple val(meta), path("*csv"), emit: cnv

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  oncoDriver.py \
    --vcf ${vcf} \
    --config ${config} \
    --outputcnv ${prefix}_cnv_oncoDriver.csv \
    --cnv ${cnv} \
    --driver_genes ${genes} \
    --canonical_ids ${gtf} \
    --output ${prefix}_oncoDriver.vcf.gz \
    ${args} 2> ${prefix}_oncoDriver.log

  oncoDriver.py --version > versions.txt
  """
}
