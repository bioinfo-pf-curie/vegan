process oncoDrivers {
  label 'minCpu'
  label 'minMem'
  label 'oncodrivers'
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
    --vcf ${vcf} 
    --outputcnv ${prefix}_cnv.csv \
    --cnv ${cnv} \
    --config ${config} \
    --driver_genes ${genes} \
    --canonical_ids ${gtf} \
    --output ${prefix}_oncoDriver.vcf.gz \
    ${args} 2> ${prefix}_oncoDriver.log

  oncoDriver.py --version > versions.txt
  """
}
