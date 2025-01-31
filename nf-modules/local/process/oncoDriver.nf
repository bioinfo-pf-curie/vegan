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
  tuple val(meta), path("*_oncodriver.vcf.gz"), emit: vcf
  tuple val(meta), path("*_oncodriver.csv"), emit: cnv
  tuple val(meta), path("*.log"), emit: log
  path ("versions.txt"), emit: versions

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  inputCnv = cnv ? "--cnv ${cnv} --outputcnv ${prefix}_cnv_oncodriver.csv" : ""
  inputVcf = vcf ? "--vcf ${vcf} --output ${prefix}_oncodriver.vcf.gz" : ""
  """
  oncoDriver.py \
    ${inputVcf} \
    ${inputCnv} \
    --config ${config} \
    --driver_genes ${genes} \
    --canonical_ids ${gtf} \
    ${args} 2> ${prefix}_oncoDriver.log

  oncoDriver.py --version > versions.txt
  """
}
