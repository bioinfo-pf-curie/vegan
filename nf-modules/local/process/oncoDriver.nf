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
  tuple val(meta), path("*_snv_oncodriver.vcf.gz"), emit: vcf, optional: true
  tuple val(meta), path("*_cnv_oncodriver.csv"), emit: cnv, optional: true
  tuple val(meta), path("*.log"), emit: log
  path ("versions.txt"), emit: versions

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  inputCnv = cnv ? "--cnv ${cnv} --outputcnv ${prefix}_cnv_oncodriver.csv" : ""
  inputVcf = vcf ? "--vcf ${vcf} --output ${prefix}_snv_oncodriver.vcf.gz" : ""
  inputGtf = gtf ?  "--canonical_ids ${gtf}" : ""
  """
  oncoDriver.py \
    ${inputVcf} \
    ${inputCnv} \
    ${inputGtf} \
    --config ${config} \
    --driver_genes ${genes} \
    ${args} 2> ${prefix}_oncoDriver.log

  oncoDriver.py --version > versions.txt
  """
}
