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
  tuple val(meta), path("*oncodriver.vcf.gz"), emit: vcf, optional: true
  tuple val(meta), path("*oncodriver.tsv"), emit: cnv, optional: true
  tuple val(meta), path("*.log"), emit: log
  path ("versions.txt"), emit: versions

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  inputCnv = cnv ? "--cnv ${cnv}" : ""
  inputVcf = vcf ? "--vcf ${vcf}" : ""
  inputGtf = gtf ?  "--canonical_ids ${gtf}" : ""
  """
  oncodriver ${inputVcf} \
    ${inputCnv} \
    ${inputGtf} \
    --config ${config} \
    --driver_genes ${genes} \
    --log ${prefix}_oncodriver.log \
    ${args}

  oncodriver --version > versions.txt
  """
}
