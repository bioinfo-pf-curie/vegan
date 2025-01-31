process annotateCNV {
  label 'bedtools'
  label 'minCpu'
  label 'lowMem'

  tag "${meta.id}"

  input:
  tuple val(meta), path(cnv)
  path(bed)

  output:
  tuple val(meta), path("*.annotated.txt"), emit: cnv
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  grep -v "chrom" ${cnv} | awk -F'\t' '{print "chr"\$2, \$3, \$4, \$1, \$5, \$6, \$7, \$8, \$9, \$10}' OFS='\t' > intermediate.txt
  intersectBed -a intermediate.txt -b ${bed} -wao | cut -f1-10,14 > ${prefix}_cnv.transformed.annotated.txt
  echo \$(bedtools --version 2>&1) &> versions.txt
  """
}