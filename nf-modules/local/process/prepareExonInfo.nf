/*
 * prepareExonInfo:
 * External parameters :
 */

process prepareExonInfo {
  // tag "${meta.id}"
  label 'bedtools'
  label 'minCpu'
  label 'minMem'

  input:
  path(bed)
  path(gtf)

  output:
  path("*exon.bed"), emit: exonBed

  script:
  //def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  awk -F"\t" -v type='gene_id' 'BEGIN{OFS="\t"} \$3=="exon" {split(\$9,annot,";");for(i=1;i<=length(annot);i++){if (annot[i]~type){anntype=annot[i]}} print \$1,\$4-1,\$5,anntype}' ${gtf} | sed -e 's/gene_id//' -e 's/"//g' | sort -T './' -u -k1,1V -k2,2n ${args} > ${gtf.baseName}_exon.bed
  """
}
