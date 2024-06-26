/*
 * prepareExonInfo
 */

process prepareExonInfo {
  label 'bedtools'
  label 'minCpu'
  label 'minMem'

  input:
  path(gtf)
  path(bed)

  output:
  path("*exon.bed"), emit: exonBed

  script:
  def args = task.ext.args ?: ''
  bedCmd = bed ? " | intersectBed -a stdin -b ${bed} " : ''
  """
  awk -F"\t" -v type='gene_id' 'BEGIN{OFS="\t"} \$3=="exon" {split(\$9,annot,";");for(i=1;i<=length(annot);i++){if (annot[i]~type){anntype=annot[i]}} print \$1,\$4-1,\$5,anntype}' ${gtf} | sed -e 's/gene_id//' -e 's/"//g' | sort -T './' -u -k1,1V -k2,2n ${bedCmd} > ${gtf.baseName}_exon.bed
  """
}
