/*
 * Preseq - Saturation Curves
 * External parameters :
 * @ params.preseqDefect : run preseq in defect mode
 */

process preseq {
  tag "${bam}"
  label 'preseq'
  label 'minCpu'
  label 'highMem'

  input:
  tuple val(meta), path(bam), path(bai)

  output:
  path("*ccurve.txt"), emit: results
  path("versions.txt"), emit: versions

  script:
  defectMode = params.preseqDefect || task.attempt > 1 ? '-D' : ''
  peOpts = meta.singleEnd ? '' : '-pe'
  """
  echo \$(preseq 2>&1 | awk '\$0~"Version"{print "Preseq",\$2}') > versions.txt
  preseq lc_extrap -seed 1 -v -B ${bam} ${peOpts} ${defectMode} -o ${meta.id}_extrap_ccurve.txt -e 200e+06 -seg_len 100000000
  """
}

