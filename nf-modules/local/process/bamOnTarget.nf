process bamOnTarget {
  tag "${meta.id}"
  label 'bedtools'
  label 'minCpu'
  label 'lowMem'

  input:
  tuple val(meta), path(bam), path(bai)
  path(targetBed)

  output:
  tuple val(meta), path("*_onTarget.bam"), path("*_onTarget.bam.bai"), emit: bam
  path("*_onTarget.flagstats"), emit: metrics
  path("versions.txt"), emit: versions

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  echo \$(bedtools --version 2>&1) &> versions.txt
  intersectBed -abam ${bam} -b ${targetBed} > ${prefix}_onTarget.bam
  samtools index ${prefix}_onTarget.bam
  samtools flagstat ${prefix}_onTarget.bam > ${prefix}_onTarget.flagstats
  """
}
