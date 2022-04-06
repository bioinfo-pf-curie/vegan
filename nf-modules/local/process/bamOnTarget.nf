process bamOnTarget {
  tag "${meta.id}"
  label 'bedtools'
  label 'minCpu'
  label 'lowMem'

  input:
  tuple val(meta), path(bam), path(bai)
  path(targetBed)

  output:
  tuple val(meta), path("*.onTarget.bam"), emit: bam
  path("*.onTarget.flagstats"), emit: metrics
  path("versions.txt"), emit: versions

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  echo \$(bedtools --version 2>&1) &> versions.txt
  intersectBed -abam ${bam} -b ${targetBed} > ${prefix}.onTarget.bam
  samtools index ${prefix}.onTarget.bam
  samtools flagstat ${prefix}.onTarget.bam > ${prefix}.onTarget.flagstats
  """
}
