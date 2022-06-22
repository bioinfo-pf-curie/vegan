process sambambaMarkdup {
  tag "${meta.id}"
  label 'sambamba'
  label 'highCpu'
  label 'highMem'

  input:
  tuple val(meta), path(bam), path(bai)

  output:
  tuple val(meta), path("*.md.bam"),path("*.md.bam.bai"), emit: bam
  path ("*.md.flagstats"), emit: metrics
  path('versions.txt'), emit: versions

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  sambamba --version &> versions.txt 2>&1 || true
  sambamba markdup --nthreads ${task.cpus} --tmpdir . ${bam} ${prefix}.md.bam
  sambamba flagstat --nthreads ${task.cpus} ${prefix}.md.bam > ${prefix}.md.flagstats
  """
}
