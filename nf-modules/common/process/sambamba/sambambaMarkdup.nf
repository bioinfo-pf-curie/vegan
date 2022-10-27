process sambambaMarkdup {
  tag "${meta.id}"
  label 'sambamba'
  label 'highCpu'
  label 'highMem'

  input:
  tuple val(meta), path(bam), path(bai)

  output:
  tuple val(meta), path("*.md.bam"), path("*.md.bam.bai"), emit: bam
  path('versions.txt'), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  sambamba markdup 2>&1 | grep ^sambamba > versions.txt
  sambamba markdup --nthreads ${task.cpus} --tmpdir . ${bam} ${prefix}.md.bam
  """
}
