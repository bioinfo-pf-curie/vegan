process bamFiltering {
  tag "${meta.id}"
  label 'samtools'
  label 'medCpu'
  label 'minMem'

  input:
  tuple val(meta), path(bam), path(bai)

  output:
  tuple val(meta), path("*.filtered.bam"), path("*.filtered.bam.bai"), emit: chBam
  path("*.filtered.idxstats"), emit: metrics1
  path("*.flagstats"), emit:metrics2
  path("versions.txt"), emit: versions

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  dupParams = ('duplicates' in args) ? "-F 0x0400" : ""
  mapqParams = ('mapq' in args) && (params.mapQual > 0) ? "-q ${params.mapQual}" : ""
  singleParams = ('singleton' in args) ? "-F 0x004 -F 0x008 -f 0x001": "-F 0x004"
  uniqParams =  ('multihits' in args)  ? "-F 0x100 -F 0x800" :  ""
  uniqFilter = ('multihits' in args) ? "| grep -v -e \\\"XA:Z:\\\" -e \\\"SA:Z:\\\" | samtools view -b -" : "| samtools view -b -"
  """
  samtools view -h -@ ${task.cpus} ${uniqParams} ${singleParams} ${dupParams} ${mapqParams} ${bam} ${uniqFilter} > ${prefix}.filtered.bam
  samtools index ${prefix}.filtered.bam
  samtools flagstat ${prefix}.filtered.bam > ${prefix}.filtered.flagstats
  samtools idxstats ${prefix}.filtered.bam > ${prefix}.filtered.idxstats
  samtools stats ${prefix}.filtered.bam > ${prefix}.filtered.stats
  samtools --version &> versions.txt 2>&1 || true
  """
}
