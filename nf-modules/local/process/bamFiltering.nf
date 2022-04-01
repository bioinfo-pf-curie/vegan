process bamFiltering {
  tag "${meta.id}"
  label 'samtools'
  label 'medCpu'
  label 'minMem'
  //tag "${sampleId}-${vCType}"

  publishDir "${params.filteredBamDir}", mode: params.publishDirMode

  input:
  tuple val(meta), path(bam), path(bai)
  val(vCType)

  output:
  tuple val(sampleId), val(sampleName), val(vCType), file("${sampleId}.filtered.${vCType}.bam"), file("${sampleId}.filtered.${vCType}.bam.bai") into filteredBamCh
  file("${sampleId}.filtered.${vCType}.idxstats") into bamFilterReportCh
  file("*.flagstats") into filteringReportCh
  file('v_samtools.txt') into samtoolsBamFilterVersionCh

  script:
  dupParams = (vCType == 'SNV' && 'duplicates' in SNVFilters) | (vCType == 'SV' && 'duplicates' in SVFilters) ? "-F 0x0400" : ""
  mapqParams = (vCType == 'SNV' && 'mapq' in SNVFilters) | (vCType == 'SV' && 'mapq' in SVFilters) && (params.mapQual > 0) ? "-q ${params.mapQual}" : ""
  singleParams = (vCType == 'SNV' && 'singleton' in SNVFilters) | (vCType == 'SV' && 'single' in SVFilters) ? "-F 0x004 -F 0x008 -f 0x001": "-F 0x004"
  uniqParams =  (vCType == 'SNV' && 'multihits' in SNVFilters) | (vCType == 'SV' && 'multi' in SVFilters) ? "-F 0x100 -F 0x800" :  ""
  uniqFilter = (vCType == 'SNV' && 'multihits' in SNVFilters) | (vCType == 'SV' && 'multi' in SVFilters) ? "| grep -v -e \\\"XA:Z:\\\" -e \\\"SA:Z:\\\" | samtools view -b -" : "| samtools view -b -"
  """
  samtools view -h -@ ${task.cpus} ${uniqParams} ${singleParams} ${dupParams} ${mapqParams} ${bam} ${uniqFilter} > ${sampleId}.filtered.${vCType}.bam
  samtools index ${sampleId}.filtered.${vCType}.bam
  samtools flagstat ${sampleId}.filtered.${vCType}.bam > ${sampleId}.filtered.${vCType}.flagstats
  samtools idxstats ${sampleId}.filtered.${vCType}.bam > ${sampleId}.filtered.${vCType}.idxstats
  samtools stats ${sampleId}.filtered.${vCType}.bam > ${sampleId}.filtered.${vCType}.stats
  samtools --version &> v_samtools.txt 2>&1 || true
  """
}
