/* 
 * Mapping Worflow with BWA
 */

include { bwaMem } from '../../common/process/bwa/bwaMem'
include { bwaMem2 } from '../../common/process/bwamem2/bwaMem2'
include { dragmap } from '../../common/process/dragmap/dragmap'
include { samtoolsSort } from '../../common/process/samtools/samtoolsSort'
include { samtoolsIndex } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsMerge } from '../../common/process/samtools/samtoolsMerge'
include { samtoolsFlagstat } from '../../common/process/samtools/samtoolsFlagstat'

workflow mappingFlow {

  take:
  reads
  index

  main:
  chVersions = Channel.empty()

  if (params.aligner == 'bwa-mem'){
    bwaMem(
      reads,
      index.collect()
    )
    chVersions = chVersions.mix(bwaMem.out.versions)
    chBams = bwaMem.out.bam
    chMappingLogs = bwaMem.out.logs
  }else if (params.aligner == 'bwa-mem2'){
    bwaMem2(
      reads,
      index.collect()
    )
    chVersions = chVersions.mix(bwaMem2.out.versions)
    chBams = bwaMem2.out.bam
    chMappingLogs = bwaMem2.out.logs
  }else if (params.aligner == 'dragmap'){
    dragmap(
      reads,
      index
    )
    chVersions = chVersions.mix(dragmap.out.versions)
    chBams = dragmap.out.bam
    chMappingLogs = dragmap.out.logs
  }

  // Merge BAM file with the same prefix
  chBams.groupTuple(by:[0])
    .branch {
      singleCh: it[1].size() == 1
      multipleCh: it[1].size() > 1
  }.set{bamMapped}

  samtoolsMerge(
    bamMapped.multipleCh
  )
  chVersions = chVersions.mix(samtoolsMerge.out.versions)

  samtoolsSort(
    samtoolsMerge.out.bam.mix(bamMapped.singleCh)
  )
  chVersions = chVersions.mix(samtoolsSort.out.versions)

  samtoolsIndex(
    samtoolsSort.out.bam
  )
  chVersions = chVersions.mix(samtoolsIndex.out.versions)

  samtoolsFlagstat(
    samtoolsSort.out.bam
  )

  emit:
  bam = samtoolsSort.out.bam.join(samtoolsIndex.out.bai)
  logs = chMappingLogs
  stats = samtoolsFlagstat.out.stats.map{it-> it[1]}
  versions = chVersions
}
