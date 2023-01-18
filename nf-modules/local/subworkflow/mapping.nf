/* 
 * Mapping Worflow with BWA
 */

include { bwamem } from '../../common/process/bwa/bwamem'
include { bwamem2 } from '../../common/process/bwamem2/bwamem2'
include { dragmap } from '../../common/process/dragmap/dragmap'
include { samtoolsView } from '../../common/process/samtools/samtoolsView'
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
    bwamem(
      reads,
      index.collect()
    )
    chVersions = chVersions.mix(bwamem.out.versions)
    chBams = bwamem.out.bam
    chMappingLogs = bwamem.out.logs
  }else if (params.aligner == 'bwa-mem2'){
    bwamem2(
      reads,
      index.collect()
    )
    chVersions = chVersions.mix(bwamem2.out.versions)
    chBams = bwamem2.out.bam
    chMappingLogs = bwamem2.out.logs
  }else if (params.aligner == 'dragmap'){
    dragmap(
      reads,
      index.collect()
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
