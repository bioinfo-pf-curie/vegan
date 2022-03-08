/* 
 * Mapping Worflow with BWA
 */

include { bwaMem } from '../../common/process/bwaMem'
include { samtoolsSort } from '../../common/process/samtoolsSort'
include { samtoolsIndex } from '../../common/process/samtoolsIndex'
include { samtoolsMerge } from '../../common/process/samtoolsMerge'

workflow bwaMapping {

  take:
  reads
  index

  main:
  chVersions = Channel.empty()

  bwaMem(
    reads,
    index.collect()
  )
  chVersions = chVersions.mix(bwaMem.out.versions)

  // Merge BAM file with the same prefix
  bwaMem.out.bam.groupTuple(by:[0])
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

  emit:
  bam = samtoolsSort.out.bam.join(samtoolsIndex.out.bai)
  logs = bwaMem.out.logs
  versions = chVersions
}
