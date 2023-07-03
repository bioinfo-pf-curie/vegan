/* 
 * Load BAM files
 */

include { samtoolsSort } from '../../common/process/samtools/samtoolsSort'
include { samtoolsIndex } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsMerge } from '../../common/process/samtools/samtoolsMerge'
include { samtoolsFlagstat } from '../../common/process/samtools/samtoolsFlagstat'
include { samtoolsStats } from '../../common/process/samtools/samtoolsStats'

// Set the meta.chunk value in case of technical replicates
def setMetaChunk(row){
  def map = []
  row[1].eachWithIndex() { file, i ->
    meta = row[0].clone()
    meta.chunk = i+1
    meta.part = row[1].size()
    map += [meta, file]
  }
  return map
}

workflow loadBamFlow {

  take:
  bams

  main:
  chVersions = Channel.empty()

  chBams = bams.groupTuple()
    .flatMap { it -> setMetaChunk(it) }
    .collate(2)

  // Merge BAM file with the same prefix and use a groupKey to speed up the process
  chBamMapped = chBams
    .map{meta, bam ->
      def newMeta = [ id: meta.id, name: meta.name, singleEnd:meta.singleEnd, part:meta.part ] 
      [ groupKey(newMeta, meta.part), bam ] 
    }.groupTuple()
    .branch {
      single: it[0].part <= 1
      multiple: it[0].part > 1
    }

  samtoolsMerge(
    chBamMapped.multiple
  )
  chVersions = chVersions.mix(samtoolsMerge.out.versions)

  samtoolsSort(
    samtoolsMerge.out.bam
  )
  chVersions = chVersions.mix(samtoolsSort.out.versions)

  samtoolsIndex(
    samtoolsSort.out.bam.mix(chBamMapped.single)
  )
  chVersions = chVersions.mix(samtoolsIndex.out.versions)

  samtoolsFlagstat(
    samtoolsSort.out.bam.mix(chBamMapped.single)
  )
  chVersions = chVersions.mix(samtoolsFlagstat.out.versions)

  samtoolsStats(
    samtoolsSort.out.bam.mix(chBamMapped.single)
  )
  chVersions = chVersions.mix(samtoolsStats.out.versions)

  // Remove groupKey object
  chBamBai = samtoolsSort.out.bam
    .mix(chBamMapped.single)
    .join(samtoolsIndex.out.bai)
    .map{meta, bam, bai -> 
      def newMeta = [ id: meta.id, name: meta.name, singleEnd: meta.singleEnd ]
      [ newMeta, bam, bai ] 
    }

  emit:
  bam = chBamBai
  flagstat = samtoolsFlagstat.out.stats.map{it-> it[1]}
  stats = samtoolsStats.out.stats.map{it-> it[1]}
  versions = chVersions
}
