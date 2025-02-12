/* 
 * STAR Workflow
 */

include { starAlign } from '../../common/process/star/starAlign'
include { samtoolsSort } from '../../common/process/samtools/samtoolsSort'
include { samtoolsIndex } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsFlagstat } from '../../common/process/samtools/samtoolsFlagstat'

include {checkAlignmentPercent} from '../../../lib/functions'

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

workflow mappingStarFlow {

  take:
  reads // Channel [val(meta), [reads]]
  index // Channel path(index)
  gtf // Channel path(gtf)

  main:
  
  chVersions = Channel.empty()

  // reads always contains a chunk information
  if (params.splitFastq){
    if (params.singleEnd){
      chReads = reads.map{ it -> [it[0], it[1][0]]}
        .splitFastq(by: params.fastqChunksSize, file:true, compress:false)
        .map { it -> [it[0], [it[1]]]}
        .groupTuple()
        .flatMap { it -> setMetaChunk(it) }
        .collate(2)
    }else{
      chReads = reads.map{ it -> [it[0], it[1][0], it[1][1]]}
        .splitFastq(by: params.fastqChunksSize, pe:true, file:true, compress:false)
        .map { it -> [it[0], [it[1], it[2]]]}
        .groupTuple()
        .flatMap { it -> setMetaChunk(it) }
        .collate(2)
    }
  }else{
    chReads = reads.groupTuple()
      .flatMap { it -> setMetaChunk(it) }
      .collate(2)
  }

  starAlign(
    chReads,
    index.collect(),
    gtf.collect().ifEmpty([])
  )
  chVersions = chVersions.mix(starAlign.out.versions)

  samtoolsSort(
    starAlign.out.bam
  )
  chVersions = chVersions.mix(samtoolsSort.out.versions)
 
  samtoolsIndex(
    samtoolsSort.out.bam
  )
  chVersions = chVersions.mix(samtoolsIndex.out.versions)

  samtoolsFlagstat(
    samtoolsSort.out.bam
  )
  chVersions = chVersions.mix(samtoolsFlagstat.out.versions)

  // Filter all 'aligned' channels that fail the check
  // And remove groupKey object
  chBamBai = samtoolsFlagstat.out.stats
    .join(samtoolsSort.out.bam)
    .join(samtoolsIndex.out.bai)
    .filter { meta, logs, bam, bai -> checkAlignmentPercent(meta, logs) }
    .map { meta, logs, bam, bai -> 
      def newMeta = [ id: meta.id, name: meta.name]
      [ newMeta, bam, bai ]
    }

  emit:
  bam = chBamBai
  transcriptsBam = starAlign.out.transcriptsBam
  flagstat = samtoolsFlagstat.out.stats
  logs = starAlign.out.logs
  counts = starAlign.out.counts
  countsLogs = starAlign.out.countsLogs
  versions = chVersions
}
