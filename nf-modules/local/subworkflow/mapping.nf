/* 
 * Mapping Worflow with BWA
 */

include { bwamem } from '../../common/process/bwa/bwamem'
include { bwamem2 } from '../../common/process/bwamem2/bwamem2'
include { dragmap } from '../../common/process/dragmap/dragmap'
include { samtoolsSort } from '../../common/process/samtools/samtoolsSort'
include { samtoolsIndex } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsMerge } from '../../common/process/samtools/samtoolsMerge'
include { samtoolsFlagstat } from '../../common/process/samtools/samtoolsFlagstat'
include { samtoolsStats } from '../../common/process/samtools/samtoolsStats'

// Set the meta.chunk value in case of technical replicates
def setMetaChunk(row){
  def map = []
  row[1].eachWithIndex() { file,i ->
    meta = row[0].clone()
    meta.chunk = i
    map += [meta, file]
  }
  return map
}

// Remove meta.chunks
def removeChunks(row){
  meta = row[0].clone()
  meta.remove('chunk')
  return [meta, row[1]]
}

workflow mappingFlow {

  take:
  reads
  index

  main:
  chVersions = Channel.empty()

  if (params.splitFastq){
    if (params.singleEnd){
      chReads = reads.map{ it -> [it[0], it[1][0]]}
                     .splitFastq(by: params.fastqChunksSize, file:true, compress:true)
                     .map { it -> [it[0], [it[1]]]}
                     .groupTuple()
                     .flatMap { it -> setMetaChunk(it) }
                     .collate(2)
    }else{
      chReads = reads.map{ it -> [it[0], it[1][0], it[1][1]]}
                     .splitFastq(by: params.fastqChunksSize, pe:true, file:true, compress:true)
                     .map { it -> [it[0], [it[1], it[2]]]}
                     .groupTuple()
                     .flatMap { it -> setMetaChunk(it) }
                     .collate(2)
    }
  }else{
    chReads = reads
  }

  if (params.aligner == 'bwa-mem'){
    bwamem(
      chReads,
      index.collect()
    )
    chVersions = chVersions.mix(bwamem.out.versions)
    chBams = bwamem.out.bam
    chMappingLogs = bwamem.out.logs
  }else if (params.aligner == 'bwa-mem2'){
    bwamem2(
      chReads,
      index.collect()
    )
    chVersions = chVersions.mix(bwamem2.out.versions)
    chBams = bwamem2.out.bam
    chMappingLogs = bwamem2.out.logs
  }else if (params.aligner == 'dragmap'){
    dragmap(
      chReads,
      index.collect()
    )
    chVersions = chVersions.mix(dragmap.out.versions)
    chBams = dragmap.out.bam
    chMappingLogs = dragmap.out.logs
  }

  // Merge BAM file with the same prefix
  chBams
    .map{ it -> removeChunks(it)}
    .groupTuple()
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

  samtoolsStats(
    samtoolsSort.out.bam
  )

  emit:
  bam = samtoolsSort.out.bam.join(samtoolsIndex.out.bai)
  logs = chMappingLogs
  flagstat = samtoolsFlagstat.out.stats.map{it-> it[1]}
  stats = samtoolsStats.out.stats.map{it-> it[1]}
  versions = chVersions
}
