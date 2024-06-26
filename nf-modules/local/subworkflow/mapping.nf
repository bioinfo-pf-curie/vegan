/* 
 * Mapping Worflow
 */

include { bwaMem } from '../../common/process/bwa/bwaMem'
include { bwaMem2 } from '../../common/process/bwaMem2/bwaMem2'
include { dragmap } from '../../common/process/dragmap/dragmap'
include { samtoolsSort } from '../../common/process/samtools/samtoolsSort'
include { samtoolsIndex } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsMerge } from '../../common/process/samtools/samtoolsMerge'
include { samtoolsFlagstat } from '../../common/process/samtools/samtoolsFlagstat'
include { samtoolsStats } from '../../common/process/samtools/samtoolsStats'
include { samtoolsView as samtoolsConvert } from '../../common/process/samtools/samtoolsView'
include { samtoolsIndex as samtoolsConvertIndex } from '../../common/process/samtools/samtoolsIndex'
  
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

workflow mappingFlow {

  take:
  reads
  index
  fasta
  fai

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

  if (params.aligner == 'bwa-mem'){

    bwaMem(
      chReads,
      index.collect(),
      Channel.of(true).collect()
    )
    chVersions = chVersions.mix(bwaMem.out.versions)
    chBams = bwaMem.out.bam

  }else if (params.aligner == 'bwa-mem2'){

    bwaMem2(
      chReads,
      index.collect(),
      Channel.of(true).collect()
    )
    chVersions = chVersions.mix(bwaMem2.out.versions)
    chBams = bwaMem2.out.bam

  }else if (params.aligner == 'dragmap'){

    dragmap(
      chReads,
      index.collect(),
      Channel.of(true).collect()
    )
    chVersions = chVersions.mix(dragmap.out.versions)
    chBams = dragmap.out.bam

  }


  // Merge BAM file with the same prefix and use a groupKey to speed up the process
  chBamMapped = chBams
    .map{meta, bam ->
      def newMeta = [ id: meta.id, name: meta.name, singleEnd: meta.singleEnd, part:meta.part ] 
      [ groupKey(newMeta, meta.part), bam ] 
    }.groupTuple()
    .branch {
      single: it[0].part <= 1
      multiple: it[0].part > 1
    }

  samtoolsMerge(
    chBamMapped.multiple,
    Channel.value([])
  )
  chVersions = chVersions.mix(samtoolsMerge.out.versions)

  samtoolsIndex(
    samtoolsMerge.out.bam.mix(chBamMapped.single)
  )
  chVersions = chVersions.mix(samtoolsIndex.out.versions)

  samtoolsFlagstat(
    samtoolsMerge.out.bam.mix(chBamMapped.single)
  )
  chVersions = chVersions.mix(samtoolsFlagstat.out.versions)

  samtoolsStats(
    samtoolsMerge.out.bam.mix(chBamMapped.single),
    Channel.value([])
  )
  chVersions = chVersions.mix(samtoolsStats.out.versions)

  // Filter all 'aligned' channels that fail the check
  // And remove groupKey object
  chBamBai = samtoolsFlagstat.out.stats
    .join(samtoolsMerge.out.bam.mix(chBamMapped.single))
    .join(samtoolsIndex.out.bai)
    .filter { meta, logs, bam, bai -> checkAlignmentPercent(meta, logs) }
    .map { meta, logs, bam, bai -> 
      def newMeta = [ id: meta.id, name: meta.name, singleEnd: meta.singleEnd ]
      [ newMeta, bam, bai ]
    }

  // Convert BAM to CRAM
  chCramCrai = Channel.empty()
  if (params.cram){
    samtoolsConvert(
      chBamBai.map{it->[it[0], it[1]]},
      fasta
    )
    samtoolsConvertIndex(
      samtoolsConvert.out.cram
    )
    chCramCrai = samtoolsConvert.out.cram.join(samtoolsConvertIndex.out.crai)
  }
                                                    
  emit:
  bam = chBamBai
  cram = chCramCrai
  flagstat = samtoolsFlagstat.out.stats
  stats = samtoolsStats.out.stats
  versions = chVersions
}
