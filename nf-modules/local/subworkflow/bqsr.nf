/*
 * BQSR step
 */

include { baseRecalibrator } from '../../common/process/gatk/baseRecalibrator'
include { gatherBQSRReports } from '../../common/process/gatk/gatherBQSRReports'
include { applyBQSR } from '../../common/process/gatk/applyBqsr'
include { samtoolsIndex } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsMerge } from '../../common/process/samtools/samtoolsMerge'

workflow bqsrFlow {

  take:
  bam
  intervals
  dbsnp
  dbsnpIndex
  fasta
  fastaFai
  knownIndels
  knownIndelsIndex
  dict

  main:
  chVersions = Channel.empty()

  // Add intervals size for groupKey function
  chIntervalsNum = intervals.collect().map{it -> [it, it.size() ]}.transpose()

  if (params.noIntervals){
    chBamIntervals = bam.map{ meta, bam, bai -> 
      def newMeta = meta.clone()
      newMeta.numIntervals = 0
      [newMeta,bam,bai,[]]
    }
  }else{
    chBamIntervals = bam.combine(chIntervalsNum)
      .map { meta, bam, bai, intervals, num ->
        def newMeta = meta.clone()
        newMeta.numIntervals = num
        [newMeta, bam, bai, intervals]
      }
  }

  baseRecalibrator(
    chBamIntervals,
    dbsnp,
    dbsnpIndex,
    fasta,
    fastaFai,
    knownIndels,
    knownIndelsIndex,
    dict
  )
  chVersions = chVersions.mix(baseRecalibrator.out.versions)

  // groupKey to speed-up parallel processing
  chRecalTable = baseRecalibrator.out.table
    .map{meta, table ->
      def newMeta = meta.clone()
      [ groupKey(newMeta, meta.numIntervals), table ] 
    }.groupTuple()
    .branch {
      single: it[0].numIntervals <= 1
      multiple: it[0].numIntervals > 1
    }

  gatherBQSRReports(
    chRecalTable.multiple
  )

  chRecalAllTable = gatherBQSRReports.out.table.mix(chRecalTable.single)
  chBaseRecalTable = chBamIntervals
    .combine(chRecalAllTable, by: 0)
    .map{ meta, bam, bai, intervals, table ->
      [meta, bam, bai, table, intervals] 
    }

  applyBQSR(
    chBaseRecalTable,
    fasta,
    fastaFai,
    dict
  )
  chVersions = chVersions.mix(applyBQSR.out.versions)

  // groupKey to speed up parallel processing
  chBQSRbam = applyBQSR.out.bam.map{ meta, bam ->
    [ groupKey(meta, meta.numIntervals), bam ]
    }.groupTuple()

  samtoolsMerge(
    chBQSRbam,
    Channel.value([])
  )
  chVersions = chVersions.mix(samtoolsMerge.out.versions)
  chBamBQSR = params.noIntervals || params.targetBed ? applyBQSR.out.bam : samtoolsMerge.out.bam

  samtoolsIndex(
    chBamBQSR
  )
  chVersions = chVersions.mix(samtoolsIndex.out.versions)

  // Remove no longer necessary meta data
  chBamBai = chBamBQSR
    .join(samtoolsIndex.out.bai)
    .map{meta, bam, bai ->
      def newMeta = [ id: meta.id, name: meta.name, singleEnd: meta.singleEnd ]
      [ newMeta, bam, bai ]
    }

  emit:
  table = chRecalAllTable
  bam = chBamBai
  versions = chVersions
}
