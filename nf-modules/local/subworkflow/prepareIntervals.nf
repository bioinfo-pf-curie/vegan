/*
 * Prepare intervals for parallel processing
 */

include { buildIntervals } from '../../local/process/buildIntervals'
include { createIntervalBeds } from '../../local/process/createIntervalBeds'

workflow prepareIntervalsFlow {

  take:
  fai
  intervals

  main:
  chVersions = Channel.empty()

  //if (params.noIntervals) {
  //  file("${params.outDir}/no_intervals.bed").text = "no_intervals\n"

  //  chIntervals = Channel.fromPath(file("${params.outDir}/no_intervals.bed"))
  //                       .map{ it -> [it, 0]}

  //  chIntervalsCombined = Channel.fromPath(file("${params.outDir}/no_intervals.bed"))
                                 .map{ it -> [[id:it.simpleName], it]}
  //}else{

    //If no interval file is provided, then intervals are generated from FASTA file
    if (!params.intervals){

      buildIntervals(
        fai
      )
      chVersions = chVersions.mix(buildIntervals.out.versions)

      createIntervalBeds(
        buildIntervals.out.bed
      )
      chVersions = chVersions.mix(createIntervalBeds.out.versions)
      chIntervals = createIntervalBeds.out.bed
      chIntervalsCombined = buildIntervals.out.bed

    }else{
   
      createIntervalBeds(
        intervals
      )
      chVersions = chVersions.mix(createIntervalBeds.out.versions)
      chIntervals = createIntervalBeds.out.bed
      chIntervalsCombined = Channel.fromPath(file(params.intervals)).map{it -> [[id:it.baseName], it] }
    }
  
    // Interval file is split up into multiple bed files for scatter/gather
    chIntervals = chIntervals.flatten()
      .map{ intervalFile ->
            def duration = 0.0
            for (line in intervalFile.readLines()) {
              final fields = line.split('\t')
              if (fields.size() >= 5) duration += fields[4].toFloat()
              else {
                start = fields[1].toInteger()
                end = fields[2].toInteger()
                duration += (end - start) / params.nucleotidesPerSecond
              }
            }
            [duration, intervalFile]
          }.toSortedList({ a, b -> b[0] <=> a[0] })
            .flatten().collate(2)
            .map{duration, intervalFile -> intervalFile}
            //.collect().map{ it -> [it, it.size() ] } // Adding number of intervals as elements
            .transpose()
  //}

  emit:
  versions = chVersions
  intervals = chIntervals
  intervalsCombined = chIntervalsCombined
}
