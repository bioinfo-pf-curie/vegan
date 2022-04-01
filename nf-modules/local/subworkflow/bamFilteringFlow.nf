/*
 * Filter BAMs file
 */

include { markDuplicates } from '../../local/process/markDuplicates'
include { bamOnTarget } from '../../local/process/bamOnTarget'

workflow bamFiltering {

    take:
    bams // [prefix, bam, bai]
    bed

    main:
    chVersions = Channel.empty()

    // Remove duplicates
    markDuplicates(
      bams
    )
    chVersions = chVersions.mix(markDuplicates.out.versions)

    if(params.targetBed){
      bamOnTarget(
          markDuplicates.out.bam,
          bed

          )
          chVersions = chVersions.mix(bamOnTarget.out.versions)
          chBam = bamOnTarget.out.bam
    }else{
      chBam = markDuplicates.out.bam
    }

    emit:
    bam = markDuplicates.out.bam
    metrics  = markDuplicates.out.metrics
    versions = chVersions
}
