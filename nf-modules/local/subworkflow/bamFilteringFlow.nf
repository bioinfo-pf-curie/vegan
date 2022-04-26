/*
 * Filter BAMs file
 */

include { sambambaMarkdup } from '../../local/process/sambambaMarkdup'
include { bamOnTarget } from '../../local/process/bamOnTarget'
include { samtoolsFilter } from '../../common/process/samtoolsFilter'
include { samtoolsIndex } from '../../common/process/samtoolsIndex'
include { samtoolsFlagstat } from '../../common/process/samtoolsFlagstat'
include { samtoolsIdxstats } from '../../common/process/samtoolsIdxstats'
include { samtoolsStats } from '../../common/process/samtoolsStats'
// include { samtoolsFiltering as samtoolsFilterigSV} from ''
// include { samtoolsFiltering as samtoolsFilterigSNV} from ''
// Permet d'avoir des config different dans modulues.config

workflow bamFilters {

    take:
    bam // [prefix, bam, bai]
    bed

    main:
    chVersions = Channel.empty()

    // Remove duplicates
    sambambaMarkdup(
      bam
    )
    chVersions = chVersions.mix(sambambaMarkdup.out.versions)

    chBamTest = sambambaMarkdup.out.bam
    chBamTest.view()

    // Reduce to the Target
    if(params.targetBed){
      bamOnTarget(
          sambambaMarkdup.out.bam,
          bed.collect()
          )
          chVersions = chVersions.mix(bamOnTarget.out.versions)
          chBam = bamOnTarget.out.bam
    }else{
      chBam = sambambaMarkdup.out.bam
    }
    
    // Filter with samtools
    samtoolsFilter(
      chBam
      )
      chBam = samtoolsFilter.out.bam
      chVersions = chVersions.mix(samtoolsFilter.out.versions)

    // index
    samtoolsIndex(
      chBam
    )

    // flagstat
    samtoolsFlagstat(
      chBam
    )

    // IdxStats
    samtoolsIdxstats(
      chBam
    )

    // Stats
    samtoolsStats(
      chBam
    )

    emit:
    bam = samtoolsFilter.out.bam.join(samtoolsIndex.out.bai)
    markdupMetrics = sambambaMarkdup.out.metrics
    bamOnTargetMetrics = bamOnTarget.out.metrics
    flagstat  = samtoolsFlagstat.out.stats
    idxstats  = samtoolsIdxstats.out.stats
    stats  = samtoolsStats.out.stats
    versions = chVersions
}
