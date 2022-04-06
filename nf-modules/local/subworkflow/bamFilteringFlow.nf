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
    bams // [prefix, bam, bai]
    bed

    main:
    chVersions = Channel.empty()

    // Remove duplicates
    sambambaMarkdup(
      bams
    )
    chVersions = chVersions.mix(sambambaMarkdup.out.versions)

    // Reduce to the Target
    if(params.targetBed){
      bamOnTarget(
          sambambaMarkdup.out.bam,
          bed
          )
          chVersions = chVersions.mix(bamOnTarget.out.versions)
          chBam = bamOnTarget.out.bam
    }else{
      chBam = sambambaMarkdup.out.bam
    }

chBam.view()

    // Filter with samtools
    samtoolsFilter(
      chBam
      )
      chVersions = chVersions.mix(samtoolsFilter.out.versions)

    // flagstat
    samtoolsFlagstat(
      chBam
    )

    // IdxStats
    samtoolsIdxstats(
      chBam
      //samtoolsIndex.out.bai
    )

    // Stats
    samtoolsStats(
      chBam
    )

    emit:
    bam = sambambaMarkdup.out.bam
    markdupMetrics = sambambaMarkdup.out.metrics
    bamOnTargetMetrics = bamOnTarget.out.metrics
    flagstat  = samtoolsFlagstat.out.stats
    idxstats  = samtoolsIdxstats.out.stats
    stats  = samtoolsStats.out.stats
    versions = chVersions
}
