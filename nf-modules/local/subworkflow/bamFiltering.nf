/*
 * Filter BAMs file
 */

include { sambambaMarkdup } from '../../common/process/sambamba/sambambaMarkdup'
include { samtoolsFilter } from '../../common/process/samtools/samtoolsFilter'
include { samtoolsIndex as samtoolsIndexTarget } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsIndex as samtoolsIndexFilter } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsFlagstat as samtoolsFlagstatMarkdup } from '../../common/process/samtools/samtoolsFlagstat'
include { samtoolsStatsOnTarget } from '../../common/process/samtools/samtoolsStatsOnTarget'
include { samtoolsFlagstat as samtoolsFlagstatFilter  } from '../../common/process/samtools/samtoolsFlagstat'
include { samtoolsIdxstats } from '../../common/process/samtools/samtoolsIdxstats'
include { samtoolsStats } from '../../common/process/samtools/samtoolsStats'

workflow bamFiltersFlow {

    take:
    bam // [prefix, bam, bai]
    targetBed

    main:
    chVersions = Channel.empty()

    // Remove duplicates
    sambambaMarkdup(
      bam
    )
    chVersions = chVersions.mix(sambambaMarkdup.out.versions)

    samtoolsFlagstatMarkdup(
      sambambaMarkdup.out.bam.map{it->[it[0], it[1]]}
    )

    samtoolsStatsOnTarget(
      sambambaMarkdup.out.bam.map{it->[it[0], it[1]]},
      targetBed
    )

    chVersions = chVersions.mix(samtoolsStatsOnTarget.out.versions)

    // Filter with samtools
    samtoolsFilter(
      sambambaMarkdup.out.bam.map{it->[it[0], it[1]]},
      targetBed
    )

    chVersions = chVersions.mix(samtoolsFilter.out.versions)

    // index
    samtoolsIndexFilter(
      samtoolsFilter.out.bam
    )
    chVersions = chVersions.mix(samtoolsIndexFilter.out.versions)

    // flagstat
    samtoolsFlagstatFilter(
      samtoolsFilter.out.bam
    )
    chVersions = chVersions.mix(samtoolsFlagstatFilter.out.versions)

    // IdxStats
    samtoolsIdxstats(
      samtoolsFilter.out.bam
    )
    chVersions = chVersions.mix(samtoolsIdxstats.out.versions)

    // Stats
    samtoolsStats(
      samtoolsFilter.out.bam
    )
    chVersions = chVersions.mix(samtoolsStats.out.versions)

    emit:
    bam = samtoolsFilter.out.bam.join(samtoolsIndexFilter.out.bai)
    markdupFlagstats = samtoolsFlagstatMarkdup.out.stats.map{it -> it[1]}
    onTargetStats = params.targetBed ? samtoolsStatsOnTarget.out.stats.map{it -> it[1]} : Channel.empty()
    filteringFlagstats  = samtoolsFlagstatFilter.out.stats.map{it -> it[1]}
    idxstats  = samtoolsIdxstats.out.stats
    stats  = samtoolsStats.out.stats
    versions = chVersions
}
