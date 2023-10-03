/*
 * Filter BAMs file
 */

include { sambambaMarkdup } from '../../common/process/sambamba/sambambaMarkdup'
include { samtoolsView as samtoolsFilter } from '../../common/process/samtools/samtoolsView'
include { samtoolsIndex } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsFlagstat as samtoolsStatsMarkdup } from '../../common/process/samtools/samtoolsFlagstat'
include { samtoolsFlagstat as samtoolsStatsFilter  } from '../../common/process/samtools/samtoolsFlagstat'
include { samtoolsStats as samtoolsStatsOnTarget } from '../../common/process/samtools/samtoolsStats'

workflow bamFiltersFlow {

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

    // Stats on duplicates
    samtoolsStatsMarkdup(
      sambambaMarkdup.out.bam.map(it->[it[0], it[1]])
    )
    chVersions = chVersions.mix(samtoolsStatsMarkdup.out.versions)

    // Filter reads with samtools
    samtoolsFilter(
      sambambaMarkdup.out.bam.map(it->[it[0], it[1]])
    )
    chVersions = chVersions.mix(samtoolsFilter.out.versions)

    // Index
    samtoolsIndex(
      samtoolsFilter.out.bam
    )
    chVersions = chVersions.mix(samtoolsIndex.out.versions)

    // Stats after filtering
    samtoolsStatsFilter(
      samtoolsFilter.out.bam
    )
    chVersions = chVersions.mix(samtoolsStatsFilter.out.versions)

    // Stats on target if any
    samtoolsStatsOnTarget(
      samtoolsFilter.out.bam,
      bed
    )
    chVersions = chVersions.mix(samtoolsStatsOnTarget.out.versions)

    emit:
    bam = samtoolsFilter.out.bam.join(samtoolsIndex.out.bai)
    markdupStats = samtoolsStatsMarkdup.out.stats
    onTargetStats = samtoolsStatsOnTarget.out.stats
    filteringStats  = samtoolsStatsFilter.out.stats
    versions = chVersions
}
