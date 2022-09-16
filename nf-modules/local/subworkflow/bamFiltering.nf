/*
 * Filter BAMs file
 */

include { sambambaMarkdup } from '../../common/process/sambamba/sambambaMarkdup'
include { intersectBed } from '../../common/process/bedtools/intersectBed'
include { samtoolsFilter } from '../../common/process/samtools/samtoolsFilter'
include { samtoolsIndex as samtoolsIndexTarget } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsIndex as samtoolsIndexFilter } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsFlagstat as samtoolsFlagstatMarkdup } from '../../common/process/samtools/samtoolsFlagstat'
include { samtoolsFlagstat as samtoolsFlagstatOnTarget } from '../../common/process/samtools/samtoolsFlagstat'
include { samtoolsFlagstat as samtoolsFlagstatFilter  } from '../../common/process/samtools/samtoolsFlagstat'
include { samtoolsIdxstats } from '../../common/process/samtools/samtoolsIdxstats'
include { samtoolsStats } from '../../common/process/samtools/samtoolsStats'

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

    samtoolsFlagstatMarkdup(
      sambambaMarkdup.out.bam.map{it->[it[0], it[1]]}
    )

    // Reduce to the target for WES analysis
    intersectBed(
      sambambaMarkdup.out.bam,
      bed
    )
    chVersions = chVersions.mix(intersectBed.out.versions)

    samtoolsFlagstatOnTarget(
      intersectBed.out.bam
    )
    chVersions = chVersions.mix(samtoolsFlagstatOnTarget.out.versions)
    chBam = params.targetBed ? intersectBed.out.bam : sambambaMarkdup.out.bam

    // Filter with samtools
    samtoolsFilter(
      chBam
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
    onTargetFlagstats = samtoolsFlagstatOnTarget.out.stats.map{it -> it[1]}
    filteringFlagstats  = samtoolsFlagstatFilter.out.stats.map{it -> it[1]}
    idxstats  = samtoolsIdxstats.out.stats
    stats  = samtoolsStats.out.stats
    versions = chVersions
}
