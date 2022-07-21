/*
 * Filter BAMs file
 */

include { sambambaMarkdup } from '../../common/process/sambamba/sambambaMarkdup'
include { intersectBed } from '../../common/process/bedtools/intersectBed'
include { samtoolsFilter } from '../../common/process/samtools/samtoolsFilter'
include { samtoolsIndex as samtoolsIndexTarget } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsIndex as samtoolsIndexFilter } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsFlagstat as samtoolsFlagstatTarget } from '../../common/process/samtools/samtoolsFlagstat'
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


    // Reduce to the Target
    if(params.targetBed){

      intersectBed(
        sambambaMarkdup.out.bam,
        bed.collect()
      )

      samtoolsIndexTarget(
        intersectBed.out.bam
      )

      samtoolsFlagstatTarget(
        intersectBed.out.bam
      )

      chVersions = chVersions.mix(intersectBed.out.versions)
      chBam = intersectBed.out.bam
    }else{
      chBam = sambambaMarkdup.out.bam
    }

    // Filter with samtools
    samtoolsFilter(
      chBam
    )
    chVersions = chVersions.mix(samtoolsFilter.out.versions)

    // index
    samtoolsIndexFilter(
      samtoolsFilter.out.bam
    )

    // flagstat
    samtoolsFlagstatFilter(
      samtoolsFilter.out.bam
    )

    // IdxStats
    samtoolsIdxstats(
      samtoolsFilter.out.bam
    )

    // Stats
    samtoolsStats(
      samtoolsFilter.out.bam
    )

    emit:
    bam = samtoolsFilter.out.bam.join(samtoolsIndexFilter.out.bai)
    markdupMetrics = sambambaMarkdup.out.metrics
    flagstat  = samtoolsFlagstatFilter.out.stats
    idxstats  = samtoolsIdxstats.out.stats
    stats  = samtoolsStats.out.stats
    versions = chVersions
}
