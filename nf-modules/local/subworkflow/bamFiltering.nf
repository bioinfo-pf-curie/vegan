/*
 * Filter BAMs/CRAMs files
 */

//include { markDuplicates } from '../../common/process/gatk/markDuplicates'
include { samtoolsView as samtoolsFilter } from '../../common/process/samtools/samtoolsView'
include { samtoolsIndex } from '../../common/process/samtools/samtoolsIndex'
//include { samtoolsFlagstat as samtoolsStatsMarkdup } from '../../common/process/samtools/samtoolsFlagstat'
include { samtoolsFlagstat as samtoolsStatsFilter  } from '../../common/process/samtools/samtoolsFlagstat'
include { samtoolsStats as samtoolsStatsOnTarget } from '../../common/process/samtools/samtoolsStats'

workflow bamFiltersFlow {

    take:
    bam // [prefix, bam, bai]
    bed
    fasta
    fai

    main:
    chVersions = Channel.empty()

    // Filter reads with samtools
    samtoolsFilter(
      bam.map{it->[it[0], it[1]]},
      Channel.value([])
    )
    chVersions = chVersions.mix(samtoolsFilter.out.versions)

    // Index
    samtoolsIndex(
      samtoolsFilter.out.bam.mix(samtoolsFilter.out.cram)
    )
    chVersions = chVersions.mix(samtoolsIndex.out.versions)

    // Stats after filtering
    samtoolsStatsFilter(
      samtoolsFilter.out.bam.mix(samtoolsFilter.out.cram)
    )
    chVersions = chVersions.mix(samtoolsStatsFilter.out.versions)

    // Stats on target if any
    samtoolsStatsOnTarget(
      samtoolsFilter.out.bam.mix(samtoolsFilter.out.cram),
      bed
    )
    chVersions = chVersions.mix(samtoolsStatsOnTarget.out.versions)

    emit:
    bam = samtoolsFilter.out.bam.join(samtoolsIndex.out.bai)
    cram = samtoolsFilter.out.cram.join(samtoolsIndex.out.crai)
    onTargetStats = samtoolsStatsOnTarget.out.stats
    filteringStats  = samtoolsStatsFilter.out.stats
    versions = chVersions
}
