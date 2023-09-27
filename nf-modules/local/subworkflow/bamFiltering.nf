/*
 * Filter BAMs/CRAMs files
 */

include { markDuplicates } from '../../common/process/gatk/markDuplicates'
include { samtoolsView as samtoolsFilter } from '../../common/process/samtools/samtoolsView'
include { samtoolsIndex } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsFlagstat as samtoolsStatsMarkdup } from '../../common/process/samtools/samtoolsFlagstat'
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

    // Mark duplicates
    markDuplicates(
      bam,
      fasta,
      fai
    )
    chVersions = chVersions.mix(markDuplicates.out.versions)

    // Stats on duplicates
    samtoolsStatsMarkdup(
      markDuplicates.out.bam.mix(markDuplicates.out.cram)
    )
    chVersions = chVersions.mix(samtoolsStatsMarkdup.out.versions)

    // Stats on target if any
    samtoolsStatsOnTarget(
      markDuplicates.out.bam.mix(markDuplicates.out.cram),
      bed
    )
    chVersions = chVersions.mix(samtoolsStatsOnTarget.out.versions)

    // Filter reads with samtools
    samtoolsFilter(
      markDuplicates.out.bam.mix(markDuplicates.out.cram)
    )
    chVersions = chVersions.mix(samtoolsFilter.out.versions)

    // Index
    samtoolsIndex(
      samtoolsFilter.out.bam.mix(samtoolsFilter.out.cram)
    )
    chVersions = chVersions.mix(samtoolsIndex.out.versions)

    // Flagstat
    samtoolsStatsFilter(
      samtoolsFilter.out.bam.mix(samtoolsFilter.out.cram)
    )
    chVersions = chVersions.mix(samtoolsStatsFilter.out.versions)

    emit:
    bam = samtoolsFilter.out.bam.join(samtoolsIndex.out.bai)
    cram = samtoolsFilter.out.cram.join(samtoolsIndex.out.crai)
    markdupStats = samtoolsStatsMarkdup.out.stats
    onTargetStats = bed ? samtoolsStatsOnTarget.out.stats : Channel.empty()
    filteringStats  = samtoolsStatsFilter.out.stats
    versions = chVersions
}
