/*
 * RNA Flow
 */
include { loadDesign } from '../../../lib/functions'
include {checkTumorOnly} from '../../../lib/functions'
include { mappingStarFlow } from '../../local/subworkflow/mappingStar'
include { preseq } from '../../common/process/preseq/preseq'
include { markDuplicates } from '../../common/process/gatk/markDuplicates'
include { samtoolsIndex as samtoolsIndexMarkdup } from '../../common/process/samtools/samtoolsIndex'
include { samtoolsStats as samtoolsStatsMarkdup } from '../../common/process/samtools/samtoolsStats'
include { splitNCigarReads } from '../../common/process/gatk/splitNCigarReads'
include { identitoFlow } from '../../common/subworkflow/identito'
include { bqsrFlow } from '../../local/subworkflow/bqsr'
include { haplotypeCallerFlow } from '../../local/subworkflow/haplotypeCaller'
include { variantFiltration } from '../../common/process/gatk/variantFiltration'

workflow rnaFlow {

    take:
    reads
    index
    gtf
    fasta
    fai
    dict
    polymBed
    intervals
    dbsnp
    dbsnpIndex
    knownIndels
    knownIndelsIndex

    main:

    tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase()} : []

    mappingStarFlow(
      reads,
      index,
      gtf
    )
    
    chVersions = Channel.empty()
    chAlignedBam = mappingStarFlow.out.bam
    chFilteredBam = Channel.empty()
    chMarkedBam = Channel.empty()
    // chMappingStats = mappingStarFlow.out.stats.map{it->it[1]}
    chMappingFlagstat = mappingStarFlow.out.flagstat.map{it->it[1]}
    chVersions = chVersions.mix(mappingStarFlow.out.versions)

    // preseq(
    //   chAlignedBam.mix(chMarkedBam, chFilteredBam)
    // )

    // chPreseqMqc = preseq.out.results.collect()
    // chVersions = chVersions.mix(preseq.out.versions)
    
    markDuplicates(
      chAlignedBam,
      fasta,
      fai
    )
    chVersions = chVersions.mix(markDuplicates.out.versions)

    // Index
    samtoolsIndexMarkdup(
      markDuplicates.out.bam.mix(markDuplicates.out.cram)
    )
    chVersions = chVersions.mix(samtoolsIndexMarkdup.out.versions)

    chMarkedBam = markDuplicates.out.bam.mix(markDuplicates.out.cram)
        .join(samtoolsIndexMarkdup.out.bai) 

    // Stats on duplicates
    // samtoolsStatsMarkdup(
    //   chMarkedBam.map{it->[it[0],it[1]]}
    // )
    // chMarkdupStatsMqc = samtoolsStatsMarkdup.out.stats.map{it->it[1]}
    // chVersions = chVersions.mix(samtoolsStatsMarkdup.out.versions)

    splitNCigarReads(
        chMarkedBam,
        fasta,
        fai,
        dict
    )
    chVersions = chVersions.mix(splitNCigarReads.out.versions)
    chFilteredBam = splitNCigarReads.out.bam.join(splitNCigarReads.out.bai)

    if (!params.skipIdentito){
      identitoFlow(
        chFilteredBam,
        fasta,
        fai,
        polymBed
      )
      chIdentitoMqc = identitoFlow.out.tsv
      chVersions = chVersions.mix(identitoFlow.out.versions)
    }

    chProcBam = chFilteredBam
      bqsrFlow(
        chFilteredBam,
        intervals,
        dbsnp,
        dbsnpIndex,
        fasta,
        fai,
        knownIndels,
        knownIndelsIndex,
        dict
      )
      chVersions = chVersions.mix(bqsrFlow.out.versions)
      chProcBam = params.skipBQSR ? chFilteredBam : bqsrFlow.out.bam.mix(bqsrFlow.out.cram)

    chSingleBam = Channel.empty()
    chPairBam = Channel.empty()
    chTumorOnlyBam = Channel.empty()
    chGermlineOnlyBam = Channel.empty()

    // Load design file
    if (params.design){
    Channel
        .fromPath(params.design)
        .ifEmpty { exit 1, "Design file not found: ${params.design}" }
        .set { chDesignFile }

    chDesign = loadDesign(params.design)

    //Separate the design in paired / tumor only / germline only
    chDesign = chDesign
        .branch{
        paired: it[0] != '' && it[1] != ''
        tumorOnly: it[0] != '' && it[1] == ''
        germlineOnly: it[0] == '' && it[1] != ''
    }
    // Check if PON is defined for tumor only samples
    chDesign.tumorOnly.map{it->it[0]}.toList().map{ it -> checkTumorOnly(it, params, tools) }
    }else{
    chDesignFile = Channel.empty()
    chDesign = Channel.empty()
    }

    if (params.design && (params.step == "mapping" || params.step == "filtering" || params.step == "calling" || params.step == "rna")){

      //******************
      // PAIRED BAMS
      chDesignPaired=chDesign.paired.map{it -> [[it[0], it[1]], it]}
      chDesignTumorOnly=chDesign.tumorOnly.map{it -> [it[0], it]}
      chDesignGermlineOnly=chDesign.germlineOnly.map{it -> [it[1], it]}

      chPairBam = chProcBam
        .combine(chProcBam)
        .map{meta1, bam1, bai1, meta2, bam2, bai2 ->
          [[meta1.id, meta2.id], [meta1, meta2, bam1, bai1, bam2, bai2]]}
        .combine(chDesignPaired, by:0)
        .map{it ->
          def meta = [tumor_id:it[1][0].id, normal_id:it[1][1].id, pair_id:it[2][2], id:it[0][0]+"_vs_"+it[0][1], status:"pair", sex:it[2][3]]
          [meta, it[1][2], it[1][3], it[1][4], it[1][5]]
        }

      //[meta], tumor_bam, tumor_bai
      chTumorBam = chPairBam
        .map{ it ->
          def meta = [tumor_id:it[0].tumor_id,id:it[0].tumor_id, status: "tumor", sex:it[0].sex]
          return [meta, it[1], it[2] ]
        }

      //[meta], normal_bam, normal_bai
      chNormalBam = chPairBam
        .map{ it ->
          def meta = [id:it[0].normal_id, status: "normal", sex:it[0].sex]
          return [meta, it[3], it[4] ]
        }
      chSingleBam = chNormalBam.mix(chTumorBam)

      //******************  
      // TUMOR ONLY 

      chTumorOnlyBam = chProcBam
        .map{ meta, bam, bai -> [meta.id, [meta, bam, bai]]}
        .combine(chDesignTumorOnly, by:0)
        .map{ it ->
          def meta = [tumor_id:it[0], id:it[0], status: "tumor", sex:it[2][3]]
          return [meta, it[1][1], it[1][2] ]
        }
      chSingleBam = chSingleBam.mix(chTumorOnlyBam) 

      //*******************
      // GERMLINE ONLY

      chGermlineOnlyBam = chProcBam
        .map{ meta, bam, bai -> [meta.id, [meta, bam, bai]]}
        .combine(chDesignGermlineOnly, by:0)
        .map{ it ->
          def meta = [id:it[0], status: "normal", sex:it[2][3]]
          return [meta, it[1][1], it[1][2] ]
        }
      chSingleBam = chSingleBam.mix(chGermlineOnlyBam)
    }

    chSingleBam = chSingleBam.unique()
    chAllGermlineVcf=Channel.empty()

    if('haplotypecaller' in tools){
    haplotypeCallerFlow(
        chSingleBam,
        intervals,
        dbsnp,
        dbsnpIndex,
        fasta,
        fai,
        dict
    )
    chVersions = chVersions.mix(haplotypeCallerFlow.out.versions)
    chAllGermlineVcf = haplotypeCallerFlow.out.vcf
    }

    chFilteredRnaVcf = Channel.empty()
    variantFiltration(
      chAllGermlineVcf,
      fasta,
      fai,
      dict
    )
    chVersions = chVersions.mix(variantFiltration.out.versions)

    emit:
    vcf = variantFiltration.out.vcf
}
