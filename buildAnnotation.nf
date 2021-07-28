#!/usr/bin/env nextflow
/*
Copyright Institut Curie 2020
This software is a computer program whose purpose is to analyze high-throughput sequencing data.
You can use, modify and/ or redistribute the software under the terms of license (see the LICENSE file for more details).
The software is distributed in the hope that it will be useful, but "AS IS" WITHOUT ANY WARRANTY OF ANY KIND.
Users are therefore encouraged to test the software's suitability as regards their requirements in conditions enabling the security of their systems and/or data.
The fact that you are presently reading this means that you have had knowledge of the license and that you accept its terms.

This script is based on the nf-core guidelines. See https://nf-co.re/ for more information
*/

import common.VeganTools
@BaseScript(VeganTools)
import groovy.transform.BaseScript

/*
================================================================================
                        project : EUCANCAN/vegan
================================================================================
 Build VEGAN annotation files
--------------------------------------------------------------------------------

nextflow run buildAnnotation.nf --genome 'hg19' -profile conda --genomeAnnotationPath '/data/annotations/pipelines/'
*/

welcome()

/*
 * Genome configuration
 */

params << [
  fasta: params.genome ? params.genomes[params.genome].fasta ?: null : null,
  dbsnp: params.genome ? params.genomes[params.genome].dbsnp ?: null : null,
  germlineResource: params.genome ? params.genomes[params.genome].germlineResource ?: null : null,
  knownIndels: params.genome ? params.genomes[params.genome].knownIndels ?: null : null
]

/*
 * Channel
 */

fastaCh = params.fasta ? Channel.value(file(params.fasta)) : "null"
dbsnpCh = params.dbsnp ? Channel.value(file(params.dbsnp)) : "null"
germlineResourceCh = params.germlineResource ? Channel.value(file(params.germlineResource)) : "null"
knownIndelsCh = params.knownIndels ? Channel.value(file(params.knownIndels)) : "null"
ponCh = params.pon ? Channel.value(file(params.pon)) : "null"

/*
================================================================================
                                BUILDING INDEXES
================================================================================
*/

/*
process buildBWAindexes {
  label 'bwa'
  label 'highCpu'
  label 'highMem'

  publishDir "${params.outDir}/bwaIndex", mode: params.publishDirMode
 
  input:
  file(fasta) from fastaCh

  output:
  file("${fasta}.*") into bwaIndexesCh
  file("v_bwa.txt") into bwaVersionCh

  script:
  """
  bwa index ${fasta}
  bwa &> v_bwa.txt 2>&1 || true
  """
}
*/

process buildDict {
  label 'gatk'
  label 'minCpu'
  label 'lowMem'

  publishDir "${params.outDir}/genome/", mode: params.publishDirMode

  input:
  file(fasta) from fastaCh

  output:
  file("${fasta.baseName}.dict") into dictBuiltCh

  script:
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g" \
      CreateSequenceDictionary \
      --REFERENCE ${fasta} \
      --OUTPUT ${fasta.baseName}.dict
  """
}

process buildFastaFai {
  label 'samtools'
  label 'minCpu'
  label 'lowMem'

  tag "${fasta}"

  publishDir "${params.outDir}/genome/", mode: params.publishDirMode

  input:
  file(fasta) from fastaCh

  output:
  file("${fasta}.fai") into fastaFaiCh, fastaFaiVcfCh

  script:
  """
  samtools faidx ${fasta}
  """
}

process buildDbsnpIndex {
  label 'tabix'
  label 'minCpu'
  label 'lowMem'

  tag "${dbsnp}"

  publishDir "${params.outDir}/indexes/", mode: params.publishDirMode

  input:
  file(dbsnp) from dbsnpCh

  output:
  file("${dbsnp}.tbi") into dbsnpIndexBuiltCh

  script:
  """
  tabix -p vcf ${dbsnp}
  """
}

process updateVcfName {
  label 'bcftools'
  label 'minCpu'
  label 'lowMem'

  input:
  file(dbsnp) from dbsnpCh
  file(fai) from fastaFaiVcfCh

  output:
  file("*base.vcf.gz") into desnpBaseCh

  script:
  grep -v "_" ${fai} | awk '{x=$1; gsub("chr","",$1); print x"\t"$1}' > chrom.map
  bcftools annotate --rename-chrs chrom.map -O z ${dbsnp} > ${dbsnp.baseName}_rename.vcf.gz
  bcftools reheader -f ${fai} ${dbsnp.baseName}_rename.vcf.gz > ${dbsnp.baseName}_base.vcf.gz
}


process buildGermlineResourceIndex {
  label 'tabix'
  label 'minCpu'
  label 'lowMem'

  publishDir "${params.outDir}/indexes/", mode: params.publishDirMode

  input:
  file(germlineResource) from germlineResourceCh

  output:
  file("${germlineResource}.tbi") into germlineResourceIndexBuiltCh

  script:
  """
  tabix -p vcf ${germlineResource}
  """
}

process buildKnownIndelsIndex {
  label 'tabix'
  label 'minCpu'
  label 'lowMem'

  publishDir "${params.outDir}/indexes/", mode: params.publishDirMode

  input:
  file(knownIndels) from knownIndelsCh.flatten()

  output:
  file("${knownIndels}.tbi") into knownIndelsIndexBuiltCh

  script:
  """
  tabix -p vcf ${knownIndels}
  """
}

process buildPonIndex {
  label 'tabix'
  label 'minCpu'
  label 'lowMem'

  publishDir "${params.outDir}/pon", mode: params.publishDirMode

  when:
  params.pon

  input:
  file(pon) from ponCh

  output:
  file("${pon}.tbi") into ponIndexBuiltCh

  script:
  """
  tabix -p vcf ${pon}
  """
}

process buildIntervals {
  label 'onlyLinux'
  label 'minCpu'
  label 'lowMem'

  tag "${fastaFai}"

  publishDir params.outDir, mode: params.publishDirMode,
    saveAs: {params.saveGenomeIndex ? "reference_genome/${it}" : null }

  input:
  file(fastaFai) from fastaFaiCh

  output:
  file("${fastaFai.baseName}.bed") into intervalBuiltCh

  script:
  """
  awk -v FS='\t' -v OFS='\t' '{ print \$1, \"0\", \$2 }' ${fastaFai} > ${fastaFai.baseName}.bed
  """
}

