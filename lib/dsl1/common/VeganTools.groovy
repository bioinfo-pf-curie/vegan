/* groovylint-disable LineLength */
package common

import nextflow.Channel
import nextflow.Nextflow

/**
 * Any method specific to the actual pipeline are saved in this class as public methods
 */
abstract class VeganTools extends NFTools {

  /**
   * Has the run name been specified by the user?
   * This has the bonus effect of catching both -name and --name
   *
   * @return customRunName
   */
  def checkRunName(workflowRunName, runName) {
    return workflowRunName ==~ /[a-z]+_[a-z]+/ && runName ?
      runName : workflowRunName
  }

  /**
   * Change step to annotate if input is a vcf file
   *
   * @return step
   */
  def getStep(input, step) {
    return input && ["vcf", "vcf.gz"].collect { hasExtension(input, it) }.any() ? "annotate" : step ? step.toLowerCase() : ''
  }

  /**
   * If no input file specified, trying to get TSV files corresponding to step in the TSV directory (only for steps
   * recalibrate and variantCalling)
   *
   * @return inputPath
   */
  def getPath(step, inputParam, outDir) {
    if (!(inputParam && ["tsv", "csv", "vcf", "vcf.gz"].collect { hasExtension(inputParam, it) }.any())) {
      log.warn "Input file $inputParam extension is actually not supported (tsv, csv, vcf, vcf.gz)"
    }
    if (!inputParam && !['mapping', 'annotate'].contains(step)) {
      inputPath = step == 'recalibrate' ? "$outDir/Preprocessing/CSV/samplePlan.recal.csv" : "$outDir/Preprocessing/TSV/recalibrated.tsv"
    } else {
      inputPath = inputParam
    }
    return inputPath
  }


  /**
   * Parse first line of a FASTQ file, return the flowcell id and lane number.
   *
   * @param path
   * @return
   */
  def flowcellLaneFromFastq(path) {
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    InputStream fileStream = new FileInputStream(path.toFile())
    InputStream gzipStream = new java.util.zip.GZIPInputStream(fileStream)
    Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
    BufferedReader buffered = new BufferedReader(decoder)
    def line = buffered.readLine()
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(' ')[0].split(':')
    String fcid
    int lane
    if (fields.size() == 7) {
      // CASAVA 1.8+ format
      fcid = fields[2]
      lane = fields[3].toInteger()
    } else if (fields.size() == 5) {
      fcid = fields[0]
      lane = fields[1].toInteger()
    }
    [fcid, lane]
  }


  /**
   * Create a channel of germline FASTQs from a directory dirPath: "my_samples/*"
   * All FASTQ files in subdirectories are collected and emitted
   * they must have _R1_ and _R2_ in their names.
   *
   * @param dirPath
   * @return
   */
  // TODO: Use Channel.fromFilePairs
  def extractFastqFromDir(dirPath, singleEnd) {
    def fastq = Channel.create()
    // a temporary channel does all the work
    if (dirPath) {
      Channel
        .fromPath(dirPath, type: 'dir')
        .ifEmpty { error "No directories found matching dirPath '${dirPath}'" }
        .subscribe onNext: { 
          sampleDir ->
          // the last name of the sampleDir is assumed to be a unique sample id
          sampleID = sampleDir.getFileName().toString()

          for (path1 in file("${sampleDir}/**_R1_*.fastq.gz")) {
            assert path1.getName().contains('_R1_')
            path2 = file(path1.toString().replace('_R1_', '_R2_'))
            if (!path2.exists()) error "Path '${path2}' not found"
            (flowcell, lane) = flowcellLaneFromFastq(path1)
            sampleName = sampleID
            rgId = "${flowcell}.${sampleName}.${lane}"
            result = singleEnd ? [sampleID, sampleName, [path1]] : [sampleID, sampleName, [path1, path2]]
            fastq.bind(result)
          }
        }, onComplete: { fastq.close() }
    }
    return fastq
  }

  /**
   * Channeling the input file containing FASTQ or BAM
   * Format is: "idSample,sampleName,pathToFastq1,[pathToFastq2]"
   * or: "idSample,sampleName,pathToBam"
   *
   * @param inputPath
   * @return
   */
  def extractFastqOrBam(inputPath, sep, singleEnd, reads, readPaths) {
    if (inputPath) {
      return Channel
        .fromPath(inputPath)
        .splitCsv(sep: sep, header: false)
        .map { row ->
          def sampleID = row[0]
          def sampleName = row[1]
          def inputFile1 = returnFile(row[2])
          def inputFile2 = 'null'

          if ((!singleEnd) && (hasExtension(inputFile1, 'fastq.gz') || hasExtension(inputFile1, 'fq.gz') || hasExtension(inputFile1, 'fastq'))) {
            checkNumberOfItem(row, 4)
            inputFile2 = returnFile(row[3])
            if (!hasExtension(inputFile2, 'fastq.gz') && !hasExtension(inputFile2, 'fq.gz') && !hasExtension(inputFile2, 'fastq')) {
              Nextflow.exit(1, "File: ${inputFile2} has the wrong extension. See --help for more information")
            }
          } else if (hasExtension(inputFile1, 'bam')) {
            checkNumberOfItem(row, 3)
          } else {
            log.warn "No recognisable extention for input file: ${inputFile1}"
          }
          // ["sampleID": sampleID, "sampleName": sampleName, "inputFile1": inputFile1, "inputFile2": inputFile2]
          return singleEnd ? [sampleID, sampleName, [inputFile1]] : [sampleID, sampleName, [inputFile1, inputFile2]]
        }
    } else if (readPaths) {
      return Channel
        .fromList(readPaths)
        .map { row ->
          def sampleId = row[0]
          def inputFile1 = returnFile(row[1][0])
          def inputFile2 = singleEnd ? null: returnFile(row[1][1])
          singleEnd ? [sampleId, sampleId, [inputFile1]] : [sampleId, sampleId, [inputFile1, inputFile2]]
        }
        .ifEmpty { Nextflow.exit 1, "params.readPaths was empty - no input files supplied" }
    } else {
      return Channel
        .fromFilePairs(reads, size: singleEnd ? 1 : 2)
        .ifEmpty { Nextflow.exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .map { row -> singleEnd ? [row[0], row[0], [row[1][0]]] : [row[0], row[0], [row[1][0], row[1][1]]] }
    }
  }

  /**
   * Channeling the input file containing Filtered Bams
   * Format is: "sampleID sampleName vCType bam bai"
   *
   * @param inputPath
   * @return
   */
  def extractRecal(inputPath, sep = '\t') {
    return inputPath ? Channel
      .fromPath(inputPath)
      .splitCsv(sep: sep)
      .map { row ->
        checkNumberOfItem(row, 5)
        def sampleID = row[0]
        def sampleName = row[1]
        def vCType = row[2]
        def bamFile = returnFile(row[3])
        def baiFile = returnFile(row[4])

        if (!hasExtension(bamFile, "bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
        if (!hasExtension(baiFile, "bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"

        [sampleID, sampleName, vCType, bamFile, baiFile]
      } : Channel.empty()
  }

  /**
   * Channeling the input file containing BAM.
   * Format is: "subject gender status sample bam bai"
   *
   * @param inputPath
   * @return
   */
  def extractBam(inputPath, sep = '\t') {
    return inputPath ? Channel
      .fromPath(inputPath)
      .splitCsv(sep: sep)
      .map { row ->
        checkNumberOfItem(row, 5)
        def sampleId = row[0]
        def sampleName = row[1]
        def sampleType = row[2]
        def bamFile = returnFile(row[3])
        def baiFile = returnFile(row[4])

        if (!hasExtension(bamFile, "bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
        if (!hasExtension(baiFile, "bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"

        return [sampleId, sampleName, sampleType, bamFile, baiFile]
      } : Channel.empty()
  }

  def getDesign(Object designPath) {
    def designFile = designPath ? Nextflow.file(designPath) : null
    def designExt = designPath ? getExtension(designPath, ["tsv", "csv"]) : ""
    def separator = (designExt == 'tsv') ? '\t' : (designExt == 'csv') ? ',' : ''
    if (designExt) {
      return Channel.of(designFile)
        .splitCsv(sep: separator, header: ['germlineId', 'tumorId', 'pairId', 'sex']).map { row ->
        checkNumberOfItem(row, 4)
        [row.germlineId, row.tumorId, row.pairId, row.sex]
      }
    } else {
      return Channel.empty()
    }
  }

  /**
   * Extract information from inputPath
   *
   * @return samplePlanCh, samplePlanPathCh
   */
  def getSamplePlan(String inputPath, String step, Boolean singleEnd, String reads, List readPaths) {
    
    def samplePlanPathCh = inputPath ? Channel.fromPath(inputPath) : null

    def samplePlanCh
    def samplePlanPathTmpCh
    
    // if inputPath is a valid path or reads is a valid string
    if (inputPath || reads || readPaths) {
      def inputExt = getExtension(inputPath, ['csv', 'tsv'])
      // Define csv as the default separator
      def sep = (inputExt == 'tsv') ? '\t' : (inputExt == 'csv') ? ',' : ','
      switch (step) {
      // idSample,sampleName,pathToFastq1,[pathToFastq2]
      // idSample,sampleName,pathToBam
        case 'mapping':
          (samplePlanPathTmpCh, samplePlanCh) = extractFastqOrBam(inputPath, sep, singleEnd, reads, readPaths).into(2)
          
          if (!samplePlanPathCh) {
            samplePlanPathCh = samplePlanPathTmpCh
              .collectFile() {
                item -> singleEnd ? ["samplePlan.csv", item[0] + ',' + item[1] + ',' + item[2][0] + '\n'] : ["samplePlan.csv", item[0] + ',' + item[1] + ',' + item[2][0] + ',' + item[2][1] + '\n']
              }
          }
          break
      // idSample,sampleName,pathToBam,pathToBai,pathToRecalTable
      // [sampleID, sampleName, bamFile, baiFile, recalTable]
        case 'recalibrate':
          samplePlanCh = extractRecal(inputPath, sep)
          break
      // idSample,sampleName,pathToBam,pathToBai
      // [sampleID, sampleName, bamFile, baiFile]
        case 'variantcalling':
          samplePlanCh = extractBam(inputPath, sep)
          break
        case 'annotate':
          break
        default:
          Nextflow.exit(1, "Unknown step ${step}")
      }
    } else if (inputPath && Nextflow.file(inputPath).isDirectory()) {
      if (step != 'mapping') {
        Nextflow.exit(1, 'No other step than "mapping" support a dir as an input')
      }
      log.info "Reading $inputPath directory"
      def fastqTMP
      (samplePlanCh, fastqTMP) = extractFastqFromDir(inputPath, singleEnd).into(2)
      fastqTMP
        .toList()
        .tap{ samplePlanPathTmpCh }
        .subscribe onNext: {
          if (it.size() == 0) {
            Nextflow.exit(1, "No FASTQ files found in --input directory '${params.input}'")
          }
        }
      samplePlanPathCh = samplePlanPathTmpCh
        .collectFile(){
          item -> singlEnd ? ["samplePlan.csv", item[0] + ',' + item[1] + ',' + item[2][0] + '\n'] : ["samplePlan.csv", item[0] + ',' + item[1] + ',' + item[2][0] + ',' + item[2][1] + '\n']
        }
    } else if (inputPath && step == 'annotate') {
      log.info "Annotating ${inputPath}"
    } else if (step == 'annotate') {
      log.info "Trying automatic annotation on file in the VariantCalling directory"
    } else if (reads) {
      // TODO: probably useless since it's already done with extractFastqOrBam function (same for readsPaths)
      (samplePlanCh, samplePlanPathTmpCh) = Channel
        .fromFilePairs(reads, size: singleEnd ? 1 : 2)
        .ifEmpty {
          Nextflow.exit(1, """\
Cannot find any reads matching: ${reads}
NB: Path needs to be enclosed in quotes!
NB: Path requires at least one * wildcard!
If this is single-end data, please specify --singleEnd on the command line.""");
        }
        .map { row -> singleEnd ? [row[0], row[0], [row[1][0]]] : [row[0], row[0], [row[1][0], row[1][1]]] }.into(2)
      samplePlanPathCh = samplePlanPathTmpCh
        .collectFile(){
          item -> singlEnd ? ["samplePlan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + '\n'] : ["samplePlan.csv", item[0] + ',' + item[0] + ',' + item[1][0] + ',' + item[1][1] + '\n']
        }
    } else {
      Nextflow.exit(1, 'No input data were defined, see --help')
    }
    
    return samplePlanCh && samplePlanPathCh ? [samplePlanCh, samplePlanPathCh]: samplePlanCh ? [samplePlanCh, channel.empty()] : samplePlanPathCh ? [Channel.empty(), samplePlanPathCh] : [Channel.empty(), Channel.empty()]
  }

  /**
   * Extract gender and status from design Channel
   *
   * @param channel
   * @return
   */
  def extractInfos(designChannel) {
    def genderMap = [:]
    def statusMap = [:]
    def pairMap = [:]
    designChannel.map { it ->
      def normalSampleID = it[0]
      def tumorSampleID = it[1]
      def pairName = it[2]
      def gender = it[3]
      genderMap[normalSampleID] = gender
      genderMap[tumorSampleID] = gender
      statusMap[normalSampleID] = 0
      statusMap[tumorSampleID] = 1
      pairMap[[normalSampleID, tumorSampleID]] = pairName
    }
    return [genderMap, statusMap, pairMap]
  }

  // TODO: fork method not working
  /**
   * Fork sample plan in bam and pair channels for mapping and add an id for multiplexed files
   *
   * @param mappingSamplePlan nextflow channel
   * @return [bamSamplesCh, pairReadsSamplesCh]
   */
  def forkMappingSamplePlan(Object mappingSamplePlanCh) {
    def runIds = [:]
    def samplePlanForks = mappingSamplePlanCh.map {
      runIds[it[0]] = runIds.containsKey(it[0]) ? runIds[it[0]] + 1 : 0
      return it[0, 1] + [[it[0], runIds[it[0]].toString()].join("_")] + it[2..-1]
    }.branch {
      bam: it[3] =~ /.*bam$/
      pair: it[3] =~ /.*(fastq.gz|fq.gz|fastq|fq)$/
    }
    return [samplePlanForks.bamCh, samplePlanForks.pairCh]
  }

  /**
   * Return sample status (0 == Normal, 1 == Tumor)
   *
   * @return Integer
   */
  def returnStatus(it) {
    if (!(it in [0, 1])) exit 1, "Status is not recognized in TSV file: ${it}, see --help for more information"
    return it
  }

  /**
   * Remove .ann .gz and .vcf extension from a VCF file
   *
   * @param file
   * @return
   */
  def reduceVCF(file) {
    return file.fileName.toString().minus(".ann").minus(".vcf").minus(".gz")
  }

}
