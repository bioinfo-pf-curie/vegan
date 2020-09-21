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
    static String checkRunName(workflowRunName, runName) {
        return workflowRunName ==~ /[a-z]+_[a-z]+/ && runName ?
                runName : workflowRunName
    }

    /**
     * Change step to annotate if input is a vcf file
     *
     * @return step
     */
    static String getStep(input, step) {
        return input && ["vcf", "vcf.gz"].collect{ hasExtension(input, it) }.any() ? "annotate": step ? step.toLowerCase() : ''
    }

    /**
     * If no input file specified, trying to get TSV files corresponding to step in the TSV directory (only for steps
     * recalibrate and variantCalling)
     *
     * @return inputPath
     */
    static Object getPath(step, inputParam, outputDir) {
        def inputPath = inputParam && ["tsv", "csv", "vcf", "vcf.gz"].collect { hasExtension(inputParam, it) }.any() ?
                inputParam : null
        if (!inputParam && !['mapping', 'annotate'].contains(step)) {
            inputPath = step == 'recalibrate' ? "$outputDir/Preprocessing/TSV/duplicateMarked.tsv" : "$outputDir/Preprocessing/TSV/recalibrated.tsv"
        }
        return inputPath
    }

    /**
     * Channeling the input file containing BAM.
     * Format is: "subject gender status sample bam bai"
     *
     * @param inputFile
     * @return
     */
    def extractBam(inputFile, sep='\t') {
        Channel.of(inputFile)
                .splitCsv(sep: sep)
                .map { row ->
                    checkNumberOfItem(row, 4)
                    def platformID = row[0]
                    def sampleName   = row[1]
                    def bamFile   = returnFile(row[2])
                    def baiFile   = returnFile(row[3])

                    if (!hasExtension(bamFile, "bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
                    if (!hasExtension(baiFile, "bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"

//                    return [idPatient, gender, status, idSample, bamFile, baiFile]
                    return [platformID, sampleName, bamFile, baiFile]
                }
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
     * Create a channel of germline FASTQs from a directory pattern: "my_samples/*"
     * All FASTQ files in subdirectories are collected and emitted
     * they must have _R1_ and _R2_ in their names.
     *
     * @param pattern
     * @return
     */
    def extractFastqFromDir(pattern) {
        def fastq = Channel.create()
        // a temporary channel does all the work
        Channel
                .fromPath(pattern, type: 'dir')
                .ifEmpty { error "No directories found matching pattern '${pattern}'" }
                .subscribe onNext: { sampleDir ->
            // the last name of the sampleDir is assumed to be a unique sample id
            sampleID = sampleDir.getFileName().toString()

            for (path1 in file("${sampleDir}/**_R1_*.fastq.gz")) {
                assert path1.getName().contains('_R1_')
                path2 = file(path1.toString().replace('_R1_', '_R2_'))
                if (!path2.exists()) error "Path '${path2}' not found"
                (flowcell, lane) = flowcellLaneFromFastq(path1)
                sampleName = sampleID
                rgId = "${flowcell}.${sampleName}.${lane}"
                result = [sampleID, sampleName, path1, path2]
                fastq.bind(result)
            }
        }, onComplete: { fastq.close() }
        fastq
    }

    /**
     * Channeling the input file containing FASTQ or BAM
     * Format is: "idSample,sampleName,pathToFastq1,[pathToFastq2]"
     * or: "idSample,sampleName,pathToBam"
     *
     * @param inputFile
     * @return
     */
    static def extractFastqOrBam(inputFile, sep) {
        return Channel.of(inputFile)
                .splitCsv(sep: sep)
                .map { row ->
                    def sampleID = row[0]
                    def sampleName  = row[1]
                    def inputFile1  = returnFile(row[2])
                    def inputFile2  = "null"
                    if (hasExtension(inputFile1, "fastq.gz") || hasExtension(inputFile1, "fq.gz")) {
                        checkNumberOfItem(row, 4)
                        inputFile2 = returnFile(row[3])
                        if (!hasExtension(inputFile2, "fastq.gz") && !hasExtension(inputFile2, "fq.gz")) {
                            Nextflow.exit 1, "File: ${inputFile2} has the wrong extension. See --help for more information"
                        }
                    }
                    else if (hasExtension(inputFile1, "bam")) checkNumberOfItem(row, 3)
                    else "No recognisable extention for input file: ${inputFile1}"
                    // ["sampleID": sampleID, "sampleName": sampleName, "inputFile1": inputFile1, "inputFile2": inputFile2]
                    [sampleID, sampleName, inputFile1, inputFile2]
                }
    }

    /**
     * Channeling the input file containing Recalibration Tables
     * Format is: "subject gender status sample bam bai recalTables"
     *
     * @param inputFile
     * @return
     */
    def extractRecal(inputFile, sep='\t') {
        Channel.of(inputFile)
                .splitCsv(sep: sep)
                .map { row ->
                    checkNumberOfItem(row, 5)
                    def sampleID = row[0]
                    def sampleName   = row[1]
                    def bamFile    = returnFile(row[2])
                    def baiFile    = returnFile(row[3])
                    def recalTable = returnFile(row[4])

                    if (!hasExtension(bamFile, "bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
                    if (!hasExtension(baiFile, "bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"
                    if (!hasExtension(recalTable, "recal.table")) exit 1, "File: ${recalTable} has the wrong extension. See --help for more information"

                    [sampleID, sampleName, bamFile, baiFile, recalTable]
                }
    }

    static def getDesign(Object designPath) {
        def designFile = designPath ? Nextflow.file(designPath) : null
        def designExt = designPath ? getExtension(designPath, ["tsv", "csv"]) : ""
        def separator = (designExt == 'tsv') ? '\t' : (designExt == 'csv') ? ',' : ''
        if (designExt) {
            return Channel.of(designFile)
                .splitCsv(sep: separator).map{row ->
                    checkNumberOfItem(row, 4)
                    def normalSampleID = row[0]
                    def tumorSampleID  = row[1]
                    def pairName       = row[2]
                    def sex            = row[3]
                    [normalSampleID, tumorSampleID, pairName, sex]
            }
        }
    }

    /**
     * Extract information from inputPath
     *
     * @return inputSample
     */
    def getSamplePlan(Object inputPath) {
        inputSample = Channel.empty()
        def input = inputPath ? Nextflow.file(inputPath) : null
        if (inputPath) {
            inputExt = getExtension(inputPath, ['csv', 'tsv'])
            sep = (inputExt == 'tsv') ? '\t' : (inputExt == 'csv') ? ',' : ''
            switch (step) {
                // idSample,sampleName,pathToFastq1,[pathToFastq2]
                // idSample,sampleName,pathToBam
                case 'mapping': return extractFastqOrBam(input, sep); break
                // idSample,sampleName,pathToBam,pathToBai,pathToRecalTable
                // [sampleID, sampleName, bamFile, baiFile, recalTable]
                case 'recalibrate': return extractRecal(input, sep); break
                // idSample,sampleName,pathToBam,pathToBai
                // [sampleID, sampleName, bamFile, baiFile]
                case 'variantcalling': return extractBam(input, sep); break
                case 'annotate': break
                default: Nextflow.exit 1, "Unknown step ${step}"
            }
        } else if (inputPath && input.isDirectory()) {
            if (step != 'mapping') Nextflow.exit 1, 'No other step than "mapping" support a dir as an input'
            log.info "Reading $inputPath directory"
            inputSample = extractFastqFromDir(inputPath)
            (inputSample, fastqTMP) = inputSample.into(2)
            fastqTMP.toList().subscribe onNext: {
                if (it.size() == 0) Nextflow.exit 1, "No FASTQ files found in --input directory '${params.input}'"
            }
        } else if (inputPath && step == 'annotate') {
            log.info "Annotating ${inputPath}"
        }
        else if (step == 'annotate') {
            log.info "Trying automatic annotation on file in the VariantCalling directory"
        } else {
            Nextflow.exit 1, 'No sample were defined, see --help'
        }
        return inputSample
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
        designChannel.map{ it ->
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
        mappingSamplePlanCh.map {
            runIds[it[0]] = runIds.containsKey(it[0]) ? runIds[it[0]] + 1 : 0
            return it[0, 1] + [[it[0], runIds[it[0]].toString()].join("_")] + it[2..-1]
        }.branch {
            bam: it[3] =~ /.*bam$/
            pair: it[3] =~ /.*(fastq.gz|fq.gz|fastq|fq)$/
        }.set { samplePlanForks }
        return [samplePlanForks.bamCh, samplePlanForks.pairCh]
    }

    /**
     * Return sample status (0 == Normal, 1 == Tumor)
     *
     * @return Integer
     */
    Integer returnStatus(it) {
        if (!(it in [0, 1])) exit 1, "Status is not recognized in TSV file: ${it}, see --help for more information"
        return it
    }

    /**
     * Remove .ann .gz and .vcf extension from a VCF file
     *
     * @param file
     * @return
     */
    String reduceVCF(file) {
        return file.fileName.toString().minus(".ann").minus(".vcf").minus(".gz")
    }

}
