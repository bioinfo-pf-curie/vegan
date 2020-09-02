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
    String getRunName() {
        return workflow.runName ==~ /[a-z]+_[a-z]+/ && params.runName ?
                params.runName : workflow.runName
    }

    /**
     * Change step to annotate if input is a vcf file
     *
     * @return step
     */
    String getStep() {
        return params.input && ["vcf", "vcf.gz"].collect{ hasExtension(params.input, it) }.any() ? "annotate": params.step ? params.step.toLowerCase() : ''
    }

    /**
     * If no input file specified, trying to get TSV files corresponding to step in the TSV directory (only for steps
     * recalibrate and variantCalling)
     *
     * @return tsvPath
     */
    Object getTsvPath() {
        def tsvPath = params.input && ["tsv", "vcf", "vcf.gz"].collect { hasExtension(params.input, it) }.any() ?
                params.input : null
        def step = getStep()
        if (!params.input && step != 'mapping' && step != 'annotate') {
            tsvPath = step == 'recalibrate' ? "${params.outputDir}/Preprocessing/TSV/duplicateMarked.tsv" : "${params.outputDir}/Preprocessing/TSV/recalibrated.tsv"
        }
        return tsvPath
    }

    /**
     * Channeling the TSV file containing BAM.
     * Format is: "subject gender status sample bam bai"
     *
     * @param tsvFile
     * @return
     */
    def extractBam(tsvFile) {
        Channel.of(tsvFile)
                .splitCsv(sep: '\t')
                .map { row ->
                    checkNumberOfItem(row, 6)
                    def idPatient = row[0]
                    def gender    = row[1]
                    def status    = returnStatus(row[2].toInteger())
                    def idSample  = row[3]
                    def bamFile   = returnFile(row[4])
                    def baiFile   = returnFile(row[5])

                    if (!hasExtension(bamFile, "bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
                    if (!hasExtension(baiFile, "bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"

                    return [idPatient, gender, status, idSample, bamFile, baiFile]
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
            sampleId = sampleDir.getFileName().toString()

            for (path1 in file("${sampleDir}/**_R1_*.fastq.gz")) {
                assert path1.getName().contains('_R1_')
                path2 = file(path1.toString().replace('_R1_', '_R2_'))
                if (!path2.exists()) error "Path '${path2}' not found"
                (flowcell, lane) = flowcellLaneFromFastq(path1)
                patient = sampleId
                gender = 'ZZ'  // unused
                status = 0  // normal (not tumor)
                rgId = "${flowcell}.${sampleId}.${lane}"
                result = [patient, gender, status, sampleId, rgId, path1, path2]
                fastq.bind(result)
            }
        }, onComplete: { fastq.close() }
        fastq
    }

    /**
     * Channeling the TSV file containing FASTQ or BAM
     * Format is: "subject gender status sample lane fastq1 fastq2"
     * or: "subject gender status sample lane bam"
     *
     * @param tsvFile
     * @return
     */
    def extractFastqOrBam(tsvFile) {
        Channel.from(tsvFile)
                .splitCsv(sep: '\t')
                .map { row ->
                    def idPatient  = row[0]
                    def gender     = row[1]
                    def status     = returnStatus(row[2].toInteger())
                    def idSample   = row[3]
                    def idRun      = row[4]
                    def file1      = returnFile(row[5])
                    def file2      = "null"
                    if (hasExtension(file1, "fastq.gz") || hasExtension(file1, "fq.gz")) {
                        checkNumberOfItem(row, 7)
                        file2 = returnFile(row[6])
                        if (!hasExtension(file2, "fastq.gz") && !hasExtension(file2, "fq.gz")) exit 1, "File: ${file2} has the wrong extension. See --help for more information"
                    }
                    else if (hasExtension(file1, "bam")) checkNumberOfItem(row, 6)
                    else "No recognisable extention for input file: ${file1}"

                    [idPatient, gender, status, idSample, idRun, file1, file2]
                }
    }

    /**
     * Channeling the TSV file containing Recalibration Tables
     * Format is: "subject gender status sample bam bai recalTables"
     *
     * @param tsvFile
     * @return
     */
    def extractRecal(tsvFile) {
        Channel.from(tsvFile)
                .splitCsv(sep: '\t')
                .map { row ->
                    checkNumberOfItem(row, 7)
                    def idPatient  = row[0]
                    def gender     = row[1]
                    def status     = returnStatus(row[2].toInteger())
                    def idSample   = row[3]
                    def bamFile    = returnFile(row[4])
                    def baiFile    = returnFile(row[5])
                    def recalTable = returnFile(row[6])

                    if (!hasExtension(bamFile, "bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
                    if (!hasExtension(baiFile, "bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"
                    if (!hasExtension(recalTable, "recal.table")) exit 1, "File: ${recalTable} has the wrong extension. See --help for more information"

                    [idPatient, gender, status, idSample, bamFile, baiFile, recalTable]
                }
    }

    /**
     * Extract information from tsvPath
     *
     * @return inputSample
     */
    def getInputSample(Object inputPath) {
        inputSample = Channel.empty()
        def inputFile = inputPath ? new File(inputPath) : null
//        if (inputPath && inputFile.exists()) {
        if (inputPath) {
            inputFile = Nextflow.file(inputPath)
            switch (step) {
                case 'mapping': inputSample = extractFastqOrBam(inputFile); break
                case 'recalibrate': inputSample = extractRecal(inputFile); break
                case 'variantcalling': inputSample = extractBam(inputFile); break
                case 'annotate': break
                default: Nextflow.exit 1, "Unknown step ${step}"
            }
        } else if (params.input && !hasExtension(params.input, "tsv")) {
            log.info "No TSV file"
            if (step != 'mapping') Nextflow.exit 1, 'No other step than "mapping" support a dir as an input'
            log.info "Reading ${params.input} directory"
            inputSample = extractFastqFromDir(params.input)
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
     * Extract gender and status from Channel
     *
     * @param channel
     * @return
     */
    def extractInfos(channel) {
        def genderMap = [:]
        def statusMap = [:]
        channel = channel.map{ it ->
            def idPatient = it[0]
            def gender = it[1]
            def status = it[2]
            def idSample = it[3]
            genderMap[idPatient] = gender
            statusMap[idPatient, idSample] = status
            [idPatient] + it[3..-1]
        }
        return [genderMap, statusMap, channel]
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
