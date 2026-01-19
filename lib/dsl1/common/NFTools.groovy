package common

import colorlog.MonoChrome
import colorlog.PolyChrome
import groovy.text.GStringTemplateEngine
import nextflow.script.BaseScript
import nextflow.script.WorkflowMetadata
import utils.ParamsLinter
import utils.ParamsReader
import org.slf4j.LoggerFactory
import org.slf4j.Logger
import nextflow.Nextflow
import nextflow.Channel


// TODO: check updates on Nextflow parameter scheme
// CF https://github.com/nextflow-io/nextflow/issues/866
// CF https://github.com/nf-core/tools/issues/577
// CF https://github.com/nf-core/tools/blob/dev/nf_core/pipeline-template/%7B%7Bcookiecutter.name_noslash%7D%7D/nextflow_schema.json
abstract class NFTools extends BaseScript {

    static final Logger log = LoggerFactory.getLogger(Nextflow.class)

    private static LinkedHashMap generateLogColors(Boolean monochromeLogs = false) {
         Map colors = monochromeLogs ? new MonoChrome().palette() : new PolyChrome().palette()
         return colors as LinkedHashMap
     }

    private static Map formatParameterHelpData(params) {

        Map result = [name: params.get('name'), value: '', usage: params.get('usage')]
        // value describes the expected input for the param
        result.value =  [params.type == boolean.toString() ? '' : params.type.toString().toUpperCase(), params.choices ?: ''].join(' ')
        return result
    }

    private static String prettyFormatParamGroupWithPaddingAndIndent(List paramGroup, String groupName, Integer padding = 2, Integer indent = 4) {

        def maxParamNameLength = paramGroup.collect { it.name.size() }.max()
        def paramChoices = paramGroup.findAll { it.choices }.collect { it.choices }
        def maxChoiceStringLength = paramChoices.collect { it.toString().size() }.max() ?: 0
        def maxTypeLength = paramGroup.collect { (it.type as String).size() }.max() ?: 0
        def maxValueLength = maxChoiceStringLength + maxTypeLength
        def usagePadding = indent + maxParamNameLength + maxValueLength + 2 * padding + 2
        Integer usageLength = ((150 - usagePadding) / 100).round(1) * 100

        def paramsFormattedList = paramGroup.sort { it.name }.collect {
            Map param ->
                def paramHelpData = formatParameterHelpData(param)
                def usage = []

                "${paramHelpData.usage}".eachMatch(".{1,${usageLength - 1}}([\\s\\.]|\$)", {
                    usage << it.pop()
                })
                usage = usage.join("\n" + " " * usagePadding)
                sprintf("%${indent}s%-${maxParamNameLength + padding}s %-${maxValueLength + padding}s %s\n", "", "--${paramHelpData.name}", "${paramHelpData.value}", "${usage}")

        }
        String.format("%s:\n%s", groupName.toUpperCase(), paramsFormattedList.join()).stripIndent()
    }

    private static String prettyFormatParamsWithPaddingAndIndent(List paramsWithUsage, Integer padding = 2, Integer indent = 4) {

        def groupedParamsWithUsage = paramsWithUsage.groupBy { it.group }
        def formattedParamsGroups = groupedParamsWithUsage.collect {
            prettyFormatParamGroupWithPaddingAndIndent(it.value, it.key, padding, indent)
        }
        return formattedParamsGroups.join('\n')
    }

    static Object readParamsFromJsonSettings(String path) {

        def paramsWithUsage
        try {
            paramsWithUsage = ParamsReader.readParamsFromJsonSettings(path).get("parameters")
        } catch (Exception e) {
            log.warn "Could not read parameters settings from JSON. $e"
            paramsWithUsage = Collections.emptyMap()
        }
        return paramsWithUsage
    }

    /**
     * Check if path has the correct file extension
     * @param it
     * @param extension
     * @return
     */
    static boolean hasExtension(it, String extension) {
        it.toString().toLowerCase().endsWith(extension.toLowerCase())
    }

    /**
     * Generate Help message
     * @param paramsWithUsage
     * @param workflow
     * @return
     */
    static String helpMessage(paramsWithUsage, workflow) {
        def CLIHelpMsg = []
        paramsWithUsage.each {it.group == "Mandatory arguments" ? CLIHelpMsg << "--" + it.name << it.type.toUpperCase() : ""}
        log.info String.format("""\
            
            Usage:
            
            The typical command for running the pipeline is as follows:
            
            nextflow run main.nf ${CLIHelpMsg.join(" ")}

            %s
            """.stripIndent(), prettyFormatParamsWithPaddingAndIndent(paramsWithUsage, 2, 4))
    }

    /**
     * Check if file extension from input path is defined in allowed extensions list
     * @param path
     * @param extensions
     * @return
     */
    String getExtension(path, List extensions) {

        for (extension in extensions) {
            if (hasExtension(path, extension)) {
                return extension
            }
        }

        def colors = binding.getVariable("colors")
        log.warn """\
        Can't guess field separator with the actual input file extension. Trying with the default one (.csv)
        """.stripIndent().toString()
        return ""
    }

    /**
     * Check if a row has the expected number of item
     *
     * @param row
     * @param number
     * @return Boolean
     */
    boolean checkNumberOfItem(row, Integer number) {
        def colors = binding.getVariable("colors")
        if (row.size() != number) {
            def message = """\
                ${colors.red}[WARNING] Malformed row in input file: ${row}
                Input file should have ${number} but have ${row.size()} items. see --help for more information${colors.reset}
                """.stripIndent()
            log.info message.toString()
        }
        return true
    }

    /**
     * Return file if it exists
     *
     * @return Nextflow.file Object
     */
    def returnFile(String it) {
        def colors = binding.getVariable("colors")
        if (it =~ /(http|ftp)/) {
            return Nextflow.file(it)
        } else if (!Nextflow.file(it).exists()) {
            Nextflow.exit(
              "${colors.red}[WARNING] Input file does not exists: ${it}, see --help for more information${colors.reset}"
            )
        } else {
            return Nextflow.file(it)
        }
    }

    /**
     * Print Header in assets in the log
     * @param params
     * @param workflow
     * @return
     */
    def nfHeader(params, WorkflowMetadata workflow) {

        LinkedHashMap context = binding.getVariable('colors')
        context << workflow.properties

        def engine = new GStringTemplateEngine()
        def txtTemplate = engine.createTemplate(
          new File("${workflow.projectDir}/assets/nfHeader.txt")
        ).make(context)
        log.info txtTemplate.toString()
    }

    /**
     * Print summary dict with the logger and return formatted channel
     * @param summary
     * @param workflow
     * @return
     */
    def summarize(summary) {
        def colors = binding.getVariable("colors")
        def workflow = binding.getVariable("workflow")
        log.info summary.collect { k, v -> "${k.padRight(18)}: $v" }.join("\n") +
                "\n${colors.dim}" + "-"*80 + "${colors.reset}"
        return Channel.fromList(summary.collect{ [it.key, it.value] })
          .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
                .reduce { a, b -> return [a, b].join("\n            ") }
                .map { x -> """
    id: 'summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'Workflow Summary'
    section_href: '$workflow.manifest.homePage'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    }

    /**
     * Check if params.hostnames is valid
     * @param params
     * @param workflow
     * @return
     */
    def checkHostname(params, WorkflowMetadata workflow) {
        def colors = binding.getVariable("colors")
        //def colors = binding.getVariable("colors")
        if (params.hostnames) {
            def hostname = "hostname".execute().text.trim()
            params.hostnames.each { prof, hosts ->
                hosts.each { host ->
                    if (hostname.contains(host) && !workflow.profile.contains(prof)) {
                        Nextflow.exit(1, """\
============================================================
  ${colors.red}WARNING!${colors.reset} You are running with `-profile $workflow.profile`
  but your machine hostname is ${colors.white}'$hostname'${colors.reset}
  ${colors.yellowBold}It's highly recommended that you use `-profile $prof${colors.reset}
============================================================""")
                    }
                }
            }
        }
}

    /**
     * Generate reports at the end of the workflow
     * @param workflow
     * @param params
     * @param reportFields
     * @param multiQCOutCh
     * @return
     */
    def makeReports(workflow, params, reportFields, multiQCOutCh) {

        def colors = binding.getVariable("colors")

        def engine = new groovy.text.GStringTemplateEngine()

        // Render the TXT template
        def txtTemplate = engine.createTemplate(new File("$workflow.projectDir/assets/onCompleteTemplate.txt")).make(reportFields)
        def txtReport = txtTemplate.toString()

        // Render the HTML template
        def hf = new File("$workflow.projectDir/assets/onCompleteTemplate.html")
        def htmlTemplate = engine.createTemplate(hf).make(reportFields)
        def htmlReport = htmlTemplate.toString()

        // Write summary e-mail HTML to a file
        def outDir = new File("${params.outDir}")
        if (!outDir.exists()) outDir.mkdirs()
        def output_hf = new File(outDir, "pipelineReport.html")
        output_hf.withWriter { w -> w << htmlReport }
        def output_tf = new File(outDir, "pipelineReport.txt")
        output_tf.withWriter { w -> w << txtReport }

        // On success try attach the multiqc report
        def mqcReport = null
        try {
            if (workflow.success) {
                mqcReport = multiQCOutCh.getVal()
                if (mqcReport.getClass() == ArrayList) {
                    log.warn "[$workflow.manifest.name] Found multiple reports from process 'multiqc', will use only one"
                    mqcReport = mqcReport[0]
                }
            }
        } catch (all) {
            log.warn "[$workflow.manifest.name] Could not attach MultiQC report to summary email"
        }

        if (params.email) {
            // Set up the e-mail variables
            def subject = workflow.success? "[$workflow.manifest.name] Successful: ${workflow.runName}": "[$workflow.manifest.name] FAILED: ${workflow.runName}"
            // Render the sendmail template
            def smailFields = [
                    to: params.email,
                    subject: subject,
                    text: txtReport,
                    body: htmlReport,
                    attach: mqcReport,
            ]
            Nextflow.sendMail(smailFields)
            log.info "[$workflow.manifest.name] Sent summary e-mail to $params.email (sendmail)"
        }

        // workflowOnComplete file
        File woc = new File(outDir, "workflowOnComplete.txt")
        def endSummary = [
            'Completed on': workflow.complete,
            'Duration': workflow.duration,
            'Success': workflow.success,
            'Exit status': workflow.exitStatus,
            'Error report': workflow.errorReport ?: '-'
        ]
        
        String endWfSummary = endSummary.collect { k,v -> "${k.padRight(30, '.')}: $v" }.join("\n")
        String execInfo = "Execution summary\n${endWfSummary}\n"
        woc.withWriter { w -> w << execInfo }

        def endMessage = workflow.success ? workflow.stats.ignoredCount > 0 ? """\
            ${colors.purple}Warning, pipeline completed, but with errored process(es)${colors.reset}
            ${colors.red}Number of ignored errored process(es) : ${workflow.stats.ignoredCountFmt}${colors.reset}
            ${colors.green}Number of successfully ran process(es) : ${workflow.stats.succeedCountFmt}${colors.reset}
            [$workflow.manifest.name]${colors.green} Pipeline completed successfully${colors.reset}
            """.stripIndent() : "[$workflow.manifest.name]${colors.green} Pipeline completed successfully${colors.reset}": "[$workflow.manifest.name]${colors.red} Pipeline completed with errors${colors.reset}"
        log.info endMessage.toString()
    }

    /**
     * Lint params scope with usage from the nf-core json
     * @param params
     * @param paramsWithUsage
     * @return
     */
    def lint(params, paramsWithUsage) {
        def linter = new ParamsLinter(params, paramsWithUsage, log)
        return linter.lint()
    }

    /**
     * Welcome method which should be launch at the beginning of the workflow
     * @return
     */
    def welcome() {

        def sessionParams = binding.getParams()
        def sessionWorkflow = binding.getVariable('workflow')
        binding.setVariable('colors', generateLogColors(sessionParams.get("monochromeLogs", false) as Boolean))

        nfHeader(sessionParams, sessionWorkflow as WorkflowMetadata)
        if ("${sessionWorkflow.manifest.version}" =~ /dev/ ) {
            def devMessageFile = new File("${sessionWorkflow.projectDir}/assets/devMessage.txt")
            log.info devMessageFile.text
        }
        if (sessionParams.help) {
            def paramsWithUsage = readParamsFromJsonSettings("${sessionWorkflow.projectDir}/parameters.settings.json")
            helpMessage(paramsWithUsage, sessionWorkflow)
            Nextflow.exit(1)
        }

    }

}