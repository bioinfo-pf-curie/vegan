#!/usr/bin/env nextflow


Channel.fromPath(params.ChrList)
       .splitCsv(header:false, sep:'\t')
       .map{ row-> params.resultdir + "/" + row[0] }
       .set{ChrFileNamePrfx}

process mpibwa_process {

    cpus 1
    executor 'slurm'
    memory '1GB'
    publishDir "${params.resultdir}"

    output:
    val (bwa_mem) into channel1 

    script:

    bwa_mem=1    
    """
    sbatch --wait -p dev -N 3 -n 12 -c 8 --ntasks-per-node=4 --mem-per-cpu=2G ${params.MPIBWAWRAP} ${params.fastq1} ${params.fastq2} ${params.RefGenome} ${params.Header} ${params.MPIBWA}
    """
}

process name_mpisort_process {

    cpus 1
    executor 'slurm'
    memory '1GB'
    publishDir "${params.resultdir}"
    
    input: 
        file(file_prfx) from ChrFileNamePrfx
        val(alignment) from channel1          

    output:
	file(file_mrkdp) into channel2   

    shell:    

    file_mrkdp=file_prfx + "_mrkdp.bam"
    ok=1  
    """
    file_sam=!{file_prfx}".sam"
    sam_size=\$(du -h \${file_sam} | awk '{print \$1}')
    sam_size=\$(echo \${sam_size} | sed 's/,/./')
    unit=\${sam_size:0-1}
    sam_size=\${sam_size::-1} 
    if [[ "\${unit}" == "M" ]]; then
    	sam_size=3
    fi
         
    
    slrm_opt=\$(python ${params.PythonScript} -c 20 -m 190 -s \${sam_size} -t uniq) &&   
    sbatch --wait -p dev \${slrm_opt} !{params.MPISortQuery} \${file_sam} !{params.MPISORT} !{params.resultdir}
    file_bam=!{file_prfx}".bam" &&
    file_query_srt=!{file_prfx}"_querysort.bam" &&
    mv \${file_bam} \${file_query_srt} &&
    !{params.SAMTOOLS} fixmate -@ 8 -m -O sam \${file_query_srt} !{file_prfx}"_fixmate.sam" &&

    file_sam_fxmt=!{file_prfx}"_fixmate.sam" &&
    echo \${file_sam_fxmt} " " \${slrm_opt}
    sbatch --wait -p dev \${slrm_opt} !{params.MPISortCoord} \${file_sam_fxmt} !{params.MPISORT} !{params.resultdir}

    file_bam=!{file_prfx}".bam"
    file_mrkdp=!{file_prfx}"_mrkdp.bam"
    !{params.SAMTOOLS} index  \${file_bam} &&
    !{params.SAMTOOLS} markdup -@ 8 -O bam \${file_bam}  \${file_mrkdp} &&
    !{params.SAMTOOLS} index \${file_mrkdp}

    """
}


process merge_bam {

        cpus 1
        executor 'slurm'
        memory '1GB'
        publishDir "${params.resultdir}"

        input:
        val file_bam from channel2.collect()

        script:
        
        """
        ${params.SAMTOOLS} view -@ 8 -b -o ${params.SAMTOOLS}"/unmapped.bam" ${params.SAMTOOLS}"/unmapped.sam" &&
        ${params.SAMTOOLS} merge -@ 8 -h ${params.resultdir}"/HG001.sam" ${params.resultdir}"/HG001.bam" ${file_bam.join(" ")} ${params.SAMTOOLS}"/unmapped.bam"   &&
        ${params.SAMTOOLS} index ${params.resultdir}"/HG001.bam"
        """
}

