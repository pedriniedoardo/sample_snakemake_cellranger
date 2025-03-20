rule runCellranger:
    '''
    This is the rule to run cellranger.
    '''
    input:
        fastqs = config["FASTQ_location"] + "{sample}"
    conda:
        config['env_cellranger']
    output:
        output = config["out_location"] + "cellranger/{sample}_finished.log",
        folder = directory(config["out_location"] + "cellranger/{sample}")
    log:
        'logs/{sample}/runCellranger.log'
    benchmark:
        'benchmarks/{sample}/runCellranger.txt'
    resources:
        mem_gb = lambda wildcards, attempt: config['cellranger_RAM'] * attempt,
        # mem_mb = 16000,
        cpus = config['cellranger_CPU']
    threads:
        config['cellranger_CPU']
    params:
        transcriptome = config['transcriptome'],
        expected_cells = lambda w: config["SAMPLES"]["{}".format(w.sample)]['expectCells'],
        # sample-tailored intron flag. comment global intron flag and uncomment the line below to use the sample-tailored intron flag
        # INTRON_flag = lambda w: config["SAMPLES"]["{}".format(w.sample)]['INTRON_flag'],
        # global intron flag. comment the sample-tailored intron flag above and uncomment the line below to use the global intron flag
        INTRON_flag = config['INTRON_flag'],
        BAM_flag = config['BAM_flag'],
        cpus = config['cellranger_CPU'],
        RAM = config['cellranger_RAM'],
        INTRON_param = config['INTRON_param'],
        outs = directory(config["out_location"] + "cellranger/{sample}")
    shell:
        '''
        echo "generate the output folder" >> {log}
        mkdir -p {params.outs} >> {log}
        
        echo "start cellranger run for <{wildcards.sample}>" >> {log}
         
        cellranger count --id={wildcards.sample} \
        --transcriptome={params.transcriptome} \
        --fastqs={input.fastqs} \
        --sample={wildcards.sample} \
        --localcores={params.cpus} \
        --create-bam={params.BAM_flag} \
        {params.INTRON_param}={params.INTRON_flag} \
        --localmem={params.RAM} \
        --output-dir={params.outs}
        
        # generate the output file
        touch {output.output}

        echo "completed run cellranger for <{wildcards.sample}>" >> {log}
        '''

rule moveCellrangerSummary:
    '''
    This is the rule to move the outputs of interest.
    '''
    input:
        cellranger_summary = rules.runCellranger.output.output
    conda:
        config['env_cellranger']
    output:
        summary = config["out_location"] + "web_summaries/{sample}_web_summary.html",
    log:
        'logs/{sample}/02_MoveCellranger_summary.log'
    benchmark:
        'benchmarks/{sample}/02_MoveCellranger_summary.txt'
    resources:
        mem_mb = 250,
        cpus = 1
    threads: 1
    params:
        wd_summary = config["out_location"] + 'cellranger/{sample}/outs/web_summary.html',
        bam = config["out_location"] + "cellranger/{sample}/outs/possorted_genome_bam.bam"
    shell:
        '''
        # copy
        echo "copy summary <{wildcards.sample}>" >> {log}

        # copy the web summary from the output in an individula folder
        cp {params.wd_summary} {output.summary} >> {log}

        echo "summary copied <{wildcards.sample}>" >> {log}
        '''

rule runMultiqc1:
    '''
    This rule allow to run the multiqc on the final folder.
    Currently it is set to run after cellranger. this is the only step to grep some outputs.
    '''
    input:
        # this is needed to trigger it after the generation of the outputs
        input_file = expand(rules.moveCellrangerSummary.output.summary, sample=config['SAMPLES'])
    conda:
        config['env_bioinfo']
    output:
        html = config["out_location"] + "multiQC/multiqc_report.html"
    log:
        'logs/multiqc.log'
    benchmark:
        'benchmarks/multiqc.txt'
    resources:
        mem_gb = 4,
        cpus = 2
    threads: 2
    params:
        folder_in = config["out_location"] + "cellranger/",
        folder_out = config["out_location"] + "multiQC/",
    shell:
        '''
        echo "start multiqc" >> {log}
        multiqc {params.folder_in} -o {params.folder_out}
        echo "end multiqc" >> {log}
        '''

# ---------------------------------------------------------------------------- #
rule test:
    '''
    This is a test rule
    '''
    input:
        test = config['test01']
    conda:
        config['env_cellranger']
    output:
        test = config['test_out']
    log:
        'logs/test/test.log'
    benchmark:
        'benchmarks/test/test.txt'
    resources:
        mem_mb = 500,
        cpus = 1
    threads: 1
    params:
        annotations = config['test02']
    shell:
        '''
        echo "run the test" > {log}
        cat {input.test} > {output.test}
        cat {params.annotations} | awk 'NR<=10' >> {output.test}
        echo "test done" >> {log}
        '''
