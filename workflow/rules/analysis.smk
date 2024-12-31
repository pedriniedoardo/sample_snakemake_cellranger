rule run_cellranger:
    '''
    This is the rule to run cellranger.
    '''
    input:
        # fastqs = "/beegfs/scratch/ric.absinta/ric.absinta/runs/240214_A00626_0755_AHT3N5DMXY_1/Project_AbsintaM_2105_grant_NMSS_RFA_2203_39325/{sample}"
        # fastqs = "/beegfs/scratch/ric.absinta/ric.absinta/runs/240227_A00626_0759_BHMTFLDMXY/AbsintaM_2209_scRNAseq/{sample}"
        # fastqs = "/beegfs/scratch/ric.absinta/ric.absinta/runs/240314_A01966_0015_BHML33DMXY/AbsintaM_2209_scRNAseq/{sample}"
        # fastqs = "/beegfs/scratch/ric.absinta/ric.absinta/runs/240314_A01966_0015_BHML33DMXY/AbsintaM_2105_grant_NMSS_RFA_2203_39325/{sample}"
        # fastqs = lambda w: config["SAMPLES"]["{}".format(w.sample)]['location']
        fastqs = config["FASTQ_location"] + "{sample}"
    conda:
        config['env_cellranger']
    output:
        output = config["out_location"] + "output/{sample}_finished.log",
        # summary = "/media/edo/sandiskSSD/work/training/snakemake/test_cellranger/pipeline/{sample}/outs/web_summary.html",
        folder = directory(config["out_location"] + "{sample}")
    log:
        'logs/{sample}/run_cellranger.log'
    benchmark:
        'benchmarks/{sample}/run_cellranger.txt'
    resources:
        mem_mb = lambda wildcards, attempt: config['cellranger_RAM2'] * attempt,
        # mem_mb = 16000,
        cpus = config['cellranger_CPU']
    threads:
        config['cellranger_CPU']
    params:
        transcriptome = config['transcriptome'],
        expected_cells = lambda w: config["SAMPLES"]["{}".format(w.sample)]['expectCells'],
        # sample-tailored intron flag. comment global intron flag and uncomment the line below to use the sample-tailored intron flag
        INTRON_flag = lambda w: config["SAMPLES"]["{}".format(w.sample)]['INTRON_flag'],
        # global intron flag. comment the sample-tailored intron flag above and uncomment the line below to use the global intron flag
        # INTRON_flag = config['INTRON_flag'],
        BAM_flag = config['BAM_flag'],
        cpus = config['cellranger_CPU'],
        RAM = config['cellranger_RAM'],
        INTRON_param = config['INTRON_param'],
        out_cellranger_wd = '{sample}'
    shell:
        '''
        echo "start cellranger run for <{wildcards.sample}>" >> {log}
         
        cellranger count --id={wildcards.sample} \
        --transcriptome={params.transcriptome} \
        --fastqs={input.fastqs} \
        --sample={wildcards.sample} \
        --localcores={params.cpus} \
        --create-bam={params.BAM_flag} \
        {params.INTRON_param}={params.INTRON_flag} \
        --localmem={params.RAM}
        
        # generate the output file
        touch {output.output}

        echo "completed run cellranger for <{wildcards.sample}>" >> {log}
        '''

rule move_cellranger:
    '''
    This is the rule to move all the outputs to the output folder.
    '''
    input:
        wd_cellranger = rules.run_cellranger.output.folder
    conda:
        config['env_cellranger']
    output:
        folder = directory(config["out_location"] + "output/{sample}"),
        summary = config["out_location"] + "output/web_summaries/{sample}_web_summary.html"
    log:
        'logs/{sample}/02_MoveCellranger.log'
    benchmark:
        'benchmarks/{sample}/02_MoveCellranger.txt'
    resources:
        mem_mb = 250,
        cpus = 1
    threads: 1
    params:
        wd_summary = config["out_location"] + '{sample}/outs/web_summary.html'
    shell:
        '''
        # copy
        echo "copy summary <{wildcards.sample}>" >> {log}

        # copy the web summary from the output in an individula folder
        cp {params.wd_summary} {output.summary} >> {log}

        echo "summary copied <{wildcards.sample}>" >> {log}

        # move
        echo "start move files <{wildcards.sample}>" >> {log}
        
        # move the files to the destination folder
        mv {input.wd_cellranger} {output.folder} >> {log}

        echo "move completed <{wildcards.sample}>" >> {log}
        '''
