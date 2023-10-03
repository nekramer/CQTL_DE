#!/usr/bin/env python3

import pandas as pd
import os, shutil
import glob
from datetime import date

## Load config file
configfile: "config/config.yaml"

## Read in samplesheet
samples = pd.read_csv(config["samplesheet"],sep = ",")

## Convert samplesheet columns to strings
samples = samples.astype(str)

## Concatenate Sequencing_Directory to Read1 and Read2 for full read paths
samples['Read1'] = samples[['Sequencing_Directory', 'Read1']].apply(lambda row: os.path.join(*row), axis=1)
samples['Read2'] = samples[['Sequencing_Directory', 'Read2']].apply(lambda row: os.path.join(*row), axis=1)

## Group Seq_Reps
samples['id'] = samples[['Proj', 'Donor']].agg('_'.join, axis=1) + '_R_' + samples[['Condition', 'Time', 'Tech_Rep']].agg('_'.join, axis=1)

## Extract grouped read1 and read2s
read1 = samples.groupby(['id'])['Read1'].apply(list).to_dict()
read2 = samples.groupby(['id'])['Read2'].apply(list).to_dict()

rule all:
    input:
        'output/qc/multiqc_report.html',
        'output/' + str(date.today()) + '_gse.rda'

rule catR1:
    input:
        lambda wildcards: read1.get(wildcards.group)
    output:
        "output/{group}/fastq/{group}_R1.fastq.gz"
    threads: 1
    log:
        err = "output/{group}/logs/{group}_catR1.err"
    shell:
        """
        mkdir -p output/{wildcards.group}/fastq
        cat {input} > {output} 2> {log.err}
        """

rule catR2:
    input:
        lambda wildcards: read2.get(wildcards.group)
    output:
        "output/{group}/fastq/{group}_R2.fastq.gz"
    threads: 1
    log:
        err = "output/{group}/logs/{group}_catR2.err"
    shell:
        """
        mkdir -p output/{wildcards.group}/fastq
        cat {input} > {output} 2> {log.err}
        """

rule qc:
    input:
        R1 = lambda wildcards: ['output/{group}/fastq/{group}_R1.fastq.gz'.format(group=wildcards.group)],
        R2 = lambda wildcards: ['output/{group}/fastq/{group}_R2.fastq.gz'.format(group=wildcards.group)]
    output:
        zips = expand('output/qc/{{group}}_{R}_fastqc.zip', R=['R1', 'R2']),
        html = expand('output/qc/{{group}}_{R}_fastqc.html',R=['R1', 'R2'])
    threads: 2
    params:
        version = config['fastqcVersion']
    log:
        err = "output/{group}/logs/{group}_qc.err"
    shell:
        """
        module load fastqc/{params.version}
        mkdir -p output/qc
        fastqc -t {threads} -o output/qc {input.R1} {input.R2} 2> {log.err}
        """

rule trim:
    input:
        R1 = lambda wildcards: ['output/{group}/fastq/{group}_R1.fastq.gz'.format(group=wildcards.group)],
        R2 = lambda wildcards: ['output/{group}/fastq/{group}_R2.fastq.gz'.format(group=wildcards.group)]
    output:
        trim1 = temp("output/{group}/trim/{group}_R1_val_1.fq.gz"),
        trim2 = temp("output/{group}/trim/{group}_R2_val_2.fq.gz"),
        report1 = temp("output/{group}/trim/{group}_R1.fastq.gz_trimming_report.txt"),
        report2 = temp("output/{group}/trim/{group}_R2.fastq.gz_trimming_report.txt")
    threads: 4
    params:
        version = config['cutadaptVersion']
    log:
        err = "output/{group}/logs/{group}_trim.err"
    shell:
        """
        module load cutadapt/{params.version}
        module load python/3.6.6
        module load pigz
        mkdir -p output/{wildcards.group}/trim
        trim_galore -o output/{wildcards.group}/trim --cores {threads} --path_to_cutadapt /nas/longleaf/apps/cutadapt/2.9/venv/bin/cutadapt --paired {input.R1} {input.R2} 2> {log.err}
        """

rule quant:
    input:
        trim1 = rules.trim.output.trim1,
        trim2 = rules.trim.output.trim2
    output:
        "output/quant/{group}/quant.sf"
    params:
        version = config['salmonVersion'],
        index = config['salmon'],
        gcFlag = config['gcBias'],
        seqFlag = config['seqBias']
    log:
        out = 'output/logs/quant_{group}.out',
        err = 'output/logs/quant_{group}.err'
    shell:
        """
        module load salmon/{params.version}

        if [ {params.gcFlag} == "TRUE" ] && [ {params.seqFlag} == "TRUE" ]; then
            salmon quant --writeUnmappedNames -l A -1 {input.trim1} -2 {input.trim2} -i {params.index} -o output/quant/{wildcards.group} --threads 1 --seqBias --gcBias 1> {log.out} 2> {log.err}
        elif [ {params.gcFlag} == "TRUE" ] && [ {params.seqFlag} != "TRUE" ]; then
            salmon quant --writeUnmappedNames -l A -1 {input.trim1} -2 {input.trim2} -i {params.index} -o output/quant/{wildcards.group} --threads 1 --gcBias 1> {log.out} 2> {log.err}
        else
            salmon quant --writeUnmappedNames -l A -1 {input.trim1} -2 {input.trim2} -i {params.index} -o output/quant/{wildcards.group} --threads 1 --seqBias 1> {log.out} 2> {log.err}
        fi
        """

rule multiqc:
    input:
        [expand('output/qc/{group}_{R}_fastqc.zip', group = key, R = ['R1', 'R2']) for key in read1],
        [expand('output/{group}/trim/{group}_{R}.fastq.gz_trimming_report.txt', group = key, R =['R1', 'R2']) for key in read1],
        [expand('output/quant/{group}/quant.sf', group = key) for key in read1]
    output:
        'output/qc/multiqc_report.html'
    params:
        version = config['multiqcVersion']

    log:
        out = 'output/logs/multiqc.out',
        err = 'output/logs/multiqc.err'
    shell:
        """
        module load multiqc/{params.version}
        multiqc -f -o output/qc . 1> {log.out} 2> {log.err}
        """

rule tximport:
    input:
        [expand('output/quant/{group}/quant.sf', group = key) for key in read1]
    output:
        'output/' + str(date.today()) + '_gse.rda'
    params:
        version = config['rVersion'],
        samplesheet = config['samplesheet']
    log:
        out = 'output/logs/tximport.out',
        err = 'output/logs/tximport.err'
    shell:
        """
        module load r/{params.version}
        Rscript tximport.R {params.samplesheet} 1> {log.out} 2> {log.err}
        """