#### test
import pprint
import os
import yaml
import sys
from directories import *
import snakemake
snakemake.utils.min_version("6.0")

yaml.warnings({'YAMLLoadWarning': False})
shell.executable("/bin/bash")


SAMPLES = {}
with open(config["sample_sheet"]) as sample_sheet_file:
    SAMPLES = yaml.load(sample_sheet_file)

primerfile = config["primer_file"]
reffile = srcdir("files/MN908947.fasta")
ref_basename = os.path.splitext(os.path.basename(reffile))[0]

rule all:
    input:
        f"{res}multiqc.html",
        f"{res + seqs}" + "concat_cov_ge_50.fasta",
        f"{res}Width_of_coverage.tsv"
    
rule Prepare_ref_and_primers:
    input:
        ref = reffile,
        prm = primerfile
    output: 
        ref = f"{datadir + refdir + ref_basename}.fasta",
        prm = f"{datadir + prim}" + "primers.fasta",
        refindex = f"{datadir + refdir + ref_basename}.fasta.fai"
    conda:
        f"{conda_envs}Alignment.yaml"
    threads: config['threads']['Index']
    shell:
        """
        cat {input.ref} | seqkit replace -p "\-" -s -r "N" > {output.ref}
        cat {input.prm} | seqkit replace -p "\-" -s -r "N" > {output.prm}
        samtools faidx {output.ref} -o {output.refindex}
        """

####
#
#SECTION::
#Load in the data from the input directory and remove adapters.
#Rules are chosen based on input-type [illumina/nanopore]
#
####
if config["platform"] == "illumina":
    rule QC_raw:
        input: lambda wildcards: SAMPLES[wildcards.sample][wildcards.read]
        output:
            html    =   f"{datadir + qc_pre}" + "{sample}_{read}_fastqc.html",
            zip     =   f"{datadir + qc_pre}" + "{sample}_{read}_fastqc.zip"
        conda:
            f"{conda_envs}Clean.yaml"
        log:
            f"{logdir}" + "QC_raw_data_{sample}_{read}.log"
        benchmark:
            f"{logdir + bench}" + "QC_raw_data_{sample}_{read}.txt"
        threads: config['threads']['QC']
        params:
            output_dir  =   f"{datadir + qc_pre}",
            script = srcdir("scripts/fastqc_wrapper.sh")
        shell:
            """
            bash {params.script} {input} {params.output_dir} {output.html} {output.zip} {log}
            """

    rule RemoveAdapters_p1:
        input: 
            ref = rules.Prepare_ref_and_primers.output.ref,
            fq  = lambda wildcards: (SAMPLES[wildcards.sample][i]
                                for i in ("R1", "R2")
                                )
        output:
            bam     = f"{datadir + cln + raln}" + "{sample}.bam",
            index   = f"{datadir + cln + raln}" + "{sample}.bam.bai"
        conda:
            f"{conda_envs}Alignment.yaml"
        log:
            f"{logdir}" + "RemoveAdapters_p1_{sample}.log"
        benchmark:
            f"{logdir + bench}"+ "RemoveAdapters_p1_{sample}.txt"
        threads: config['threads']['Alignments']
        params:
            mapthreads = config['threads']['Alignments'] - 1,
            filters = config["runparams"]["alignmentfilters"]
        shell: 
            """
            minimap2 -ax sr -t {params.mapthreads} {input.ref} {input.fq[0]:q} {input.fq[1]:q} 2>> {log} |\
            samtools view -@ {threads} {params.filters} -uS 2>> {log} |\
            samtools sort -o {output.bam} >> {log} 2>&1
            samtools index {output.bam} >> {log} 2>&1
            """
    
    rule RemoveAdapters_p2:
        input: rules.RemoveAdapters_p1.output.bam
        output: f"{datadir + cln + noad}" + "{sample}.fastq"
        conda:
            f"{conda_envs}Clean.yaml"
        threads: config['threads']['AdapterRemoval']
        params:
            script = srcdir('scripts/clipper.py')
        shell: 
            """
            python {params.script} --input {input} --output {output} --threads {threads}
            """

    rule QC_filter:
        input: rules.RemoveAdapters_p2.output
        output: 
            fq = f"{datadir + cln + qcfilt}" + "{sample}.fastq",
            html = f"{datadir + cln + qcfilt + html}" + "{sample}.fastp.html",
            json = f"{datadir + cln + qcfilt + json}" + "{sample}.fastp.json"
        conda: 
            f"{conda_envs}Clean.yaml"
        log:
            f"{logdir}" + "Cleanup_{sample}.log"
        benchmark:
            f"{logdir + bench}" + "Cleanup_{sample}.txt"
        threads: config['threads']['QC']
        params:
            score = config['runparams']['qc_filter_illumina'],
            size = config['runparams']['qc_window_illumina'],
            length = config['runparams']['qc_min_readlength']
        shell:
            """
            fastp --thread {threads} -i {input} \
            --cut_right -M {params.score} -W {params.size} -l {params.length} \
            -o {output.fq} -h {output.html} -j {output.json} > {log} 2>&1
            """

if config["platform"] == "nanopore":
    rule QC_raw:
        input: lambda wildcards: SAMPLES[wildcards.sample]
        output:
            html    =   f"{datadir + qc_pre}" + "{sample}_fastqc.html",
            zip     =   f"{datadir + qc_pre}" + "{sample}_fastqc.zip"
        conda:
            f"{conda_envs}Clean.yaml"
        log:
            f"{logdir}" + "QC_raw_data_{sample}.log"
        benchmark:
            f"{logdir + bench}" + "QC_raw_data_{sample}.txt"
        threads: config['threads']['QC']
        params:
            output_dir  =   f"{datadir + qc_pre}",
            script = srcdir("scripts/fastqc_wrapper.sh")
        shell:
            """
            bash {params.script} {input} {params.output_dir} {output.html} {output.zip} {log}
            """

    rule RemoveAdapters_p1:
        input: 
            ref = rules.Prepare_ref_and_primers.output.ref,
            fq  = lambda wildcards: SAMPLES[wildcards.sample]
        output:
            bam     = f"{datadir + cln + raln}" + "{sample}.bam",
            index   = f"{datadir + cln + raln}" + "{sample}.bam.bai"
        conda:
            f"{conda_envs}Alignment.yaml"
        log:
            f"{logdir}" + "RemoveAdapters_p1_{sample}.log"
        benchmark:
            f"{logdir + bench}"+ "RemoveAdapters_p1_{sample}.txt"
        threads: config['threads']['Alignments']
        params:
            mapthreads = config['threads']['Alignments'] - 1,
            filters = config["runparams"]["alignmentfilters"]
        shell: 
            """
            minimap2 -ax sr -t {params.mapthreads} {input.ref} {input.fq} 2>> {log} |\
            samtools view -@ {threads} {params.filters} -uS 2>> {log} |\
            samtools sort -o {output.bam} >> {log} 2>&1
            samtools index {output.bam} >> {log} 2>&1
            """
    
    rule RemoveAdapters_p2:
        input: rules.RemoveAdapters_p1.output.bam
        output: f"{datadir + cln + noad}" + "{sample}.fastq"
        conda:
            f"{conda_envs}Clean.yaml"
        threads: config['threads']['AdapterRemoval']
        params:
            script = srcdir('scripts/clipper.py')
        shell: 
            """
            python {params.script} --input {input} --output {output} --threads {threads}
            """

    rule QC_filter:
        input: rules.RemoveAdapters_p2.output
        output: 
            fq = f"{datadir + cln + qcfilt}" + "{sample}.fastq",
            html = f"{datadir + cln + qcfilt + html}" + "{sample}.fastp.html",
            json = f"{datadir + cln + qcfilt + json}" + "{sample}.fastp.json"
        conda: 
            f"{conda_envs}Clean.yaml"
        log:
            f"{logdir}" + "Cleanup_{sample}.log"
        benchmark:
            f"{logdir + bench}" + "Cleanup_{sample}.txt"
        threads: config['threads']['QC']
        params:
            score = config['runparams']['qc_filter_nanopore'],
            size = config['runparams']['qc_window_nanopore'],
            length = config['runparams']['qc_min_readlength']
        shell:
            """
            fastp --thread {threads} -i {input} \
            --cut_right -M {params.score} -W {params.size} -l {params.length} \
            -o {output.fq} -h {output.html} -j {output.json} > {log} 2>&1
            """

rule RemovePrimers:
    input: 
        fq = rules.QC_filter.output.fq,
        pr = rules.Prepare_ref_and_primers.output.prm,
        ref = rules.Prepare_ref_and_primers.output.ref
    output: f"{datadir + cln + prdir}" + "{sample}.fastq"
    conda: 
        f"{conda_envs}Clean.yaml"
    log:
        f"{logdir}" + "RemovePrimers_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "RemovePrimers_{sample}.txt"
    threads: config['threads']['PrimerRemoval']
    params:
        amplicontype = config["amplicon_type"]
    shell:
        """
        AmpliGone -i {input.fq} \
        -ref {input.ref} -pr {input.pr} \
        -o {output} -at {params.amplicontype} \
        -t {threads}
        """

rule QC_clean:
    input: rules.RemovePrimers.output
    output: 
        html    =   f"{datadir + qc_post}" + "{sample}_fastqc.html",
        zip     =   f"{datadir + qc_post}" + "{sample}_fastqc.zip"
    conda: 
        f"{conda_envs}Clean.yaml"
    log:
        f"{logdir}" + "QC_clean_data_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "QC_clean_data_{sample}.txt"
    threads: config['threads']['QC']
    params:
        outdir = f"{datadir + qc_post}"
    shell:
        """
        if [ -s "{input}" ]; then
            fastqc -t {threads} --quiet --outdir {params.outdir} {input} > {log} 2>&1
        else
            touch {output.html}
            touch {output.zip}
        fi
        """


### Align the cleaned reads to the reference
if config["platform"] == "illumina":
    rule Alignment:
        input:
            fq = rules.RemovePrimers.output,
            ref = rules.Prepare_ref_and_primers.output.ref
        output:
            bam = f"{datadir + aln + bf}" + "{sample}.bam",
            index = f"{datadir + aln + bf}" + "{sample}.bam.bai"
        conda: 
            f"{conda_envs}Alignment.yaml"
        log:
            f"{logdir}" + "Alignment_{sample}.log"
        benchmark:
            f"{logdir + bench}" + "Alignment_{sample}.txt"
        threads: config['threads']['Alignments']
        params:
            mapthreads = config['threads']['Alignments'] - 1,
            filters = config["runparams"]["alignmentfilters"]
        shell: 
            """
            minimap2 -ax sr -t {params.mapthreads} {input.ref} {input.fq} 2>> {log} |\
            samtools view -@ {threads} {params.filters} -uS 2>> {log} |\
            samtools sort -o {output.bam} >> {log} 2>&1
            samtools index {output.bam} >> {log} 2>&1
            """

if config["platform"] == "nanopore":
    rule Alignment:
        input:
            fq = rules.RemovePrimers.output,
            ref = rules.Prepare_ref_and_primers.output.ref
        output:
            bam = f"{datadir + aln + bf}" + "{sample}.bam",
            index = f"{datadir + aln + bf}" + "{sample}.bam.bai"
        conda: 
            f"{conda_envs}Alignment.yaml"
        log:
            f"{logdir}" + "Alignment_{sample}.log"
        benchmark:
            f"{logdir + bench}" + "Alignment_{sample}.txt"
        threads: config['threads']['Alignments']
        params:
            mapthreads = config['threads']['Alignments'] - 1,
            filters = config["runparams"]["alignmentfilters"]
        shell: 
            """
            minimap2 -ax map-ont -t {params.mapthreads} {input.ref} {input.fq} 2>> {log} |\
            samtools view -@ {threads} {params.filters} -uS 2>> {log} |\
            samtools sort -o {output.bam} >> {log} 2>&1
            samtools index {output.bam} >> {log} 2>&1
            """

rule Consensus:
    input: 
        bam = rules.Alignment.output.bam,
        gff = srcdir("files/MN908947.gff"),
        ref = rules.Prepare_ref_and_primers.output.ref
    output: 
        cons_1 = f"{datadir + cons + seqs}" + "{sample}_cov_ge_1.fa",
        cons_5 = f"{datadir + cons + seqs}" + "{sample}_cov_ge_5.fa",
        cons_10 = f"{datadir + cons + seqs}" + "{sample}_cov_ge_10.fa",
        cons_50 = f"{datadir + cons + seqs}" + "{sample}_cov_ge_50.fa",
        cons_100 = f"{datadir + cons + seqs}" + "{sample}_cov_ge_100.fa",
        cov = f"{datadir + cons + covs}" + "{sample}_coverage.tsv",
        gff = f"{datadir + cons + features}" + "{sample}.gff",
        vcf_1 = f"{datadir + aln + vf}" + "{sample}_cov_ge_1.vcf",
        vcf_5 = f"{datadir + aln + vf}" + "{sample}_cov_ge_5.vcf",
        vcf_10 = f"{datadir + aln + vf}" + "{sample}_cov_ge_10.vcf",
        vcf_50 = f"{datadir + aln + vf}" + "{sample}_cov_ge_50.vcf",
        vcf_100 = f"{datadir + aln + vf}" + "{sample}_cov_ge_100.vcf"
    params:
        mincov = "1 5 10 50 100",
        outdir = f"{datadir + cons + seqs}",
        vcfdir = f"{datadir + aln + vf}"
    conda: 
        f"{conda_envs}Consensus.yaml"
    log:
        f"{logdir}" + "Consensus_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Consensus_{sample}.txt"
    threads: config['threads']['Consensus']
    shell:
        """
        TrueConsense --input {input.bam} \
        --reference {input.ref} --features {input.gff} \
        --coverage-levels {params.mincov} \
        --samplename {wildcards.sample} \
        --output {params.outdir} \
        --variants {params.vcfdir} \
        --output-gff {output.gff} \
        --depth-of-coverage {output.cov} \
        --noambiguity --threads {threads}
        """

rule concat_seqs:
    input:
        cons = expand( "{p}{sample}_cov_ge_{t}.fa",
                p = f"{datadir + cons + seqs}",
                t = [1, 5, 10, 50, 100],
                sample = SAMPLES
                )
    output:
        cons_1 = f"{res + seqs}" + "concat_cov_ge_1.fasta",
        cons_5 = f"{res + seqs}" + "concat_cov_ge_5.fasta",
        cons_10 = f"{res + seqs}" + "concat_cov_ge_10.fasta",
        cons_50 = f"{res + seqs}" + "concat_cov_ge_50.fasta",
        cons_100 = f"{res + seqs}" + "concat_cov_ge_100.fasta",
    params:
        wcar_1 = f"{datadir + cons + seqs}" + "*_cov_ge_1.fa",
        wcar_5 = f"{datadir + cons + seqs}" + "*_cov_ge_5.fa",
        wcar_10 = f"{datadir + cons + seqs}" + "*_cov_ge_10.fa",
        wcar_50 = f"{datadir + cons + seqs}" + "*_cov_ge_50.fa",
        wcar_100 = f"{datadir + cons + seqs}" + "*_cov_ge_100.fa",
    shell:
        """
        cat {params.wcar_1} >> {output.cons_1}
        cat {params.wcar_5} >> {output.cons_5}
        cat {params.wcar_10} >> {output.cons_10}
        cat {params.wcar_50} >> {output.cons_50}
        cat {params.wcar_100} >> {output.cons_100}
        """


###
#> add some rules here for concatenating the vcf files to a single tsv
###

rule Get_Breadth_of_coverage:
    input:
        reference = rules.Prepare_ref_and_primers.output.ref,
        coverage = rules.Consensus.output.cov,
    output:
        temp(f"{datadir + boc}" + "{sample}.tsv")
    conda:
        f"{conda_envs}Consensus.yaml"
    threads: 1
    params:
        script = srcdir("scripts/boc.py")
    shell:
        """
        python {params.script} {input.reference} {wildcards.sample} {input.coverage} {output}
        """

rule concat_boc:
    input:
        pct = expand("{p}{sample}.tsv",
                    p = f"{datadir + boc}",
                    sample = SAMPLES)
    output: f"{res}Width_of_coverage.tsv"
    threads: 1
    shell:
        """
        echo -e "Sample_name\tWidth_at_mincov_1\tWidth_at_mincov_5\tWidth_at_mincov_10\tWidth_at_mincov_50\tBoC_at_coverage_threshold_100" > {output}
        cat {input} >> {output}
        """

if config['platform'] == "illumina":
    rule MultiQC_report:
        input: 
            expand("{p}{sample}_{read}_fastqc.zip",
                p = f"{datadir + qc_pre}",
                sample = SAMPLES,
                read = "R1 R2".split()
                ),
            expand("{p}{sample}_fastqc.zip",
                p = f"{datadir + qc_post}",
                sample = SAMPLES
                )
        output: 
            f"{res}multiqc.html",
            expand("{p}multiqc_{program}.txt",
                p = f"{res + mqc_data}",
                program = ['fastqc'])
        conda:
            f"{conda_envs}Clean.yaml"
        log:
            f"{logdir}" + "MultiQC_report.log"
        benchmark:
            f"{logdir + bench}" + "MultiQC_report.txt"
        threads: 1
        params:
            conffile = srcdir('files/multiqc_config.yaml'),
            outdir = f"{res}"
        shell:
            """
            multiqc -d --force --config {params.conffile} -o {params.outdir} -n multiqc.html {input} > {log} 2>&1
            """

if config['platform'] == "nanopore":
    rule MultiQC_report:
        input: 
            expand("{p}{sample}_fastqc.zip",
                p = f"{datadir + qc_pre}",
                sample = SAMPLES
                ),
            expand("{p}{sample}_fastqc.zip",
                p = f"{datadir + qc_post}",
                sample = SAMPLES
                )
        output: 
            f"{res}multiqc.html",
            expand("{p}multiqc_{program}.txt",
                p = f"{res + mqc_data}",
                program = ['fastqc'])
        conda:
            f"{conda_envs}Clean.yaml"
        log:
            f"{logdir}" + "MultiQC_report.log"
        benchmark:
            f"{logdir + bench}" + "MultiQC_report.txt"
        threads: 1
        params:
            conffile = srcdir('files/multiqc_config.yaml'),
            outdir = f"{res}"
        shell:
            """
            multiqc -d --force --config {params.conffile} -o {params.outdir} -n multiqc.html {input} > {log} 2>&1
            """

