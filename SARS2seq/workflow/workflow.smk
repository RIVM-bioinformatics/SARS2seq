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
        f"{res + seqs}concat_cov_ge_50.fasta",
        f"{res + seqs}concat_mutations_cov_ge_50.tsv",
        f"{res}Width_of_coverage.tsv",
        f"{res}Typing_results.tsv",
        f"{res}" + "annotation_happened.txt",
    
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

#TODO: Add one of the following lines to vcf headers (see https://samtools.github.io/hts-specs/VCFv4.2.pdf chapter 1.2.7): 
    ##contig=<ID=MN908947.3>
    ##contig=<ID=MN908947.3,length=29903>
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
        vcf_100 = f"{datadir + aln + vf}" + "{sample}_cov_ge_100.vcf",
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

rule fetch_problematic_sites:
    output:
        vcf         = f"{datadir + ann}" + "problematic_sites.vcf",
        vcf_gz      = f"{datadir + ann}" + "problematic_sites.vcf.gz",
        vcf_gz_tbi  = f"{datadir + ann}" + "problematic_sites.vcf.gz.tbi",
    conda:
        f"{conda_envs}Mutations.yaml"
    params:
        url = 'https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf'
    log:
        f"{logdir}fetch_problematic_sites.log"
    shell: 
        """
        curl -vs {params.url} 2> {log} > {output.vcf} 2>> {log}
        bgzip -c {output.vcf} > {output.vcf_gz} 2>> {log}
        tabix {output.vcf_gz} 2>> {log}
        """

rule annotate_problematic_sites:
    input:
        vcf = f"{datadir + aln + vf}" + "{sample}_cov_ge_{cov}.vcf",
        problematic_sites = rules.fetch_problematic_sites.output.vcf_gz,
    output:
        vcf = f"{datadir + aln + vf + ann}" + "{sample}_cov_ge_{cov}_annotated.vcf",
        vcf_gz = temp(f"{datadir + aln + vf + ann}" + "{sample}_cov_ge_{cov}_annotated.vcf.gz"),
        vcf_gz_tbi = temp(f"{datadir + aln + vf + ann}" + "{sample}_cov_ge_{cov}_annotated.vcf.gz.tbi"),
    conda: f"{conda_envs}Mutations.yaml"
    shadow: "shallow"
    log:
        f"{logdir}" + "annotate_{sample}_cov_ge_{cov}.log"
    shell:
        """
        bgzip -c {input.vcf} > {output.vcf_gz} 2> {log}
        tabix -f {output.vcf_gz} 2>> {log}
        bcftools annotate -Ov -a {input.problematic_sites} {output.vcf_gz} \
            -c FILTER,EXC,GENE,AA_POS,AA_REF,AA_ALT > {output.vcf} 2>> {log}
        """

rule annotation_happened:
    input:
        expand(f"{datadir + aln + vf + ann}" + "{sample}_cov_ge_{cov}_annotated.vcf",
            sample = SAMPLES, cov = [1,5,10,50,100]
        ),
    output:
        temp(f"{res}" + "annotation_happened.txt")
    shell: "touch {output}"

rule vcf_to_tsv:
    input:
        f"{datadir + aln + vf}" + "{sample}_cov_ge_{cov}.vcf",
    output:
        temp(f"{datadir + aln + vf}" + "{sample}_cov_ge_{cov}.tsv"),
    conda: f"{conda_envs}Mutations.yaml"
    log:
        f"{logdir}" + "vcf_to_tsv_{sample}_cov_ge_{cov}.log"
    shell:
        """
        bcftools query {input} -f '{wildcards.sample}\t%CHROM\t%POS\t%REF\t%ALT\t%DP\n' -e 'ALT="N"' > {output} 2>> {log}
        """

rule concat_tsvs:
    input:
        expand(
            "{p}{sample}_cov_ge_50.tsv",
            p = datadir + aln + vf,
            sample = SAMPLES
            )
    output:
        f"{res + seqs}concat_mutations_cov_ge_50.tsv",
    log:
        f"{logdir}" + "concat_tsvs.log"
    run:
        shell("echo -e 'Sample\tReference_Chromosome\tPosition\tReference\tAlternative\tDepth' > {output} 2> {log}")
        shell("cat {input} >> {output} 2>> {log}")


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

rule Catch_typing_versions:
    output:
        pangolin = temp(f"{datadir}" + "pangolin.version"),
        nextclade= temp(f"{datadir}" + "nextclade.version")
    conda:
        f"{conda_envs}Typing.yaml"
    threads: 1
    shell:
        """
        pangolin -v > {output.pangolin}
        nextclade --version > {output.nextclade}
        """

rule Typing:
    input: rules.Consensus.output.cons_1
    output: 
        pango = temp(f"{datadir + cons + tbl}" + "{sample}_pangolin.csv"),
        nextc = temp(f"{datadir + cons + tbl}" + "{sample}_nextclade.csv")
    conda:
        f"{conda_envs}Typing.yaml"
    threads: config['threads']['Typing']
    log:
        f"{logdir}" + "Typing_{sample}.log"
    params:
        pango_dir = f"{datadir + cons + tbl}",
    shell:
        """
        nextclade -i {input} -c {output.nextc} -j {threads}
        pangolin {input} -o {params.pango_dir} --outfile {wildcards.sample}_pangolin.csv > {log} 2>&1
        """

rule format_typing:
    input:
        pango = rules.Typing.output.pango,
        nextc = rules.Typing.output.nextc,
        pangv = rules.Catch_typing_versions.output.pangolin,
        nextv = rules.Catch_typing_versions.output.nextclade
    output: f"{datadir + cons + tbl}" + "{sample}_typingresults.tsv"
    conda:
        f"{conda_envs}Typing.yaml"
    params:
        script = srcdir('scripts/typingagg.py')
    shell:
        """
        python {params.script} {wildcards.sample} {input.nextv} {input.pangv} {input.nextc} {input.pango} {output}
        """

rule combine_typing_results:
    input: 
        expand( "{p}{sample}_typingresults.tsv",
                    p = f"{datadir + cons + tbl}",
                    sample = SAMPLES)
    output: f"{res}Typing_results.tsv"
    shell:
        """
        echo -e "Sample_name\tTyping_date\tPangolin version\tNextClade version\tPangolin lineages version\tPangolin Lineage\tNextClade Clade\tProbability\tPangolin status\tNextClade QC" > {output}
        cat {input} >> {output}
        """