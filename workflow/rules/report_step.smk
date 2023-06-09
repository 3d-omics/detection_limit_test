rule report_step_reads:
    """Collect all reports for the reads step"""
    input:
        rules.reads_fastqc_all.input,
    output:
        html=REPORT_STEP / "reads.html",
    log:
        REPORT_STEP / "reads.log",
    conda:
        "../envs/report.yml"
    params:
        dir=REPORT_STEP,
    shell:
        """
        multiqc \
            --filename reads \
            --title reads \
            --force \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report_step_fastp:
    """Collect all reports for the fastp step"""
    input:
        rules.fastp_report_all.input,
    output:
        html=REPORT_STEP / "fastp.html",
    log:
        REPORT_STEP / "fastp.log",
    conda:
        "../envs/report.yml"
    params:
        dir=REPORT_STEP,
    shell:
        """
        multiqc \
            --title fastp \
            --force \
            --filename fastp \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report_step_bowtie2:
    """Collect all reports for the bowtie2 step"""
    input:
        rules.bowtie2_report_all.input,
    output:
        html=REPORT_STEP / "bowtie2.html",
    log:
        REPORT_STEP / "bowtie2.log",
    conda:
        "../envs/report.yml"
    params:
        dir=REPORT_STEP,
    shell:
        """
        multiqc \
            --title bowtie2 \
            --force \
            --filename bowtie2 \
            --outdir {params.dir} \
            --dirs \
            --dirs-depth 1 \
            {input} \
        2> {log} 1>&2
        """


rule report_step_bowtie2_human:
    """Collect all reports for the bowtie2 step"""
    input:
        reports=[
            BOWTIE2_HUMAN / f"{sample}.{library}.{report}"
            for sample, library in SAMPLE_LIB
            for report in BAM_REPORTS
        ],
    output:
        html=REPORT_STEP / "bowtie2_human.html",
    log:
        REPORT_STEP / "bowtie2_human.log",
    conda:
        "../envs/report.yml"
    params:
        dir=REPORT_STEP,
    shell:
        """
        multiqc \
            --title bowtie2_human \
            --force \
            --filename bowtie2_human \
            --outdir {params.dir} \
            --dirs \
            --dirs-depth 1 \
            {input.reports} \
        2> {log} 1>&2
        """


rule report_step_bowtie2_chicken:
    """Collect all reports for the bowtie2 step"""
    input:
        reports=[
            BOWTIE2_CHICKEN / f"{sample}.{library}.{report}"
            for sample, library in SAMPLE_LIB
            for report in BAM_REPORTS
        ],
    output:
        html=REPORT_STEP / "bowtie2_chicken.html",
    log:
        REPORT_STEP / "bowtie2_chicken.log",
    conda:
        "../envs/report.yml"
    params:
        dir=REPORT_STEP,
    shell:
        """
        multiqc \
            --title bowtie2_chicken \
            --force \
            --filename bowtie2_chicken \
            --outdir {params.dir} \
            --dirs \
            --dirs-depth 1 \
            {input.reports} \
        2> {log} 1>&2
        """


rule report_step_bowtie2_mags:
    """Collect all reports for the bowtie2 step"""
    input:
        reports=[
            BOWTIE2_MAGS / f"{sample}.{library}.{report}"
            for sample, library in SAMPLE_LIB
            for report in BAM_REPORTS
        ],
    output:
        html=REPORT_STEP / "bowtie2_mags.html",
    log:
        REPORT_STEP / "bowtie2_mags.log",
    conda:
        "../envs/report.yml"
    params:
        dir=REPORT_STEP,
    shell:
        """
        multiqc \
            --title bowtie2_mags \
            --force \
            --filename bowtie2_mags \
            --outdir {params.dir} \
            --dirs \
            --dirs-depth 1 \
            {input.reports} \
        2> {log} 1>&2
        """


rule report_step:
    """Collect all per step reports for the pipeline"""
    input:
        rules.report_step_reads.output,
        rules.report_step_fastp.output,
        rules.report_step_bowtie2.output,
        rules.report_step_bowtie2_human.output,
        rules.report_step_bowtie2_chicken.output,
        rules.report_step_bowtie2_mags.output,


localrules:
    report_step_reads,
    report_step_fastp,
    report_step_bowtie2,
    report_step_bowtie2_human,
    report_step_bowtie2_chicken,
    report_step_bowtie2_mags,
