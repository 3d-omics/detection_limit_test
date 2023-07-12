rule crai:
    """Generate a cram index"""
    input:
        cram="{prefix}/{sample}.{library}.cram",
    output:
        crai="{prefix}/{sample}.{library}.cram.crai",
    log:
        "{prefix}/{sample}.{library}.cram.crai.log",
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools index {input} 2> {log} 1>&2"


rule samtools_flagstats_cram:
    """Compute flagstats for a cram"""
    input:
        cram="{prefix}/{sample}.{library}.cram",
        crai="{prefix}/{sample}.{library}.cram.crai",
    output:
        txt="{prefix}/{sample}.{library}.flagstats.txt",
    log:
        "{prefix}/{sample}.{library}.flagstats.log",
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools flagstats {input.cram} > {output.txt} 2> {log}"


rule samtools_idxstats_cram:
    """Compute idxstats for a cram"""
    input:
        cram="{prefix}/{sample}.{library}.cram",
        crai="{prefix}/{sample}.{library}.cram.crai",
    output:
        tsv="{prefix}/{sample}.{library}.idxstats.tsv",
    log:
        "{prefix}/{sample}.{library}.idxstats.log",
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools idxstats {input.cram} > {output.tsv} 2> {log}"


rule samtools_stats_cram:
    """Compute stats for a cram"""
    input:
        cram=BOWTIE2 / "{genome}/{sample}.{library}.cram",
        crai=BOWTIE2 / "{genome}/{sample}.{library}.cram.crai",
        reference=REFERENCE / "{genome}.fa.gz",
    output:
        tsv=BOWTIE2 / "{genome}/{sample}.{library}.stats.tsv",
    log:
        BOWTIE2 / "{genome}/{sample}.{library}.stats.log",
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools stats \
            --reference {input.reference} \
            {input.cram} \
        > {output.tsv} \
        2> {log}
        """
