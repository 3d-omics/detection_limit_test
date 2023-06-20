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


rule samtools_stats_cram:
    """Compute stats for a cram"""
    input:
        cram="{prefix}/{sample}.{library}.cram",
        crai="{prefix}/{sample}.{library}.cram.crai",
        # reference=REFERENCE / "fa.gz",
    output:
        tsv="{prefix}/{sample}.{library}.stats.tsv",
    log:
        "{prefix}/{sample}.{library}.stats.log",
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools stats \
            {input.cram} \
        > {output.tsv} \
        2> {log}
        """


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
