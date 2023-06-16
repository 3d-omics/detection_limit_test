rule crai:
    """Generate a cram index"""
    input:
        cram="{prefix}/{sample}.{library}.{genome}.cram",
    output:
        crai="{prefix}/{sample}.{library}.{genome}.cram.crai",
    log:
        "{prefix}/{sample}.{library}.{genome}.cram.crai.log",
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools index {input} 2> {log} 1>&2"


rule samtools_stats_cram:
    """Compute stats for a cram"""
    input:
        cram="{prefix}/{sample}.{library}.{genome}.cram",
        crai="{prefix}/{sample}.{library}.{genome}.cram.crai",
        reference=REFERENCE / "{genome}.fa.gz",
    output:
        tsv="{prefix}/{sample}.{library}.{genome}.stats.tsv",
    log:
        "{prefix}/{sample}.{library}.{genome}.stats.log",
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


rule samtools_flagstats_cram:
    """Compute flagstats for a cram"""
    input:
        cram="{prefix}/{sample}.{library}.{genome}.cram",
        crai="{prefix}/{sample}.{library}.{genome}.cram.crai",
    output:
        txt="{prefix}/{sample}.{library}.{genome}.flagstats.txt",
    log:
        "{prefix}/{sample}.{library}.{genome}.flagstats.log",
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools flagstats {input.cram} > {output.txt} 2> {log}"


rule samtools_idxstats_cram:
    """Compute idxstats for a cram"""
    input:
        cram="{prefix}/{sample}.{library}.{genome}.cram",
        crai="{prefix}/{sample}.{library}.{genome}.cram.crai",
    output:
        tsv="{prefix}/{sample}.{library}.{genome}.idxstats.tsv",
    log:
        "{prefix}/{sample}.{library}.{genome}.idxstats.log",
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools idxstats {input.cram} > {output.tsv} 2> {log}"
