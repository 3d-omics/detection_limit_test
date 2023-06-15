def compose_prefix_for_nonpareil(wildcards):
    return DIVERSITY / f"{wildcards.sample}.{wildcards.library}"


rule diversity_nonpareil_one:
    """
    Note: Nonpareil only ask for one of the pair-end reads
    Note2: it has to be fastq. The process substitution trick does not work
    """
    input:
        forward_=BOWTIE2 / "{sample}.{library}.nonchicken_1.fq.gz",
    output:
        forward_fq=temp(DIVERSITY / "{sample}.{library}_nonchicken_1.fq"),
        npa=DIVERSITY / "{sample}.{library}.npa",
        npc=DIVERSITY / "{sample}.{library}.npc",
        npl=DIVERSITY / "{sample}.{library}.npl",
        npo=DIVERSITY / "{sample}.{library}.npo",
    log:
        DIVERSITY / "{sample}.{library}_nonchicken_1.log",
    conda:
        "../envs/diversity.yml"
    params:
        prefix=compose_prefix_for_nonpareil,
    shell:
        """
        gzip -dc {input.forward_} > {output.forward_fq} 2> {log}

        nonpareil \
            -s {output.forward_fq} \
            -T kmer \
            -b {params.prefix} \
            -f fastq \
            -t {threads} \
        2>> {log} 1>&2
        """


rule diversity_nonpareil_all:
    input:
        [
            DIVERSITY / f"{sample}.{library}.{extension}"
            for extension in ["npa", "npc", "npl", "npo"]
            for sample, library in SAMPLE_LIB
        ],


rule diversity_singlem_one:
    """
    Note: SingleM asks in the documentation for the raw reads. Here we are
    passing it the non-host and trimmed ones.
    """
    input:
        forward_=BOWTIE2 / "{sample}.{library}.nonchicken_1.fq.gz",
        reverse_=BOWTIE2 / "{sample}.{library}.nonchicken_2.fq.gz",
    output:
        otu_table=DIVERSITY / "{sample}.{library}.otu_table.tsv",
    log:
        DIVERSITY / "{sample}.{library}.singlem.log",
    conda:
        "../envs/diversity.yml"
    threads: 24
    shell:
        """
        singlem pipe \
            --forward {input.forward_} \
            --reverse {input.reverse_} \
            --otu_table {output.otu_table} \
            --threads {threads} \
        2> {log} 1>&2
        """


rule diversity_singlem_all:
    input:
        [
            DIVERSITY / f"{sample}.{library}.otu_table.tsv"
            for sample, library in SAMPLE_LIB
        ],


rule diversity:
    input:
        rules.diversity_nonpareil_all.input,
        rules.diversity_singlem_all.input,
