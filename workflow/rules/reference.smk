rule reference_recompress_human:
    """Extract the fasta.gz on config.yaml into genome.fa,gz with bgzip"""
    input:
        fa_gz=features["reference"]["human"],
    output:
        fa_gz=REFERENCE / "human.fa.gz",
    log:
        REFERENCE / "human.log",
    conda:
        "../envs/reference.yml"
    threads: 8
    shell:
        """
        (gzip \
            --decompres \
            --stdout {input.fa_gz} \
        | bgzip \
            --compress-level 9 \
            --threads {threads} \
            --stdout \
            /dev/stdin \
        > {output.fa_gz} \
        ) 2> {log}
        """


rule reference_recompress_chicken:
    """Extract the fasta.gz on config.yaml into genome.fa,gz with bgzip"""
    input:
        fa_gz=features["reference"]["chicken"],
    output:
        fa_gz=REFERENCE / "chicken.fa.gz",
    log:
        REFERENCE / "chicken.log",
    conda:
        "../envs/reference.yml"
    threads: 8
    shell:
        """
        (gzip \
            --decompres \
            --stdout {input.fa_gz} \
        | bgzip \
            --compress-level 9 \
            --threads {threads} \
            --stdout \
            /dev/stdin \
        > {output.fa_gz} \
        ) 2> {log}
        """


rule reference_recompress_mags:
    """Extract the fasta.gz on config.yaml into genome.fa,gz with bgzip"""
    input:
        fa_gz=features["reference"]["mags"],
    output:
        fa_gz=REFERENCE / "mags.fa.gz",
    log:
        REFERENCE / "mags.log",
    conda:
        "../envs/reference.yml"
    threads: 8
    shell:
        """
        (gzip \
            --decompres \
            --stdout {input.fa_gz} \
        | bgzip \
            --compress-level 9 \
            --threads {threads} \
            --stdout \
            /dev/stdin \
        > {output.fa_gz} \
        ) 2> {log}
        """


rule reference:
    """Re-bgzip the reference genome and known variants"""
    input:
        rules.reference_recompress_human.output,
        rules.reference_recompress_chicken.output,
        rules.reference_recompress_mags.output,
