rule index_human:
    input:
        genome=features["references"]["human_genome"]
    output:
        index=features["references"]["human_genome"] + ".rev.2.bt2l"
    conda:
        "../envs/index_human.yml"
    threads:
        16
    resources:
        load=1,
        mem_gb=96,
        time='03:00:00'
    log:
        "resources/index_human/{sample}.log",
    benchmark:
        "resources/index_human/{sample}.bmk"
    shell:
        """
        # Index catted genomes
        bowtie2-build \
        --large-index \
        --threads {threads} \
        {input.genome} {output.index} \
        &> {log}
        """
