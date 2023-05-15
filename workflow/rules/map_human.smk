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
        "resources/references/human_genome.log"
    benchmark:
        "resources/references/human_genome.bmk"
    shell:
        """
        # Index catted genomes
        bowtie2-build \
        --large-index \
        --threads {threads} \
        {input.genome} {output.index} \
        &> {log}
        """

rule map_human:
    input:
        fwd=FASTP / "{sample}_1.fq.gz",
        rev=FASTP / "{sample}_2.fq.gz",
        genome=features["references"]["human_genome"] + ".rev.2.bt2l"
    output:
        all=MAP_HUMAN / "{sample}_all.bam",
        host=MAP_HUMAN / "{sample}_host.bam",
        fwd=MAP_HUMAN / "{sample}_1.fq.gz",
        rev=MAP_HUMAN / "{sample}_2.fq.gz",
    conda:
        "../envs/map_human.yml"
    threads:
        16
    resources:
        load=1,
        mem_gb=96,
        time='03:00:00'
    log:
        MAP_HUMAN / "{sample}.log",
    benchmark:
        MAP_HUMAN / "{sample}.bmk"
    shell:
        """
        # Map reads to catted reference using Bowtie2
        bowtie2 \
            --time \
            --threads {threads} \
            -x {input.genome} \
            -1 {input.fwd} \
            -2 {input.rev} \
        | samtools view -b -@ {threads} - | samtools sort -@ {threads} -o {output.all} - &&

        # Extract non-host reads (note we're not compressing for nonpareil)
        samtools view -b -f12 -@ {threads} {output.all} \
        | samtools fastq -@ {threads} -1 {output.fwd} -2 {output.rev} - &&

        # Send host reads to BAM
        samtools view -b -F12 -@ {threads} {output.all} \
        | samtools sort -@ {threads} -o {output.host} -
        """

rule map_human_all_samples:
    input:
        [
            MAP_HUMAN / f"{sample}_host.bam"
            for sample in SAMPLES
        ],


rule map_human_all:
    input:
        rules.map_human_all_samples.input
