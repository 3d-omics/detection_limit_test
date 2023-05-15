rule map_human:
    input:
        fwd=FASTP / "{sample}_1.fq.gz",
        rev=FASTP / "{sample}_2.fq.gz",
        genome=features["references"]["human_genome"] + ".rev.2.bt2l"
    output:
        all=INDEX_HUMAN / "{sample}_all.bam",
        host=INDEX_HUMAN / "{sample}_host.bam",
        fwd=INDEX_HUMAN / "{sample}_1.fq.gz",
        rev=INDEX_HUMAN / "{sample}_2.fq.gz",
    conda:
        "../envs/map_human.yml"
    threads:
        16
    resources:
        load=1,
        mem_gb=96,
        time='03:00:00'
    log:
        INDEX_HUMAN / "{sample}.log",
    benchmark:
        INDEX_HUMAN / "{sample}.bmk"
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
