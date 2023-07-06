rule kraken2_assign_one:
    input:
        forward_=FASTP / "{sample}.{library}_1.fq.gz",
        reverse_=FASTP / "{sample}.{library}_2.fq.gz",
        database=features["kraken2_db"],
    output:
        out_gz=KRAKEN2 / "{sample}.{library}.out.gz",
        report=KRAKEN2 / "{sample}.{library}.report",
    log:
        log=KRAKEN2 / "{sample}.{library}.log",
    conda:
        "../envs/kraken2.yml"
    threads: 24
    resources:
        mem_mb=8 * 1024,
        runtime=60,
    shell:
        """
        kraken2 \
            --db {input.database} \
            --threads {threads} \
            --paired \
            --gzip-compressed \
            --output >(pigz -11 > {output.out_gz}) \
            --report {output.report} \
            {input.forward_} \
            {input.reverse_} \
        > {log} 2>&1
        """


rule kraken2_assign_all:
    input:
        [KRAKEN2 / f"{sample}.{library}.report" for sample, library in SAMPLE_LIB],


rule kraken2_report_one:
    input:
        KRAKEN2 / "{sample}.{library}.report",


rule kraken2_report_all:
    input:
        [KRAKEN2 / f"{sample}.{library}.report" for sample, library in SAMPLE_LIB],


rule kraken2:
    input:
        rules.kraken2_assign_all.input,
        rules.kraken2_report_all.input,
