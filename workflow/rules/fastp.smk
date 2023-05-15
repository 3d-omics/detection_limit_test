rule fastp:
    """Run fastp on one library"""
    input:
        read1=READS / "{sample}_1.fq.gz",
        read2=READS / "{sample}_2.fq.gz"
    output:
        fwd=temp(FASTP / "{sample}_1.fq.gz"),
        rev=temp(FASTP / "{sample}_2.fq.gz"),
        html=FASTP / "{sample}.html",
        json=FASTP / "{sample}_fastp.json"
    log:
        FASTP / "{sample}.log"
    benchmark:
        FASTP / "{sample}.bmk"
    params:
        adapter_fwd=params["fastp"]["adapter_fwd"],
        adapter_rev=params["fastp"]["adapter_rev"]
    threads:
        8
    resources:
        load=1,
        mem_gb=10,
        time='01:00:00'
    conda:
        "../envs/fastp.yml"
    shell:
        """
        fastp \
            --in1 {input.read1} \
            --in2 {input.read2} \
            --out1 {output.read1} \
            --out2 {output.read2} \
            --html {output.html} \
            --json {output.json} \
            --compression 1 \
            --verbose \
            --adapter_sequence {params.adapter_fwd} \
            --adapter_sequence_r2 {params.adapter_rev} \
            --thread {threads} \
        """

rule fastp_all_samples:
    """Collect fastp files"""
    input:
        [
            FASTP / f"{sample}_{end}.fq.gz"
            for sample in SAMPLES
            for end in "1 2".split()
        ],


rule fastp_all:
    input:
        rules.fastp_all_samples.input
