rule fastp_trim_one:
    """Run fastp on one library"""
    input:
        fwd=READS / "{sample}_1.fq.gz",
        rev=READS / "{sample}_2.fq.gz",
    output:
        fwd=temp(FASTP / "{sample}_1.fq.gz"),
        rev=temp(FASTP / "{sample}_2.fq.gz"),
        html=FASTP / "{sample}.html",
        json=FASTP / "{sample}_fastp.json",
    log:
        FASTP / "{sample}.log",
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
            --in1 {input.fwd} \
            --in2 {input.rev} \
            --out1 {output.fwd} \
            --out2 {output.rev} \
            --html {output.html} \
            --json {output.json} \
            --compression 1 \
            --verbose \
            --adapter_sequence {params.adapter_fwd} \
            --adapter_sequence_r2 {params.adapter_rev} \
            --thread {threads} \
        """
