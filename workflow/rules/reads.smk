rule reads_link_pe:
    """Make a link to the original file, with a prettier name than default"""
    input:
        forward_=get_forward,
        reverse_=get_reverse
    output:
        forward_=temp(READS / "{sample}_1.fq.gz"),
        reverse_=temp(READS / "{sample}_2.fq.gz")
    log:
        READS / "{sample}.log"
    conda:
        "../envs/empty.yml"
    shell:
        """
        ln --symbolic $(readlink --canonicalize {input.forward_}) {output.forward_}
        ln --symbolic $(readlink --canonicalize {input.reverse_}) {output.reverse_}
        """


rule reads_link:
    input:
        [
            READS / f"{sample}_{end}.fq.gz"
            for sample in SAMPLES
            for end in ["1", "2"]
        ],


rule reads_fastqc:
    """Collect fasqtc reports from the reads"""
    input:
        [
            READS / f"{sample}_{end}_fastqc.html"
            for sample in SAMPLES
            for end in ["1", "2"]
        ],