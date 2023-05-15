rule reads_link_one:
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


rule reads_link_all:
    """Link all reads in the samples.tsv"""
    input:
        [
            READS / f"{sample}_{end}.fq.gz"
            for sample in SAMPLES
            for end in ["1", "2"]
        ],


rule reads_fastqc_all:
    """Run fastqc on all raw reads"""
    input:
        [
            READS / f"{sample}_{end}_fastqc.{extension}"
            for sample in SAMPLES
            for end in ["1", "2"]
            for extension in ["html", "zip"]
        ],


rule reads:
    """Link all reads and run fastqc on them"""
    input:
        rules.reads_link_all.input,
        rules.reads_fastqc_all.input

localrules:
    reads_link_one
