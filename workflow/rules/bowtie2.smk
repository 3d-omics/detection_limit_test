rule bowtie2_build:
    """Build bowtie2 index for the human reference

    Let the script decide to use a small or a large index based on the size of
    the reference genome.
    """
    input:
        reference=REFERENCE / "{genome}.fa.gz",
    output:
        mock=touch(REFERENCE / "{genome}"),
    log:
        BOWTIE2 / "build_{genome}.log",
    benchmark:
        BOWTIE2 / "build_{genome}.bmk"
    conda:
        "../envs/bowtie2.yml"
    params:
        extra=params["bowtie2"]["extra"],
    threads: 8
    resources:
        mem_mb=32 * 1024,
        runtime=24 * 60,
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {params.extra} \
            {input.reference} \
            {output.mock} \
        2> {log} 1>&2
        """


rule bowtie2_map_human_one:
    """Map one library to reference genome using bowtie2

    Output SAM file is piped to samtools sort to generate a CRAM file.
    """
    input:
        forward_=FASTP / "{sample}.{library}_1.fq.gz",
        reverse_=FASTP / "{sample}.{library}_2.fq.gz",
        mock=REFERENCE / "human",
        reference=REFERENCE / "human.fa.gz",
    output:
        cram=protected(BOWTIE2 / "{sample}.{library}.human.cram"),
    log:
        BOWTIE2 / "{sample}.{library}.human.log",
    benchmark:
        BOWTIE2 / "{sample}.{library}.human.bmk"
    params:
        index_prefix=REFERENCE / "human",
        extra=params["bowtie2"]["extra"],
        samtools_mem=params["bowtie2"]["samtools"]["mem_per_thread"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
        unpaired_template=lambda wildcards: BOWTIE2
        / f"{wildcards.sample}.{wildcards.library}.unpaired_%.fq.gz",
    threads: 24
    conda:
        "../envs/bowtie2.yml"
    resources:
        mem_mb=64 * 1024,
        runtime=1440,
    shell:
        """
        (bowtie2 \
            -x {input.mock} \
            -1 {input.forward_} \
            -2 {input.reverse_} \
            --threads {threads} \
            --rg-id '{params.rg_id}' \
            --rg '{params.rg_extra}' \
            {params.extra} \
        | samtools sort \
            -l 9 \
            -M \
            -m {params.samtools_mem} \
            -o {output.cram} \
            --reference {input.reference} \
            --threads {threads} \
        ) 2> {log} 1>&2
        """


rule bowtie2_map_human_all:
    input:
        [BOWTIE2 / f"{sample}.{library}.human.cram" for sample, library in SAMPLE_LIB],


rule bowtie2_extract_nonhuman_one:
    """

    I wish there were a way to do this without having to merge the BAM files
    """
    input:
        cram=BOWTIE2 / "{sample}.{library}.human.cram",
        reference=REFERENCE / "human.fa.gz",
    output:
        forward_=BOWTIE2 / "{sample}.{library}.nonhuman_1.fq.gz",
        reverse_=BOWTIE2 / "{sample}.{library}.nonhuman_2.fq.gz",
    log:
        BOWTIE2 / "{sample}.{library}.nonhuman.log",
    conda:
        "../envs/bowtie2.yml"
    threads: 8
    shell:
        """
        (samtools merge \
            -o /dev/stdout \
            <(samtools view -u -f 4  -F 264 {input.cram}) \
            <(samtools view -u -f 8  -F 260 {input.cram}) \
            <(samtools view -u -f 12 -F 256 {input.cram}) \
        | samtools fastq \
            -1 {output.forward_} \
            -2 {output.reverse_} \
            -0 /dev/null \
            -c 9 \
            --threads {threads}) \
        2> {log} 1>&2
        """


rule bowtie2_extract_nonhuman_all:
    input:
        [
            BOWTIE2 / f"{sample}.{library}.nonhuman_{end}.fq.gz"
            for sample, library in SAMPLE_LIB
            for end in ["1", "2"]
        ],


rule bowtie2_map_chicken_one:
    """Map one library to reference genome using bowtie2

    Output SAM file is piped to samtools sort to generate a CRAM file.
    """
    input:
        forward_=BOWTIE2 / "{sample}.{library}.nonhuman_1.fq.gz",
        reverse_=BOWTIE2 / "{sample}.{library}.nonhuman_2.fq.gz",
        mock=REFERENCE / "chicken",
        reference=REFERENCE / "chicken.fa.gz",
    output:
        cram=protected(BOWTIE2 / "{sample}.{library}.chicken.cram"),
    log:
        BOWTIE2 / "{sample}.{library}.chicken.log",
    benchmark:
        BOWTIE2 / "{sample}.{library}.chicken.bmk"
    params:
        index_prefix=REFERENCE / "chicken",
        extra=params["bowtie2"]["extra"],
        samtools_mem=params["bowtie2"]["samtools"]["mem_per_thread"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    threads: 24
    conda:
        "../envs/bowtie2.yml"
    resources:
        mem_mb=32 * 1024,
        runtime=1440,
    shell:
        """
        (bowtie2 \
            -x {input.mock} \
            -1 {input.forward_} \
            -2 {input.reverse_} \
            --threads {threads} \
            --rg-id '{params.rg_id}' \
            --rg '{params.rg_extra}' \
            {params.extra} \
        | samtools sort \
            -l 9 \
            -M \
            -m {params.samtools_mem} \
            -o {output.cram} \
            --reference {input.reference} \
            --threads {threads} \
        ) 2> {log} 1>&2
        """


rule bowtie2_map_chicken_all:
    input:
        [BOWTIE2 / f"{sample}.{library}.chicken.cram" for sample, library in SAMPLE_LIB],


rule bowtie2_extract_nonchicken_one:
    input:
        cram=BOWTIE2 / "{sample}.{library}.chicken.cram",
        reference=REFERENCE / "chicken.fa.gz",
    output:
        forward_=BOWTIE2 / "{sample}.{library}.nonchicken_1.fq.gz",
        reverse_=BOWTIE2 / "{sample}.{library}.nonchicken_2.fq.gz",
    log:
        BOWTIE2 / "{sample}.{library}.nonchicken.log",
    conda:
        "../envs/bowtie2.yml"
    threads: 8
    shell:
        """
        (samtools merge \
            -o /dev/stdout \
            <(samtools view -u -f 4  -F 264 {input.cram}) \
            <(samtools view -u -f 8  -F 260 {input.cram}) \
            <(samtools view -u -f 12 -F 256 {input.cram}) \
        | samtools fastq \
            -1 {output.forward_} \
            -2 {output.reverse_} \
            -0 /dev/null
            -c 9 \
            --threads {threads}) \
        2> {log} 1>&2
        """


rule bowtie2_extract_nonchicken_all:
    input:
        [
            BOWTIE2 / f"{sample}.{library}.nonchicken_{end}.fq.gz"
            for sample, library in SAMPLE_LIB
            for end in ["1", "2"]
        ],


rule bowtie2_map_mags_one:
    """Map one library to reference genome using bowtie2

    Output SAM file is piped to samtools sort to generate a CRAM file.
    """
    input:
        forward_=BOWTIE2 / "{sample}.{library}.nonchicken_1.fq.gz",
        reverse_=BOWTIE2 / "{sample}.{library}.nonchicken_2.fq.gz",
        mock=REFERENCE / "mags",
        reference=REFERENCE / "mags.fa.gz",
    output:
        cram=protected(BOWTIE2 / "{sample}.{library}.mags.cram"),
    log:
        BOWTIE2 / "{sample}.{library}.mags.log",
    benchmark:
        BOWTIE2 / "{sample}.{library}.mags.bmk"
    params:
        index_prefix=REFERENCE / "mags",
        extra=params["bowtie2"]["extra"],
        samtools_mem=params["bowtie2"]["samtools"]["mem_per_thread"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    threads: 24
    conda:
        "../envs/bowtie2.yml"
    resources:
        mem_mb=64 * 1024,
        runtime=1440,
    shell:
        """
        (bowtie2 \
            -x {input.mock} \
            -1 {input.forward_} \
            -2 {input.reverse_} \
            --threads {threads} \
            --rg-id '{params.rg_id}' \
            --rg '{params.rg_extra}' \
            {params.extra} \
        | samtools sort \
            -l 9 \
            -M \
            -m {params.samtools_mem} \
            -o {output.cram} \
            --reference {input.reference} \
            --threads {threads} \
        ) 2> {log} 1>&2
        """


rule bowtie2_map_mags_all:
    input:
        [BOWTIE2 / f"{sample}.{library}.mags.cram" for sample, library in SAMPLE_LIB],


rule bowtie2_report_all:
    """Generate bowtie2 reports for all libraries:
    - samtools stats
    - samtools flagstats
    - samtools idxstats
    """
    input:
        [
            BOWTIE2 / f"{sample}.{library}.{genome}.{report}"
            for sample, library in SAMPLE_LIB
            for report in BAM_REPORTS
            for genome in GENOMES
        ],


rule bowtie2:
    """Run bowtie2 on all libraries and generate reports"""
    input:
        rules.bowtie2_report_all.input,
        rules.bowtie2_map_mags_all.input,
