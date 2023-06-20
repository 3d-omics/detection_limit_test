def get_bowtie2_for_library_reports(wildcards):
    sample = wildcards.sample
    library = wildcards.library
    return [
        f"{folder}/{sample}.{library}.{report}"
        for folder in [BOWTIE2_HUMAN, BOWTIE2_CHICKEN, BOWTIE2_MAGS]
        for report in BAM_REPORTS
    ]
