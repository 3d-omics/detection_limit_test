def get_bowtie2_for_library_reports(wildcards):
    """Compose the paths for the bowtie2 reports"""
    sample = wildcards.sample
    library = wildcards.library
    return [
        f"{folder}/{sample}.{library}.{report}"
        for folder in [BOWTIE2_HUMAN, BOWTIE2_CHICKEN, BOWTIE2_MAGS]
        for report in BAM_REPORTS
    ]


def get_kraken2_for_library_reports(wildcards):
    """Compose the paths for the kraken2 reports"""
    sample = wildcards.sample
    library = wildcards.library
    return [
        f"{KRAKEN2}/{kraken2_db}/{sample}.{library}.report"
        for kraken2_db in KRAKEN2_DBS
    ]
