def compose_prefix_for_nonpareil(wildcards):
    """Compose prefix for nonpareil output files"""
    return STATS_NONPAREIL / f"{wildcards.sample}.{wildcards.library}"
