def compose_prefix_for_nonpareil(wildcards):
    return STATS_NONPAREIL / f"{wildcards.sample}.{wildcards.library}"
