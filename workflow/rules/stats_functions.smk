def compose_prefix_for_nonpareil(wildcards):
    return STATS / f"{wildcards.sample}.{wildcards.library}"
