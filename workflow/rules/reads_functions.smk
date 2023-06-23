def get_forward(wildcards):
    """Get the forward read for a given sample and library"""
    return samples[(samples["sample"] == wildcards.sample)]["forward"].tolist()[0]


def get_reverse(wildcards):
    """Get the reverse read for a given sample and library"""
    return samples[(samples["sample"] == wildcards.sample)]["reverse"].tolist()[0]
