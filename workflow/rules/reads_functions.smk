def get_forward(wildcards):
    return samples[(samples["sample"] == wildcards.sample)]["forward"].tolist()[0]


def get_reverse(wildcards):
    return samples[(samples["sample"] == wildcards.sample)]["reverse"].tolist()[0]
