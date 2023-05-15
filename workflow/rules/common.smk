def get_reads(wildcards):
    forward_, reverse_ = samples[
        (samples["sample"] == wildcards.sample)
    ][["forward", "reverse"]].values[0]
    return forward_, reverse_


def get_forward(wildcards):
    return samples[
        (samples["sample"] == wildcards.sample)
    ]["forward"].tolist()[0]


def get_reverse(wildcards):
    return samples[
        (samples["sample"] == wildcards.sample)
    ]["reverse"].tolist()[0]
