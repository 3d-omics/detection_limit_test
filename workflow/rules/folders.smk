READS = Path("results/reads/")
REFERENCE = Path("results/reference/")
FASTP = Path("results/fastp/")
KRAKEN2 = Path("results/kraken2/")

BOWTIE2 = Path("results/bowtie2/")
BOWTIE2_INDEX = BOWTIE2 / "index"
BOWTIE2_HUMAN = BOWTIE2 / "human"
BOWTIE2_NONHUMAN = BOWTIE2 / "nonhuman"
BOWTIE2_CHICKEN = BOWTIE2 / "chicken"
BOWTIE2_NONCHICKEN = BOWTIE2 / "nonchicken"
BOWTIE2_MAGS = BOWTIE2 / "mags"

STATS = Path("results/stats/")
STATS_NONPAREIL = STATS / "nonpareil"
STATS_SINGLEM = STATS / "singlem"
STATS_COVERM = STATS / "coverm"

REPORT_STEP = Path("reports/by_step/")
REPORT_LIBRARY = Path("reports/by_library/")
