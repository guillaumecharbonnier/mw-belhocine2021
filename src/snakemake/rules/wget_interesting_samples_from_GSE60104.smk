INTERESTING_SAMPLES_FROM_GSE60104 = [
    "GSM1464990_20140211_733",
    "GSM1464999_20140509_1283",
    "GSM1465003_20140509_1287",
    "GSM1465019_20140211_746"
]

rule wget_interesting_samples_from_GSE60104:
    output:
        expand("out/wget_interesting_samples_from_GSE60104/{sample}.spikein.hg19.bedgraph.tdf", sample = INTERESTING_SAMPLES_FROM_GSE60104)
    shell:
        "src/bash/wget_interesting_samples_from_GSE60104.sh"

