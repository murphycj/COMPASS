
rule all:
    input:
        "AML-59-001_cellAssignments.tsv"


rule test:
    output:
        "AML-59-001_cellAssignments.tsv"
    container:
        "ghcr.io/murphycj/compass:dev"
    script:
        "COMPASS -i data/AML-59-001 -o AML-59-001 --nchains 4 --chainlength 50 --CNV 1"