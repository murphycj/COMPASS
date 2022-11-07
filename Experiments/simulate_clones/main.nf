samples = Channel.from([
    [2, '0.9,0.1'],
    [2, '0.95,0.05'],
    [2, '0.98,0.02'],
    [2, '0.99,0.01'],
    [2, '0.995,0.005'],
    [2, '0.999,0.001'],
    [3, '0.9,0.05,0.05'],
    [3, '0.95,0.025,0.025'],
    [3, '0.98,0.01,0.01'],
    [3, '0.99,0.005,0.005'],
    [3, '0.995,0.0025,0.0025'],
    [3, '0.999,0.0005,0.0005'],
]).spread([1,2,3,4,5,6,7,8,9,10]).map{[[id: it[0].toString() + '_' + it[1].toString().replace(',', '-') + '_rep' + it[2].toString()], it[0], it[1]]}
// samples.view()

process generate {
    tag "$meta.id"
    storeDir "generate/${prefix}"

    input:
        tuple val(meta), val(nnodes), val(nodeprobconcentration) from samples
    output:
        tuple val(meta), path("${prefix}_variants.csv"), path("${prefix}_regions.csv") into generate_out

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    python /Users/charliemurphy/Desktop/github/forks/COMPASS/Experiments/simulate_clones/generate_synthetic_data_main.py \
        -o ${prefix} \
        -m COMPASS \
        --nSNVs 2 \
        --nCNVs 0 \
        --nCNLOHs 0 \
        --depth 50 \
        --ncells 2500 \
        --nnodes ${nnodes} \
        --nodeprobconcentration ${nodeprobconcentration}
    """
}

process COMPASS {
    tag "$meta.id"
    maxForks 2
    storeDir "COMPASS_out/${prefix}"

    input:
        tuple val(meta), path(variants), path(regions) from generate_out
    output:
        path("${prefix}_cellAssignmentProbs.tsv")
        path("${prefix}_cellAssignments.tsv")
        path("${prefix}_nodes_copynumbers.tsv"), optional:true
        path("${prefix}_nodes_genotypes.tsv")
        path("${prefix}_tree.gv")
        path("${prefix}_tree.json")
        path("${prefix}_tree.png")

    storeDir 'COMPASS_out'
    container 'ghcr.io/murphycj/compass:dev'
    stageInMode 'copy'

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    COMPASS -i ${prefix} -o ${prefix} --nchains 4 --chainlength 5000 --CNV 1
    dot -Tpng -o ${prefix}_tree.png ${prefix}_tree.gv
    """
}


// process 