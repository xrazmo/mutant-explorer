 process MMSEQS_CLUSTERNT{
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::mmseqs2=14.7e284-0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mmseqs2%3A14.7e284--pl5321hf1761c0_0' :
        'quay.io/microbiome-informatics/mmseqs' }"

    input:
        tuple val(meta), path(genes_fa)
        val method
        val similarity
    output:
    tuple val(meta), path('*.tsv'),  emit: tsv

    script:
    def args = task.ext.args ?: ''
    def dbname = "${meta.id}_db"
    def clustname = "${meta.id}_clust"
    def mem = task.memory.toGiga()
    def label_f = "${meta.id}.mmseqs.lbl.tsv"
    """
    mmseqs \\
            $method \\
            $genes_fa \\
            $clustname \\
            tmp \\
            $args \\
            --min-seq-id $similarity \\
            --threads $task.cpus \\
            --split-memory-limit ${mem}G \\
    """
 }