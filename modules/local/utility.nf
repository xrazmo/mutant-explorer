 
 process SAVE_TO_DB{
    tag "$meta.id"
    label "vshort"

    maxForks 1
    input:
    tuple val(type),val(meta), path(intsv)
    val(db_path)

    script:
    """
    python $baseDir/bin/save_to_db.py --db $db_path --in $intsv --type $type
    """
}

process CREATD_REF_DB{
    label "vshort"
    input:
        val db_path
        val replace

    script:
    """
    python $baseDir/bin/create_ref_db.py --db $db_path --replace $replace
    """
}

process SCAFFOLD_CONTIGS{
    tag "$meta.id"
    label "vshort"
    input:
    tuple val(meta), path(contig), path(orfs_aa), path(db)

    output:

    tuple val(meta), path("*.gbk"), emit: gbk
    tuple val(meta), path("*.fa"), emit: fa

    script:
    """
    python $baseDir/bin/prepare_ref_seqs.py --sample ${meta.id} --ann_db $db --contig $contig --protein $orfs_aa
    """
}

process EXTRACT_CONTIGS{
     tag "$meta.id"
    label "vshort"

    conda (params.enable_conda ? 'bioconda::pigz=2.6' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz%3A2.3.4' : 
        'quay.io/biocontainers/pigz' }"
        
    input:
    tuple val(meta), path(contig)

    output:
    tuple val(meta), path("*.fa"), emit: fa

    script:
    """
    pigz -df $contig
    """
}

process ANNOTATE_VCF{
    tag "$meta.id"
    label "vshort"
    input:
    tuple val(meta), path(vcf), path(db)

    output:
    tuple val(meta), path("*.csv"), emit: csv

    script:
    """
    python $baseDir/bin/vcf_parser.py --parent ${meta.parent} --database $db --vcf $vcf
    """
}
process ANNOTATE_VCF_GBK{
    tag "$meta.id"

    input:
    tuple val(meta), path(vcf), path(gbk)

    output:
    tuple val(meta), path("*.csv"), emit: csv

    script:
    """
    python $baseDir/bin/vcf_parser_gbk.py --genbank $gbk --vcf $vcf
    """
}

process MERGE_GENES{
    tag "merge-genes"
    label "vshort"

    input:
    tuple path(gene_dir)

    output:
    tuple path("*.fasta"), emit: fasta

    script:
    """
    python $baseDir/bin/mmseqs_handler.py --task merge --fasta_dir $gene_dir
    """
}

process SAVE_CLUSTERS{
    tag "save_clusters"
    
    input:
     tuple val(meta), path(tsv), path(db)

    
    script:
    """
    python $baseDir/bin/mmseqs_handler.py --task save --labels $tsv --database $db
    """
}