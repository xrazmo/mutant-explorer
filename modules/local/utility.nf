 
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
    label "vshort"
    input:
    tuple val(meta), path(contig), path(orfs_aa), path(db)

    output:

    tuple val(meta), path("*.gbk"), emit: gbk
    tuple val(meta), path("*.fa"), emit: fa

    script:
    """
    python ./bin/prepare_ref_seqs.py --sample ${meta.id} --ann_db $db --contig_dir $contig --protein_dir $orfs_aa
    """
}

