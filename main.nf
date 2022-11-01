nextflow.enable.dsl=2

include {ANNOTATE_ISOLATES} from "$baseDir/modules/local/annotate_isolate"
include {SCAFFOLD_CONTIGS} from "$baseDir/modules/local/utility"

params.data_dir = ""
params.input_gz = ""


workflow{

    // ANNOTATE_ISOLATES(input_ch,params.data_dir,True)
    contigs_ch = Channel.fromPath("${params.data_dir}/contigs/*.gz").map{it->[[id:it.simpleName],it]}
    proteins_ch = Channel.fromPath("${params.data_dir}/orfs_aa/*.faa").map{it->[[id:it.simpleName],it]}
    db_ch = Channel.fromPath("${params.data_dir}/db.sqlite3")
    
    SCAFFOLD_CONTIGS(contigs_ch.join(proteins_ch).combine(db_ch))
    
    
}
