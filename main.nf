nextflow.enable.dsl=2

include {ANNOTATE_ISOLATES} from "$baseDir/modules/local/annotate_isolate"
include {CREATD_REF_DB} from "$baseDir/modules/local/db_manager"
include {SAVE_TO_DB} from "$baseDir/modules/local/db_manager"

params.data_dir = "$baseDir/data"
params.input_gz = "$baseDir/gz"
def ref_db = "${params.data_dir}/db.sqlite3"
def contigs_dir ="${params.data_dir}/contigs"

workflow{

    CREATD_REF_DB(ref_db,true)

    input_ch = Channel.fromFilePairs("${params.input_gz}/*/*{1,2}.fastq.gz");
    input_ch.view()

    ANNOTATE_ISOLATES(input_ch,params.data_dir)

    ch_into_db = ANNOTATE_ISOLATES.out.gene_annotations.map{it->['orfs',it[0],it[1]]}
                .concat(ANNOTATE_ISOLATES.out.diamond_txt.map{it->['annotations',it[0],it[1]]})

    file(contigs_dir).mkdir()
    ANNOTATE_ISOLATES.out.contigs.map{it-> it[1]}.flatten().collectFile(storeDir:contigs_dir)
    
    SAVE_TO_DB(ch_into_db,ref_db)
 
}