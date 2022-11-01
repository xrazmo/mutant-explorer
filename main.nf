nextflow.enable.dsl=2

include {ANNOTATE_ISOLATES} from "$baseDir/modules/local/annotate_isolate"

params.data_dir = "$baseDir/data"
params.input_gz = "$baseDir/gz"


workflow{

    // ANNOTATE_ISOLATES(input_ch,params.data_dir,True)
 

    
}
