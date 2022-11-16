nextflow.enable.dsl=2

include {ANNOTATE_ISOLATES} from "$baseDir/modules/local/annotate_isolate"
include {EXTRACT_CONTIGS} from "$baseDir/modules/local/utility"
include {INDEX_REFERENCE} from "$baseDir/modules/local/cryptic"
include {FAKE_REMOVE_CONTAM} from "$baseDir/modules/local/cryptic"
include {VARIANT_CALL} from "$baseDir/modules/local/cryptic"

include {SNP_CALL} from "$baseDir/modules/local/snp_calling"
include {MMSEQS_CLUSTERNT} from "$baseDir/modules/local/mmseqs"
include {ANNOTATE_VCF} from "$baseDir/modules/local/utility"
include {MERGE_GENES} from "$baseDir/modules/local/utility"
include {SAVE_CLUSTERS} from "$baseDir/modules/local/utility"


params.data_dir = ""
params.mutant_csv = ""

// DELETE
params.parent_gz = ""

workflow cryptic_index{

    ref_ch = Channel.fromPath("${params.mutant_csv}")
                .splitCsv(header: true)
                .map{it -> [[id:it.parent],it.p_fa]}.unique().view()
   
    def reffa_dir = "${params.data_dir}/ref_fa"
    file(reffa_dir).mkdir()
   
    EXTRACT_CONTIGS(ref_ch)
    INDEX_REFERENCE(EXTRACT_CONTIGS.out.fa,reffa_dir)
}

workflow cryptic_snp{

    db_ch = Channel.fromPath("${params.data_dir}/db.sqlite3")
    def reffa_dir = "${params.data_dir}/ref_fa"
    def vcf_dir = "${params.data_dir}/out_vcf/cryptic"
    file(vcf_dir).mkdir()


    ref_dir_ch = Channel.fromPath("${reffa_dir}/*/*.fai").map{it->[it.parent.name,it.parent]}
    mutant_ch = Channel.fromPath("${params.mutant_csv}").splitCsv(header: true).map{it -> [[id:it.child,parent:it.parent],[it.c_read_1,it.c_read_2]]}
 
    FAKE_REMOVE_CONTAM(mutant_ch)
    
    decomp_ch = FAKE_REMOVE_CONTAM.out.reads.map{it-> [it[0].parent,it[0],it[1]]}
                                  .join(ref_dir_ch).map{it->[it[1],it[2],it[3]]}
    
    VARIANT_CALL(decomp_ch)
    ANNOTATE_VCF(VARIANT_CALL.out.final_vcf.combine(db_ch))

    VARIANT_CALL.out.final_vcf.map{it-> it[1]}.flatten().collectFile(storeDir:vcf_dir)
    ANNOTATE_VCF.out.csv.map{it-> it[1]}.flatten().collectFile(storeDir:vcf_dir)
}

workflow cryptic_snp_single_ref{

    db_ch = Channel.fromPath("${params.data_dir}/db.sqlite3")
    def reffa_dir = "${params.data_dir}/ref_fa"
    def vcf_dir = "${params.data_dir}/out_vcf/cryptic_CP009072"
    file(vcf_dir).mkdir()


    ref_dir_ch = Channel.fromPath("${reffa_dir}/CP009072_1/*.fai").map{it->it.parent}
    mutant_ch = Channel.fromPath("${params.mutant_csv}")
                        .splitCsv(header: true)
                        .map{it -> [[id:it.child,parent:it.parent],[it.c_read_1,it.c_read_2]]}
 
    FAKE_REMOVE_CONTAM(mutant_ch)
    
    decomp_ch = FAKE_REMOVE_CONTAM.out.reads.combine(ref_dir_ch)
    decomp_ch.view()
    VARIANT_CALL(decomp_ch)
    // ANNOTATE_VCF(VARIANT_CALL.out.final_vcf.combine(db_ch))

    VARIANT_CALL.out.final_vcf.map{it-> it[1]}.flatten().collectFile(storeDir:vcf_dir)
    // ANNOTATE_VCF.out.csv.map{it-> it[1]}.flatten().collectFile(storeDir:vcf_dir)
}

workflow annotate{

    // ANNOTATE_ISOLATES(params.parent_gz,params.data_dir,true) 
    def orfs_dir = "${params.data_dir}/orfs_dir"
    db_ch = Channel.fromPath("${params.data_dir}/db.sqlite3")

    MERGE_GENES(orfs_dir)
    MMSEQS_CLUSTERNT(MERGE_GENES.out.fasta.map{it->[[id:it.simpleName],it]},'easy-cluster',0.90)
    SAVE_CLUSTERS(MMSEQS_CLUSTERNT.out.tsv.combine(db_ch))
}

workflow snp{

    db_ch = Channel.fromPath("${params.data_dir}/db.sqlite3")

    mutant_ch = Channel.fromPath("${params.mutant_csv}")
                .splitCsv(header: true)
                .map{it -> [[id:it.child,parent:it.parent],[it.c_read_1,it.c_read_2],it.p_fa]}

    SNP_CALL(mutant_ch,params.data_dir)
    ANNOTATE_VCF(SNP_CALL.out.vcf.combine(db_ch))
    
    def vcf_dir = "${params.data_dir}/out_vcf/bcf"
    file(vcf_dir).mkdir()

    SNP_CALL.out.vcf.map{it-> it[1]}.flatten().collectFile(storeDir:vcf_dir)
    ANNOTATE_VCF.out.csv.map{it-> it[1]}.flatten().collectFile(storeDir:vcf_dir)
}
