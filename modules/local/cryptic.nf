 process map_reads {
    tag "$meta.id"
    label 'process_medium'
    errorStrategy 'ignore'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/iqbal-lab-org/clockwork/releases/download/v0.11.3/clockwork_v0.11.3.img':
        'ghcr.io/iqbal-lab-org/clockwork:latest' }"

    input:
    tuple val(meta), path(reads), val(ref_fa)
       

    output:
    tuple val(meta), path("*.sam"), emit: sam

    script:
     def threads = task.cpus
    """
     clockwork map_reads --threads $threads --unsorted_sam ${meta.id} $ref_fa ${meta.id}.sam ${reads[0]} ${reads[1]}
    """
}
process remove_contam {
    tag "$meta.id"
    label 'process_medium'
   
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/iqbal-lab-org/clockwork/releases/download/v0.11.3/clockwork_v0.11.3.img':
        'ghcr.io/iqbal-lab-org/clockwork:latest' }"

    input:
    tuple val(meta), path(in_sam), path(ref_tsv)
       

    output:
    tuple val(meta), path("*{1,2}.fq.gz"), emit: reads
    tuple val(meta), path("*.counts.tsv"), emit: tsv

    script:
     
    """
     clockwork remove_contam $ref_tsv $in_sam ${meta.id}.decontam.counts.tsv ${meta.id}.decontam_1.fq.gz ${meta.id}.decontam_2.fq.gz
    """
}
process variant_call {
    tag "$meta.id"
    label 'process_medium'
   
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://github.com/iqbal-lab-org/clockwork/releases/download/v0.11.3/clockwork_v0.11.3.img':
        'ghcr.io/iqbal-lab-org/clockwork:latest' }"

    input:
    tuple val(meta), path(reads)
    val h37Rv_dir
       

    output:
    tuple val(meta), path("*final.vcf"), emit: final_vcf
    tuple val(meta), path("*cortex.vcf"),optional:true, emit: cortex_vcf
    tuple val(meta), path("*samtools.vcf"),optional:true, emit: samtools_vcf

    script:
     
    """
     clockwork variant_call_one_sample --sample_name ${meta.id} $h37Rv_dir var_call ${reads[0]} ${reads[1]}
     cp ./var_call/final.vcf ${meta.id}.final.vcf
     cp ./var_call/samtools.vcf ${meta.id}.samtools.vcf
     cp ./var_call/cortex.vcf ${meta.id}.cortex.vcf
    """
}
process predict_DST{
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(vcf), path(catalog), path(ref_pkl) 
    output:
    tuple val(meta), path("*.effects.csv"), emit: effects
    tuple val(meta), path("*.mutations.csv"),optional:true, emit: mutations
    tuple val(meta), path("*.variants.csv"),optional:true, emit: variants

    script:
    
    """
    python ${baseDir}/bin/sp3predict.py --prefix ${meta.id} --vcf_file $vcf --catalogue_file $catalog --genome_object $ref_pkl --progress --ignore_vcf_status --ignore_vcf_filter

    """
}

params.ref_dir = "${baseDir}/ref_db"
params.input_gz = "${baseDir}/gz"
params.prefix = ""
params.input_catalog = "${baseDir}/catalogues"
def vcf_dir = "${baseDir}/out_vcf"
def out_dir = "${baseDir}/out_csv"

workflow{

   reads_ch = Channel.fromFilePairs("${params.input_gz}/*/*_{1,2}.fastq.gz").map{it->[[id:it[0]],it[1]]};
   ref_fa = Channel.fromPath("${params.ref_dir}/Ref.remove_contam/*.fa");

   map_reads(reads_ch.combine(ref_fa)) 
   ref_tsv = Channel.fromPath("${params.ref_dir}/Ref.remove_contam/*.tsv");
   remove_contam(map_reads.out.sam.combine(ref_tsv))

    def h37Rv_dir = "${params.ref_dir}/Ref.H37Rv";

   variant_call(remove_contam.out.reads,h37Rv_dir)
   
   file(vcf_dir).mkdir()
   variant_call.out.final_vcf.map{it-> it[1]}.flatten().collectFile(storeDir:vcf_dir)

    // vcf_ch = Channel.fromPath("${vcf_dir}/*.vcf").map{it-> [[id:it.simpleName],it]}
    vcf_ch = variant_call.out.final_vcf
    catalog_ch = Channel.fromPath("${params.input_catalog}/*/*.csv")
    refpkl_ch = Channel.fromPath("${params.input_catalog}/*/*.gz")
    ch = vcf_ch.combine(catalog_ch).combine(refpkl_ch)
   
    predict_DST(ch)
    file(vcf_dir).mkdir()
    predict_DST.out.effects.map{it->it[1]}.flatten()
        .collectFile(name:"merged.effects.${params.prefix}.csv",keepHeader:true,storeDir:out_dir)

}