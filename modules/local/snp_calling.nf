include { BOWTIE2_BUILD } from "$baseDir/modules/local/bowtie2/build/main"
include { SAMTOOLS_FAIDX } from "$baseDir/modules/nf-core/samtools/faidx/main"
include { BOWTIE2_ALIGN } from "$baseDir/modules/local/bowtie2/align/main"
include { TRIMGALORE } from "$baseDir/modules/nf-core/trimgalore/main"
include { EXTRACT_CONTIGS } from "$baseDir/modules/local/utility"


process BCFTOOLS_MPILEUP {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::bcftools=1.15.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.15.1--h0ea216a_0':
        'quay.io/biocontainers/bcftools:1.15.1--h0ea216a_0' }"

    input:
    tuple val(meta), path(bam), path(fasta)
       

    output:
    tuple val(meta), path("*.vcf"), emit: vcf

    script:
     
    """

    bcftools mpileup --output-type u --fasta-ref ${fasta[0]} $bam | bcftools call -vm -O v -o ${meta.id}.vcf

    """
}

workflow SNP_CALL{
    take:
        parent_child_ch
        out_dir
    main:
        ref_ch = parent_child_ch.map{it->[[id:it[0].parent],it[2]]}
        mutant_ch = parent_child_ch.map{it->[[id:it[0].id,parent:it[0].parent,single_end:false],it[1]]}

        BOWTIE2_BUILD(ref_ch)
        // EXTRACT_CONTIGS(ref_ch)
        // SAMTOOLS_FAIDX(EXTRACT_CONTIGS.out.fa)
        indexed_ch = BOWTIE2_BUILD.out.index
        mutant_ch.map{it->[it[0].parent,it[0],it[1]]}.join(indexed_ch.map{it->[it[0].id,it[1]]}).view()
        // TRIMGALORE(mutant_ch) 

    //     BOWTIE2_ALIGN(TRIMGALORE.out.reads, BOWTIE2_BUILD.out.index, false, true)

    //     ref_indexed_ch = EXTRACT_CONTIGS.out.fa.join(SAMTOOLS_FAIDX.out.fai)
    //     bam_ch = BOWTIE2_ALIGN.out.bam.map{it->[it[0].parent,it[0],it[1]]}
    //                               .join(ref_indexed_ch.map{it-> [it[0].id,[it[1],it[2]]]}).map{it->[it[1],it[2],it[3]]}
    //     // bam_ch.view()
    //     BCFTOOLS_MPILEUP(bam_ch)

    emit:
        // vcf = BCFTOOLS_MPILEUP.out.vcf
        vcf = BOWTIE2_BUILD.out.index

}