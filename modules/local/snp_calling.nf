include { BOWTIE2_BUILD as bwb } from "$baseDir/modules/nf-core/bowtie2/build/main"
include { SAMTOOLS_FAIDX as smf } from "$baseDir/modules/nf-core/samtools/faidx/main"
include { BOWTIE2_ALIGN as align} from "$baseDir/modules/nf-core/bowtie2/align/main"
include { EXTRACT_CONTIGS as eco} from "$baseDir/modules/local/utility"


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

        bwb(ref_ch)
        eco(ref_ch)
        smf(eco.out.fa)
        align(mutant_ch, bwb.out.index, false, true)

        ref_indexed_ch = eco.out.fa.join(smf.out.fai)
        bam_ch = align.out.bam.map{it->[it[0].parent,it[0],it[1]]}
                                  .join(ref_indexed_ch.map{it-> [it[0].id,[it[1],it[2]]]}).map{it->[it[1],it[2],it[3]]}
        // bam_ch.view()
        BCFTOOLS_MPILEUP(bam_ch)

    emit:
        vcf = BCFTOOLS_MPILEUP.out.vcf

}