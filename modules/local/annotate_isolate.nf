include { PRODIGAL } from "$baseDir/modules/nf-core/prodigal/main"
include { DIAMOND_BLASTX } from "$baseDir/modules/nf-core/diamond/blastx/main"
include { TRIMGALORE } from "$baseDir/modules/nf-core/trimgalore/main"
include { SPADES } from "$baseDir/modules/nf-core/spades/main"

 workflow ANNOTATE_ISOLATES{
    take:
        input_ch
        data_dir
    
    main:
       
        paired_ch = input_ch.filter{it.seq_rep=='paired'}
                                .map{it->[[id:it.sample.replace('.','-v-'),single_end:false],[it.file1,it.file2]]}
       
        
        TRIMGALORE(paired_ch)  

        SPADES(TRIMGALORE.out.reads.map{it -> [it[0],it[1],[],[]]},[])

        PRODIGAL(SPADES.out.contigs , "gff")

        ch_orfs = PRODIGAL.out.nucleotide_fasta
        
        def diamond_cols = "qseqid sseqid pident qcovhsp scovhsp mismatch gaps evalue bitscore length qlen slen qstart qend sstart send stitle"
        dmnd_ch = Channel.fromPath("${data_dir}/dmnd/*.dmnd")
        ch_orfs.combine(dmnd_ch).combine(Channel.from('txt')).combine(Channel.from(tuple(diamond_cols)))
                    .map{it -> [[id:it[0].id+"__"+it[2].simpleName],it[1],it[2],it[3],it[4]]}
                    .multiMap {it -> fasta: tuple(it[0],it[1])
                                    db: it[2]
                                    ext: it[3]
                                    col: it[4]}
                    .set{diaParam}

        DIAMOND_BLASTX(diaParam.fasta,diaParam.db,diaParam.ext,diaParam.col)
        
        

    emit:
       contigs = SPADES.out.contigs
       gene_annotations = PRODIGAL.out.gene_annotations
       diamond_txt = DIAMOND_BLASTX.out.txt
       nucleotide_fasta = PRODIGAL.out.nucleotide_fasta
       amino_acid_fasta = PRODIGAL.out.amino_acid_fasta
      
}
