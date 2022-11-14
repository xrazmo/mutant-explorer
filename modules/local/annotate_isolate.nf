include { PRODIGAL } from "$baseDir/modules/nf-core/prodigal/main"
include { DIAMOND_BLASTX } from "$baseDir/modules/nf-core/diamond/blastx/main"
include { TRIMGALORE } from "$baseDir/modules/nf-core/trimgalore/main"
include { SPADES } from "$baseDir/modules/nf-core/spades/main"
include {CREATD_REF_DB} from "$baseDir/modules/local/utility"
include {SAVE_TO_DB} from "$baseDir/modules/local/utility"

 workflow ANNOTATE_ISOLATES{
    take:
        input_gz
        data_dir
        replace_db
    main:
     
        def ref_db = "${params.data_dir}/db.sqlite3"
        def contigs_dir ="${params.data_dir}/contigs"
        def orfs_dir = "${params.data_dir}/orfs_dir"

        input_ch = Channel.fromFilePairs("${input_gz}/*/*/*{1,2}.fastq.gz");

        CREATD_REF_DB(ref_db,replace_db)
        
        paired_ch = input_ch.map{it->[[id:it[0],single_end:false],it[1]]}
       
        
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
        
        ch_into_db =  PRODIGAL.out.gene_annotations.map{it->['orfs',it[0],it[1]]}
                .concat(DIAMOND_BLASTX.out.txt.map{it->['annotations',it[0],it[1]]})

        file(contigs_dir).mkdir()
        SPADES.out.contigs.map{it-> it[1]}.flatten().collectFile(storeDir:contigs_dir)
        
        file(orfs_dir).mkdir()
        PRODIGAL.out.amino_acid_fasta.map{it-> it[1]}.flatten().collectFile(storeDir:orfs_dir) 
        PRODIGAL.out.nucleotide_fasta.map{it-> it[1]}.flatten().collectFile(storeDir:orfs_dir) 
        
        SAVE_TO_DB(ch_into_db,ref_db)
        
}
