import getopt, os, sys
import pandas as pd
import sqlite3 as lit
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from glob import glob

def merge_fa(fa_dir):
    fas = glob(os.path.join(fa_dir,'*.fna'))
    records = []
    for ff in fas:
        sample = os.path.basename(ff).split('.')[0]
        with open(ff) as hdl:
            for rec in SeqIO.parse(hdl,'fasta'):
                records.append(SeqRecord(seq=rec.seq,id=f"{sample}__{rec.id}",description="",name=""))
    
    with open('merged_genes.fasta','w') as hdl:
        SeqIO.write(records,hdl,'fasta')

def save_cluster(in_tsv, refdb):
    conn = lit.connect(refdb)
    cur = conn.cursor()
    try:
       cur.execute('ALTER TABLE orfs ADD cluster VARCHAR') 
       conn.commit()
    except:
        pass

    GetIds = lambda x: x.split('__')
    cluster_id = 0
    prog = 0
    with open(in_tsv) as hdl:
        prv_seed = ''
        for ll in hdl:
            ll = ll.strip('\r\n')
            seed, memeber = ll.split('\t')
            # ss,so = GetIds(seed)
            ms,mo = GetIds(memeber)
            
            if prv_seed != seed:
                cluster_id += 1
            prv_seed = seed
            cur.execute(f"UPDATE orfs SET cluster=? where sample = ? and accession = ?",(cluster_id,ms,mo))
            
            prog += 1
            if prog%10==1:
                print(prog)
    print(prog)
    
    conn.commit()
    conn.close()

  
def main(argv):
    task =None
    refdb = None
    fa_dir = None
    labels = None
    
    try:
        opts, _ = getopt.getopt(argv, "t:d:l:f:",["task=","database=","labels=",'fasta_dir='])
    except getopt.GetoptError:
        print(getopt.GetoptError.e)
        sys.exit(2)

    for opt, arg in opts:
        if opt in ["--database","-d"]:
            refdb = arg
        if opt in ["--labels","-l"]:
            labels = arg
        if opt in ["--task","-t"]:
            task = arg
        if opt in ["--fasta_dir","-f"]:
            fa_dir = arg
    
    if task=='merge':
        merge_fa(fa_dir)
    elif task=="save":
        save_cluster(labels,refdb)
    

if __name__ == '__main__':
    main(sys.argv[1:])