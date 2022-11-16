import getopt, os, sys
import pandas as pd
import sqlite3 as lit
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from glob import glob
from time import time

def collect_sequences(parent_csv,clusters_csv,refdb,fadir):

    pdf = pd.read_csv(parent_csv,header=0)
    cdf = pd.read_csv(clusters_csv,header=0)
    conn = lit.connect(refdb)
    cur = conn.cursor()

    pairs_idx = {}
    
    for i,row in pdf.iterrows():
        pairs_idx[row.parent] = {"index":i,"name":f"{row.parent}__{row.child}"}
        pairs_idx[row.child] = {"index":i,"name":f"{row.parent}__{row.child}"}
    
    fa_seqs = {}
    for fa in glob(os.path.join(fadir,'*.fna')):
        sample = os.path.basename(fa).split('.')[0]
        fa_seqs[sample] = SeqIO.index(fa,'fasta')
    for _,row in cdf.iterrows():
        seqs_groups = [{"name":None,"seqs":[]} for i in range(len(pairs_idx)) ]
        cmd = f"select sample,accession from orfs where cluster = {row.cluster}"
        for sa,ac in cur.execute(cmd):
            if sa not in pairs_idx:
                continue
            if sa not in fa_seqs:    
                print(f'Missing fasta file: {sa}')
                continue
            if ac not in fa_seqs[sa]:
                print(f'Missing ORF file: {ac}')
                continue
            idx = pairs_idx[sa]["index"]
            name = pairs_idx[sa]["name"]
            rec = fa_seqs[sa][ac]   
            seqs_groups[idx]["seqs"].append(SeqRecord(seq=rec.seq,id=f"{sa}_{ac}",description="",name=""))
            seqs_groups[idx]["name"] =f"{row.stitle.strip()}[{row.cluster}]_{name}"

        for r in seqs_groups:
            if len(r["seqs"])>0:
                prefix = __has_compare(r["seqs"])
                if prefix !="identical":
                    with open(f"{prefix}_{r['name']}.fa",'w') as hdl:
                         SeqIO.write(r["seqs"],hdl,'fasta')   

def collect_sequence_all(parent_csv,refdb,fadir):
    pdf = pd.read_csv(parent_csv,header=0)

    conn = lit.connect(refdb)
    cur = conn.cursor()

    pairs_idx = {}
    
    for i,row in pdf.iterrows():
        name = f"{row.parent.split('_')[-1]}__{row.child.split('_')[-1]}"
        pairs_idx[row.parent] = {"index":i,"name":name}
        pairs_idx[row.child] = {"index":i,"name":name}
    
    
    fa_seqs = {}
    for fa in glob(os.path.join(fadir,'*.fna')):
        sample = os.path.basename(fa).split('.')[0]
        fa_seqs[sample] = SeqIO.index(fa,'fasta')

    annoations = {}
    cmd = "SELECT cast(cluster as integer),stitle,sseqid,max(pident*scov) FROM orfs LEFT JOIN annotations AS ann ON orfs.accession = ann.orfs_accession where ref_db = 'nr' group by cluster;"
    for cl,st,sid,_ in cur.execute(cmd):
        #clean title
        st = st.replace(sid,'').strip()
        st = st.replace("MULTISPECIES:",'').strip()
        st = st.replace("/",' ')
        st = st.split('[')[0].strip()
        annoations[cl] = st.replace(' ','-')

    max_clust = cur.execute("select max(cast(cluster as integer)) from orfs").fetchone()[0]
    tic = time()
    for cluster in range(1,max_clust+1):
        
        if time() - tic > 1:
            print(f'\r{cluster} of {max_clust}',end='')
            tic = time()

        seqs_groups = [{"name":None,"seqs":[]} for i in range(len(pairs_idx)) ]
        cmd = f"select sample,accession from orfs where cluster = {cluster}"
        for sa,ac in cur.execute(cmd):
            if sa not in pairs_idx:
                continue
            if sa not in fa_seqs:    
                print(f'Missing fasta file: {sa}')
                continue
            if ac not in fa_seqs[sa]:
                print(f'Missing ORF file: {ac}')
                continue
            idx = pairs_idx[sa]["index"]
            name = pairs_idx[sa]["name"]
            rec = fa_seqs[sa][ac]   
            seqs_groups[idx]["seqs"].append(SeqRecord(seq=rec.seq,id=f"{sa}_{ac}",description="",name=""))
            seqs_groups[idx]["name"] =f"[{cluster}]_{name}"

        for r in seqs_groups:
            if len(r["seqs"])>0:
                prefix = __has_compare(r["seqs"])
                if prefix !="identical":
                    ann = 'unkown'
                    if cluster in annoations: 
                        ann = annoations[cluster]
                    with open(f"{prefix}_{ann}_{r['name']}.fa",'w') as hdl:
                         SeqIO.write(r["seqs"],hdl,'fasta')   

def __has_compare(seqrecords):
    hash_lst = []
    if len(seqrecords)==1:
        return 'single'

    for rec in seqrecords:
        hash_lst.append(hash(str(rec.seq)))

    diff_element = len(set(hash_lst))
    
    if  diff_element == 1 or diff_element == len(seqrecords)/2:
         return "identical"  
    return "diff"


def main(argv):
    parent_csv =None
    refdb = None
    fa_dir = None
    clusters_csv = None
    
    try:
        opts, _ = getopt.getopt(argv, "c:p:d:f:",["clusters_csv=","parent_csv=","database=",'fasta_dir='])
    except getopt.GetoptError:
        print(getopt.GetoptError.e)
        sys.exit(2)

    for opt, arg in opts:
        if opt in ["--database","-d"]:
            refdb = arg
        if opt in ["--parent_csv","-p"]:
            parent_csv = arg
        if opt in ["--clusters_csv","-c"]:
            clusters_csv = arg
        if opt in ["--fasta_dir","-f"]:
            fa_dir = arg
    
    collect_sequences(parent_csv,clusters_csv,refdb,fa_dir)


if __name__ == '__main__':
    # main(sys.argv[1:]) 
    parent_csv="/crex/proj/snic2022-23-507/private/mutant/ref_dbs/input_mutants.csv"
    clusters_csv="/crex/proj/snic2022-23-507/private/mutant/ref_dbs/out_align/selected_clusters.csv"
    refdb="/crex/proj/snic2022-23-507/private/mutant/ref_dbs/db.sqlite3"
    fa_dir="/crex/proj/snic2022-23-507/private/mutant/ref_dbs/orfs_dir"
    os.chdir("/crex/proj/snic2022-23-507/private/mutant/ref_dbs/out_align")
    collect_sequences(parent_csv,clusters_csv,refdb,fa_dir)
    # collect_sequence_all(parent_csv,refdb,fa_dir)