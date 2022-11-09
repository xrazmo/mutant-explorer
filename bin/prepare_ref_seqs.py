import os, argparse, gzip
from Bio import SeqIO
from glob import glob
from Bio.SeqRecord import SeqRecord
import sqlite3 as lit
from Bio.SeqFeature import SeqFeature, FeatureLocation 
from Bio.Seq import Seq

def get_genbank(sample,db,contig_f,protein_f):
    sep_len = 50
    scaffold = ''
    offset = 0
    contigs_offset = {} 
    
    proteins = SeqIO.index(protein_f,'fasta')

    with gzip.open(contig_f,'rt') as hdl:
        for rec in SeqIO.parse(hdl,'fasta'):
            contigs_offset[rec.id] = offset
            scaffold = scaffold + str(rec.seq) + 'N'*sep_len
            offset = len(rec.seq) + sep_len - 1

    annotations = __fetch_features(sample,db)

    record = SeqRecord(Seq(scaffold),id=sample,annotations={"molecule_type": "DNA"})

    Get_Contig_ID = lambda x: '_'.join(x.split('_')[0:-1])
    for _,o in  annotations.items():
       
        contig_id = Get_Contig_ID(o["accession"])
        offset = contigs_offset[contig_id]
        
        sf = SeqFeature(FeatureLocation(o["sidx"]+offset-1,o["eidx"]+offset,strand=o["strand"]),type='CDS')
        ptr = o['annotations']
        sf.qualifiers['id'] = [o["accession"]]
        if len(ptr)>0:
            annot_list = __flatten_annotations(ptr)

            if len(annot_list)>1:
                sf.qualifiers['gene'] = [annot_list[1]["product"]]
                sf.qualifiers['db_2'] = [annot_list[1]["ref_db"]]
                sf.qualifiers['idty_cov_2']=[f"{annot_list[1]['identity']}%-{annot_list[1]['coverage']}%"]
            
            sf.qualifiers["protein_id"] = [annot_list[0]["protein_id"]]
            sf.qualifiers["product"] = [annot_list[0]["product"]]
            sf.qualifiers["ref_db"] = [annot_list[0]["ref_db"]]
            sf.qualifiers['idty_cov']=[f"{annot_list[0]['identity']}%-{annot_list[0]['coverage']}%"]
            
        sf.qualifiers['translation'] = [proteins[o["accession"]].seq]

        record.features.append(sf)

    with open(f'{sample}.gbk','w') as hdl:
        SeqIO.write(record,hdl,'genbank')

    with open(f'{sample}.fa','w') as hdl:
        SeqIO.write(record,hdl,'fasta')


def __flatten_annotations(annot):
    dbs = ['nr','card','bacmet','vfdb']
    annot_list = []
    for db in dbs:
        if db in annot:
            annot_list.append(annot[db])
    return annot_list

def __fetch_features(sample,db):    
    conn = lit.connect(db)
    cur = conn.cursor()
    
    orfs = {}
    cmd = f"select sample,accession,start_index,end_index,strand from orfs where sample = '{sample}'"
    for s,a,si,ei,sd in cur.execute(cmd):
        orfs[f"{s}__{a}"] = {"sample":s,"accession":a,"sidx":si,"eidx":ei,"strand":[-1,1][sd=='+'],"annotations":{}}
    
    cmd = "select sample,orfs_accession,sseqid,stitle,pident,scov,ref_db " \
          f"from annotations where sample = '{sample}' and pident>70 and scov>80"

    for s,a,sid,sti,pi,sc,r in cur.execute(cmd):
        sti = sti.replace(sid,'').strip()
        key = f"{s}__{a}"
        if key not in orfs:
            print(f'annottion of missing orf, How? {key}')
            continue

        ptr = orfs[key]["annotations"]    
        
        if r not in ptr:
            ptr[r] = {'ref_db':r,'identity':pi,"coverage":sc,'product':sti,"protein_id":sid}
    return orfs

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--sample",required=True,help="name of the sample")
    parser.add_argument("--ann_db",required=True,help="SQLite DB containing all annotations")
    parser.add_argument("--contig",required=True,help="fasta.gz containing contigs")
    parser.add_argument("--protein",required=True,help="multifasta containing the predicted porteins")
    params = parser.parse_args()

    get_genbank(params.sample,params.ann_db,params.contig,params.protein)