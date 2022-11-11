
import getopt, os, sys
import pandas as pd
import sqlite3 as lit

def parse_vcf(sample,refdb,vcf):
   
    orfs = __fetch_orfs(sample,refdb)
    mdf = __vcf2df(vcf)    
    tmp_lst= []
    for _,row in mdf.iterrows():
        r = __locate_muation(row["#CHROM"],orfs,row["POS"])
        r["contig_id"] = row["#CHROM"]
        r["genome_position"] = row["POS"]
        r["ref_nt"] = row["REF"]
        r["alt_nt"] = row["ALT"]
        r["filter"] = row["FILTER"]
        tmp_lst.append(r)

    df = pd.DataFrame(tmp_lst)
    sname = os.path.basename(vcf).split('.')[0]
    ofname = f'{sname}.annotated_mutations.csv'
    df.to_csv(ofname)
    return ofname

def __vcf2df(vcf):
     # find the comments lines
    comment_lines = 0
    with open(vcf) as hdl:
        for ll in hdl:    
            if ll.startswith('##'):
               comment_lines +=1
            else:
                break

    df = pd.read_csv(vcf,sep='\t',header=0,skiprows=comment_lines)
    return df        

def __fetch_orfs(sample,refdb):
    conn = lit.connect(refdb)
    cur = conn.cursor()

    cmd = "select orfs.accession,start_index,end_index,strand,sseqid,stitle,pident,scov" \
          " from orfs left join annotations as ann on orfs.accession = ann.orfs_accession" \
          f" where ref_db = 'nr' and orfs.sample = '{sample}';"
    orfs = {}
    for ac,si,ei,snd,ss,st,pi,sc in cur.execute(cmd):
        contig_id= '_'.join(ac.split('_')[0:-1])
        if contig_id not in orfs:
            orfs[contig_id] = []
        ptr = orfs[contig_id]
        ptr.append({"accession":ac,"sidx":si,"eidx":ei,
                      "strand":snd,"sseqid":ss,"stitle":st,"pident":pi,"scov":sc})
    conn.close()
    return orfs

def __locate_muation(id,orfs,pos):
    if id not in orfs:
        return {}
    type = None
    for orf in orfs[id]:
        if orf["strand"] == '+':
            if orf["sidx"]<= pos and orf["eidx"]>= pos: # CDS
                type = 'CDS'
            elif orf["sidx"]-100<= pos and  orf["sidx"]> pos: # promoter     
                type = 'promoter'
            ge_pos = pos - orf["sidx"] + 1  
        else:
            if orf["sidx"]<= pos and orf["eidx"]>= pos: # CDS
                type = 'CDS'
            elif orf["eidx"]+100 >= pos and  orf["eidx"]< pos: # promoter     
                type = 'promoter' 
            
            ge_pos = orf["eidx"]-pos + 1 

        if type != None:
            return {**orf,**{"type":type,"gene_pos":ge_pos}}

    return {}

def main(argv):
    sample =None
    refdb = None
    vcf = None
    
    try:
        opts, _ = getopt.getopt(argv, "p:d:v:",["parent=","database=","vcf="])
    except getopt.GetoptError:
        print(getopt.GetoptError.e)
        sys.exit(2)

    for opt, arg in opts:
        if opt in ["--database","-d"]:
            refdb = arg
        if opt in ["--vcf","-v"]:
            vcf = arg
        if opt in ["--parent","-p"]:
            sample = arg
    
    
    parse_vcf(sample,refdb,vcf)

if __name__ == '__main__':
    # main(sys.argv[1:])
    os.chdir("/crex/proj/snic2022-23-507/private/mutant/ref_dbs/out_vcf/")
    df = pd.read_csv("../input_mutants.csv",header=0)
    for _,row in df.iterrows():
        parse_vcf(row["parent"],"../db.sqlite3",f"{row['id']}.final.vcf")