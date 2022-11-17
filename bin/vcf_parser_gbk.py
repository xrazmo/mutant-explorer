
import getopt, os, sys, re
import Gemucator as ge
import pandas as pd
from time import time

def parse_vcf(genbank,vcf):
    
    ref_genome =  ge.gemucator(genbank_file=genbank)
    
    mdf = __vcf2df(vcf)

    col_order = ["pos_genome","filter","ref_nt","ref_nt2","alt_nt","gene_name","protein_id",
                "product","pos_gene","ref_aa","alt_aa","strand"]
    tmp_lst = []
    tic = time()
    for i,row in mdf.iterrows():
        
        if time()- tic>1:
            print(f"\r{i}/{len(mdf)}",end='')
            tic = time()
            # if len(tmp_lst)>100:
            #     break

        rlen = len(row["REF"])
        alen = len(row["ALT"])
        if(rlen==1 and alen==1):
            r = ref_genome.identify_gene(row["POS"]-1,row["ALT"])
        else:
            r = ref_genome.identify_gene(row["POS"]-1)
        tmp_lst.append({**{"pos_genome":row["POS"],"ref_nt":row["REF"],"alt_nt":row["ALT"],"filter":row["FILTER"]},**r})

    df = pd.DataFrame(tmp_lst)
    sname = os.path.basename(vcf).split('.')[0]
    ofname = f'{sname}.annotated_mutations.csv'
    df = df[col_order]
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

def main(argv):
    genbank = None
    vcf = None
    
    try:
        opts, _ = getopt.getopt(argv, "g:v:",["genbank=","vcf="])
    except getopt.GetoptError:
        print(getopt.GetoptError.e)
        sys.exit(2)

    for opt, arg in opts:
       
        if opt in ["--genbank","-g"]:
            genbank = arg
        if opt in ["--vcf","-v"]:
            vcf = arg
    
    parse_vcf(genbank,vcf)
    

if __name__ == '__main__':
    main(sys.argv[1:])
    # os.chdir("/crex/proj/snic2022-23-507/private/mutant/ref_dbs/out_vcf/cryptic_CP009072")
    # gbk = "/crex/proj/snic2022-23-507/private/mutant/ref_dbs/ref_fa/CP009072.1.gbk"
    # parse_vcf(gbk,"220822_M7.final.vcf")