
import getopt, os, sys, re
import Gemucator as ge
import pandas as pd

def parse_vcf(genbank,vcf):
    ref_genome =  ge.gemucator(genbank_file=genbank)
    mdf = pd.read_csv(vcf,sep='\t',header=0,skiprows=18)  
    tmp_lst = []
    for _,row in mdf.iterrows():
        rlen = len(row["REF"])
        alen = len(row["ALT"])
        r = {}
        
        if(rlen==1 and alen==1):
           r=ref_genome.identify_gene(row["POS"]-1,row["ALT"])
           if(len(r)>0):
            if r["ref_base"] == row["REF"].lower():
                # coding region or promoter?
                alt = [r['alt_residue'], row["ALT"].lower()][r["gene_pos"]<0]
                r["grammar"] = f"{r['gene_name']}@{r['ref_residue']}{r['gene_pos']}{alt}"
                r["grammar_nt"] = f"{r['gene_name']}@{row['REF'].lower()}{r['gene_pos']}{row['ALT'].lower()}"
                r["genome_pos"]=row["POS"]
            else:
                print("ERROR: POINT")
                print(row["POS"],r)
                r={}
        elif(rlen>alen):# Deletion
            r = ref_genome.identify_gene(row["POS"]-1)
            if(len(r)>0):
                if r["ref_base"] == row["REF"][0].lower():
                   r["grammar"] = f"{r['gene_name']}@{r['gene_pos']}_del_{row['REF'].lower()}"
                   r["genome_pos"]=row["POS"]
                    # print(row["POS"],grammar)       
                else:
                    print("ERROR:DEL")
                    print(row["POS"],r)
                    r={}

        elif(rlen<alen):# Insertion
            r = ref_genome.identify_gene(row["POS"]-1)
            
            if(len(r)>0):
                if r["ref_base"] == row["REF"].lower():
                    r["grammar"] = f"{r['gene_name']}@{r['gene_pos']}_ins_{row['ALT'].lower()}"
                    r["genome_pos"]=row["POS"]
                    # print(row["POS"],grammar)       
                else:
                    print("ERROR:INS")
                    print(row["POS"],r)
                    r={}

        if len(r)>0:
           tmp_lst.append(r)     

    sname = os.path.basename(vcf).split('.')[0]   
    df = pd.DataFrame(tmp_lst)
    ofname = f'{sname}.grammars.csv'
    df.to_csv(ofname)
    return ofname

def parse_grammars(catalog,grammars):

    cdf = pd.read_csv(catalog,header=0)
    cdf["gene"] = cdf['MUTATION'].apply(lambda x:x.split('@')[0])
    cataloged_genes = set(cdf["gene"].tolist())
    
    gdf = pd.read_csv(grammars,header=0,index_col=0)
    gdf["gene"] = gdf['grammar'].apply(lambda x:x.split('@')[0])
    gdf = gdf[gdf['gene'].isin(cataloged_genes)]
   
    for _,crow in cdf.iterrows():
        grammars = crow["MUTATION"].split('&')
        for gr in grammars:
            cdec = decode(gr)
            
            if(len(cdec)==0):
                # print(gr,crow["PREDICTION"])
                continue

            for _,grow in gdf.iterrows():
                gdec = decode(grow["grammar"])
                if(len(gdec)==0): continue      
                
                if(verify(gdec,cdec)):
                    print(gr,crow["PREDICTION"],grow["grammar"],crow["DRUG"])

def decode(grammar):
    out={}
    patterns = {"mutation":'(.*)@([a-z|!]+)(-?[0-9]+)([a-z|!|\*|\-\*|\?|=]+)',
                "indel1":  '(.*)@(-?[0-9]+)\_([del|ins]+)\_([a-z]+)',
                "indel2": '(.*)@(-?[0-9|\*]+)\_([indel]+)' }

    for t,ptr in patterns.items():
        group = re.findall(ptr,grammar,flags=re.IGNORECASE)
        if(len(group)>0):
            if t=="mutation":
                out["gene_name"], out["ref"], out["pos"], out["alt"] = group[0]
                out["type"] = t
            elif t=="indel1":
              out["gene_name"], out["pos"], out["type"], out["indel_base"] = group[0]
            elif t=="indel2":
              out["gene_name"], out["pos"], out["type"] = group[0]
    return out

def verify(gdec,cdec):
    if (gdec["gene_name"]== cdec["gene_name"] 
        and gdec["pos"]== cdec["pos"]) :
        if  gdec["type"]== cdec["type"] or (gdec["type"]!='mutation' and cdec["type"]=='indel'):

            if gdec['type'] == "mutation" and (gdec["alt"] == cdec["alt"] or cdec["alt"] in ['?','*'] ):
                return True
            
            # indels
            #stop codons
            #promoters

    return False        
            

def main(argv):
    catalog =None
    genbank = None
    vcf = None
    
    try:
        opts, _ = getopt.getopt(argv, "c:g:v:",["catalog=","genbank=","vcf="])
    except getopt.GetoptError:
        print(getopt.GetoptError.e)
        sys.exit(2)

    for opt, arg in opts:
        if opt in ["--catalog","-c"]:
            catalog = arg
        if opt in ["--genbank","-g"]:
            genbank = arg
        if opt in ["--vcf","-v"]:
            vcf = arg
    
    grammars_csv = parse_vcf(genbank,vcf)
    #parse_grammars(catalog,'ERR4810943.grammars.csv')

if __name__ == '__main__':
    main(sys.argv[1:])