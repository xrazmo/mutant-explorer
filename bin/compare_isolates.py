import pandas as pd
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from glob import glob


def filter_monosnp(p_csv, m_csv):
    pdf = pd.read_csv(p_csv, header=0, index_col=0)
    mdf = pd.read_csv(m_csv, header=0, index_col=0)
    mp = pd.merge(mdf, pdf, how='left', on=["pos_genome", "ref_nt"], suffixes=("", "_y"))
    mp = mp[mp["filter_y"].isnull()]
    mp["where"] = "mutant"

    pm = pd.merge(pdf, mdf, how='left', on=["pos_genome", "ref_nt"], suffixes=("", "_y"))
    pm = pm[pm["filter_y"].isnull()]
    pm["where"] = "parent"

    diff_df = pd.concat([pm, mp])
    cols = list(pdf.columns)
    diff_df = diff_df[cols + ["where"]]
    sample = os.path.basename(p_csv).split('.')[0]
    diff_df.to_csv(f"{sample}.diff.csv")


def extract_insertions(csv_dir):
    csv_files = glob(os.path.join(csv_dir, '*.P.csv')) + glob(os.path.join(csv_dir, '*.M.csv'))
    csv_files.sort()
    seq_list = []
    for csv in csv_files:
        pp = os.path.basename(csv).split('.')
        id = '_'.join(pp[:2])
        df = pd.read_csv(csv, header=0, index_col=0)
        df['ins'] = df.apply(lambda row: len(row["ref_nt"]) < len(row["alt_nt"]), axis=1)
        tdf = df[df["ins"]].sort_values(by="alt_nt", key=lambda x: x.str.len(), ascending=False)
        for _, row in tdf.iterrows():
            ppseq = row["alt_nt"].split(',')
            ppseq.sort(key=lambda x: len(x), reverse=True)
            if len(ppseq[0]) > 75:
                seq_list.append(SeqRecord(seq=Seq(ppseq[0]), id=f"{id}_{row['pos_genome']}", description='', name=''))

    with open(f"all.insertion.fa", 'w') as hdl:
        SeqIO.write(seq_list, hdl, 'fasta')


if __name__ == '__main__':
    # main(sys.argv[1:])
    wd = "/crex/proj/snic2022-23-507/private/mutant/ref_dbs/out_vcf/cryptic_CP009072/rename"
    os.chdir(wd)

    # for id in ["CR07", "CR52", "CR13", "CR12", "CR26"]:
    #     filter_monosnp(f"{id}.P.csv", f"{id}.M.csv")

    extract_insertions(wd)
