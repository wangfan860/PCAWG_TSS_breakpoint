import pandas as pd
import os
import glob


def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]

input=glob.glob('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/CNV_bedtool/CNV/*merge.csv')
total_input_list=list(chunks(input, 10))

for i in list(range(0,len(total_input_list))):
    input_list = total_input_list[i]
    print(i)
    ref=pd.read_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/gencode_hg38_56354gene.csv')
    ref1=ref[['gene_id2']]
    ref1.columns=['gene']
    for df_ in input_list:
        print(os.path.basename(df_))
        df_csv=pd.read_csv(df_,index_col=0)
        df_csv.columns=['gene',os.path.basename(df_)]
        df_csv_nodup=df_csv.drop_duplicates(subset='gene', keep='first')
        ref1=pd.merge(ref1,df_csv_nodup, on='gene', how='outer')
    ref1.to_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/merged_cnv_segment.{}.csv'.format(i), index=False)
##LUAD_driverless
import os
import glob

def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]

input=glob.glob('*merge.csv')
total_input_list=list(chunks(input, 11))

for i in list(range(0,len(total_input_list))):
    input_list = total_input_list[i]
    print(i)
    ref=pd.read_csv('gencode_hg19.csv')
    ref1=ref[['transcript_name']]
    ref1.columns=['gene']
    for df_ in input_list:
        print(os.path.basename(df_))
        df_csv=pd.read_csv(df_,index_col=0)
        df_csv.columns=['gene',os.path.basename(df_)]
        df_csv_nodup=df_csv.drop_duplicates(subset='gene', keep='first')
        ref1=pd.merge(ref1,df_csv_nodup, on='gene', how='outer')
    ref1.to_csv('merged_cnv_segment.{}.csv'.format(i), index=False)
