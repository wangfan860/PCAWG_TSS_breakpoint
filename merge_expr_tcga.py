import pandas as pd
import os
import glob


def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]

input = glob.glob('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/RNA/*.txt')
total_input_list=list(chunks(input, 10))



for i in range(0,len(total_input_list)):
    input_list = total_input_list[i]
    print(i)
    test1=pd.read_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/RNA/7ea980bb-9a4f-4adf-bebd-9c7c56cc5651.FPKM-UQ.txt', sep='\t', header=None)
    test1.columns=['gene','fpkm']
    test=test1[['gene']]
    for df_ in input_list:
        df_csv=pd.read_csv(df_, sep='\t', header=None)
        df_csv.columns=['gene',os.path.basename(df_)]
        test=pd.concat([test,df_csv[[os.path.basename(df_)]]], axis=1)
    test.to_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/merged_fpkm.{}.csv'.format(i), index=False)
