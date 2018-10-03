import pandas as pd
import os
import glob
import argparse
def get_args():
    '''
    Loads the parser
    '''
    # Main parser
    parser = argparse.ArgumentParser(description="merge tcga fpkm")
    # Args
    required = parser.add_argument_group("Required input parameters")
    # Metadata from input table
    required.add_argument('-g', '--gene_block', required=True, help='file block.')
    return parser.parse_args()

gb = int(args.gene_block)

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

input = glob.glob('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/RNA/*.txt')
total_input_list=list(chunks(input, 10))
input_list = total_input_list[gb]
test1=pd.read_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/RNA/7ea980bb-9a4f-4adf-bebd-9c7c56cc5651.FPKM-UQ.txt', sep='\t', header=None)
test1.columns=['gene','fpkm']
test=test1[['gene']]


for df_ in input_list:
    df_csv=pd.read_csv(df_, sep='\t', header=None)
    df_csv.columns=['gene',os.path.basename(df_)]
    test=pd.concat([test,df_csv[[os.path.basename(df_)]]], axis=1, ignore_index=True)
test.to_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/merged_fpkm.{}.csv'.format(gb), index=False)
