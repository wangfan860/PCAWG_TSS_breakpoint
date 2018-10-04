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

def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]

if __name__ == '__main__':
    args = get_args()
    gb = int(args.gene_block)
    input=glob.glob('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/CNV_bedtool/CNV/*merge.csv')
    total_input_list=list(chunks(input, 10))
    input_list = total_input_list[gb]
    ref=pd.read_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/gencode_hg38_56354gene.csv')
    ref1=ref[['gene_id2']]
    ref1.columns=['gene']

    for df_ in input_list:
        print(os.path.basename(df_))
        df_csv=pd.read_csv(df_,index_col=0)
        df_csv.columns=['gene',os.path.basename(df_)]
        ref1=pd.merge(ref1,df_csv, on='gene', how='outer')
        print(os.path.basename(df_))
    ref1.to_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/merged_cnv_segment.{}.csv'.format(gb), index=False)
