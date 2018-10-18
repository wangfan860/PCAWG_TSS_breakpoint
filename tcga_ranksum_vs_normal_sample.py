'''wilcoxon by cancer type'''

import argparse
import pandas as pd
import scipy.stats as ss

def get_args():
    '''
    Loads the parser
    '''
    # Main parser
    parser = argparse.ArgumentParser(description="wilcoxon by cancer type")
    # Args
    required = parser.add_argument_group("Required input parameters")
    # Metadata from input table
    required.add_argument('-d', '--disease_type', required=True, help='Disease type.')
    return parser.parse_args()

if __name__ == '__main__':
    args = get_args()
    disease = args.disease_type

    distance_table=pd.read_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/final_dis_56268x8990.csv',index_col=0)
    distance_table0=distance_table.reset_index()
    distance_table1=distance_table0.drop(['index'], axis=1)
    part_expr=pd.read_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/new_norm_exp.csv',index_col=0)
    expr=pd.concat([distance_table1[['barcode','project_id']],part_expr],axis=1)
    normal_expr=pd.read_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/101818_norm_final_normal_675x56268_fpkm.csv',index_col=0)

    disease_sv=distance_table1[distance_table1.project_id == '{}'.format(disease)]
    disease_expr=expr[expr.project_id == '{}'.format(disease)]
    normal_disease_expr=normal_expr[normal_expr.project_id == '{}'.format(disease)]

    gene_list=list(expr.columns)[2:]
    table = pd.DataFrame()
    for j in gene_list:
        print(j)
        smaller=disease_expr[disease_sv[j] < 100000][j].astype(float).dropna()
        bigger=normal_disease_expr[j].astype(float).dropna()
        results = ss.ranksums(smaller, bigger)
        table = table.append({'cancer_type':'{}'.format(disease),'gene':j, 'pvalue':results[1]/2,'z_statistic':results[0],'samplesize_sv':len(smaller), 'samplesize_normal':len(bigger)}, ignore_index=True)
    table.to_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/tcga_ranksum_vs_normal_{}.csv'.format(disease))
