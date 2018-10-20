import pandas as pd
import scipy.stats as ss

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
    mean_normal=normal_disease_expr.mean(axis=0)
    sd_normal=normal_disease_expr.std(axis=0)
    cutoff_normal=mean_normal + sd_normal*3
    cutoff_frame=cutoff_normal.to_frame(name='cutoff')


    gene_list=list(expr.columns)[2:]
    table = pd.DataFrame()
    for j in gene_list:
        print(j)
        mean_ex=cutoff_frame.loc[j,'cutoff']
        smaller=disease_expr[disease_sv[j] < 100000][j].astype(float).dropna()
        s_high=smaller[smaller>= mean_ex]
        s_low=smaller[smaller< mean_ex]
        bigger=disease_expr[(disease_sv[j] >= 100000) | disease_sv[j].isnull()][j].astype(float).dropna()
        b_high=bigger[bigger>= mean_ex]
        b_low=bigger[bigger< mean_ex]
        results = ss.fisher_exact([[len(s_high), len(b_high)], [len(s_low), len(b_low)]])
        table = table.append({'gene':j, 'pvalue':results[1],'oddsratio':results[0]}, ignore_index=True)
    table.to_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/tcga_fisher_{}.csv'.format(disease))
