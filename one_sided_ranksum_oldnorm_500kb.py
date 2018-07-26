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

    distance_table1=pd.read_csv('/gpfs/data/lyang-lab/users/fan/distance84149x1214_0627.csv',index_col=0)
    expr=pd.read_csv('/gpfs/data/lyang-lab/users/fan/correct_oldnorm_expression84149x1214_0630.csv', index_col=0)

    disease_sv=distance_table1[distance_table1.dcc_project_code == '{}'.format(disease)]
    disease_expr=expr[expr.dcc_project_code == '{}'.format(disease)]
    gene_list=list(expr.columns)[2:]
    table = pd.DataFrame()
    for j in gene_list:
        print(j)
        smaller=disease_expr[disease_sv[j] < 500000][j].astype(float).dropna()
        bigger=disease_expr[(disease_sv[j] >= 500000) | disease_sv[j].isnull()][j].astype(float).dropna()
        results = ss.ranksums(smaller, bigger)
        table = table.append({'cancer_type':'{}'.format(disease),'gene':j, 'pvalue':results[1]/2,'z_statistic':results[0],'samplesize_sv':len(smaller), 'samplesize_wosv':len(bigger)}, ignore_index=True)
    table.to_csv('/gpfs/data/lyang-lab/users/fan/one_sided_oldnorm_500kb/ranksum_{}.csv'.format(disease))
