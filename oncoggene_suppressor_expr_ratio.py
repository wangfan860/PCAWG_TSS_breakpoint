'''get expression fold change by cancer type, expr_with_sv divide by expr_w/o_sv, oldnorm for oncogene, newnorm for suppressor'''

import argparse
import pandas as pd
import statistics

def get_args():
    '''
    Loads the parser
    '''
    # Main parser
    parser = argparse.ArgumentParser(description="get expression fold change by cancer type")
    # Args
    required = parser.add_argument_group("Required input parameters")
    # Metadata from input table
    required.add_argument('-d', '--disease_type', required=True, help='Disease type.')
    return parser.parse_args()

if __name__ == '__main__':
    args = get_args()
    disease = args.disease_type

    distance_table1=pd.read_csv('/gpfs/data/lyang-lab/users/fan/distance84149x1214_0627.csv',index_col=0)
    old_expr=pd.read_csv('/gpfs/data/lyang-lab/users/fan/correct_oldnorm_expression84149x1214_0630.csv', index_col=0)
    new_expr=pd.read_csv('/gpfs/data/lyang-lab/users/fan/newnorm_expression84149x1214_0627.csv', index_col=0)
    onco=pd.read_csv('/gpfs/data/lyang-lab/users/fan/top5000_oncogene.csv')
    suppr=pd.read_csv('/gpfs/data/lyang-lab/users/fan/top5000_suppressor.csv')

    oncolist=onco.gene_y.tolist()
    supprlist=suppr.gene_y.tolist()
    oncolist.insert(0,'tumor_wgs_aliquot_id')
    oncolist.insert(0,'dcc_project_code')
    supprlist.insert(0,'tumor_wgs_aliquot_id')
    supprlist.insert(0,'dcc_project_code')

    disease_sv=distance_table1[distance_table1.dcc_project_code == '{}'.format(disease)]
    disease_old_expr=old_expr[old_expr.dcc_project_code == '{}'.format(disease)]
    disease_new_expr=new_expr[new_expr.dcc_project_code == '{}'.format(disease)]

    disease_sv_onco=disease_sv[oncolist]
    disease_old_expr_onco=disease_old_expr[oncolist]

    disease_sv_suppr=disease_sv[supprlist]
    disease_new_expr_suppr=disease_new_expr[supprlist]

    onco_median_list = []
    onco_list=list(disease_old_expr_onco.columns)[2:]
    for j in onco_list:
        print(j)
        bigger=disease_old_expr_onco[(disease_sv_onco[j] >= 100000) | disease_sv_onco[j].isnull()][j].astype(float).dropna()
        median = statistics.median(bigger)
        onco_median_list.append(median)
    old_expr_onco_fc=disease_old_expr_onco[onco.gene_y.tolist()].divide(onco_median_list, axis='columns')
    old_expr_onco_fc1=old_expr_onco_fc.set_index(disease_old_expr_onco.tumor_wgs_aliquot_id).transpose()
    old_expr_onco_fc1.reset_index(drop=True, inplace=True)
    onco_table = pd.concat([onco, old_expr_onco_fc1 ], axis=1)
    onco_table.to_csv('/gpfs/data/lyang-lab/users/fan/expr_fold_change_by_cancer_type/onco_foldchange_{}.csv'.format(disease))

    suppr_median_list = []
    suppr_list=list(disease_new_expr_suppr.columns)[2:]
    for j in suppr_list:
        print(j)
        bigger=disease_new_expr_suppr[(disease_sv_suppr[j] >= 100000) | disease_sv_suppr[j].isnull()][j].astype(float).dropna()
        median = statistics.median(bigger)
        suppr_median_list.append(median)
    new_expr_suppr_fc=disease_new_expr_suppr[suppr.gene_y.tolist()].divide(suppr_median_list, axis='columns')
    new_expr_suppr_fc1=new_expr_suppr_fc.set_index(disease_new_expr_suppr.tumor_wgs_aliquot_id).transpose()
    new_expr_suppr_fc1.reset_index(drop=True, inplace=True)
    suppr_table = pd.concat([suppr, new_expr_suppr_fc1 ], axis=1)
    suppr_table.to_csv('/gpfs/data/lyang-lab/users/fan/expr_fold_change_by_cancer_type/suppr_foldchange_{}.csv'.format(disease))
