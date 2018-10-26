import pandas as pd
import scipy.stats as ss
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import argparse

'''plot scatter plot using tcga data'''

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
    required.add_argument('-g', '--gene_name', required=True, help='Gene name.')
    return parser.parse_args()

if __name__ == '__main__':
    args = get_args()
    disease = args.disease_type
    gene_id = args.gene_name

    distance_table=pd.read_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/final_dis_56268x8990.csv',index_col=0)
    distance_table0=distance_table.reset_index()
    distance_table1=distance_table0.drop(['index'], axis=1)
    part_expr=pd.read_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/new_norm_exp.csv',index_col=0)
    expr=pd.concat([distance_table1[['barcode','project_id']],part_expr],axis=1)
    normal_expr=pd.read_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/101818_norm_final_normal_675x56268_fpkm.csv',index_col=0)
    map1=pd.read_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/gencode_hg38_56354gene.csv')

    disease_sv=distance_table1[distance_table1.project_id == '{}'.format(disease)]
    disease_expr=expr[expr.project_id == '{}'.format(disease)]
    normal_disease_expr=normal_expr[normal_expr.project_id == '{}'.format(disease)]
    mean_normal=normal_disease_expr.mean(axis=0)
    sd_normal=normal_disease_expr.std(axis=0)
    cutoff_normal=mean_normal + sd_normal*3
    cutoff_frame=cutoff_normal.to_frame(name='cutoff')

    with PdfPages('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/tcga_scatter_{}.pdf'.format(gene_id)) as pdf:
        gene=gene_id
        official_gene=map1[map1.gene_id2==str(gene)]['gene_name'].values[0]
        mean_ex=cutoff_frame.loc[gene,'cutoff']
        gene_dis=disease_sv[gene]
        gene_ex=disease_expr[gene]
        gene_dis1=gene_dis.fillna(1000000)
        fig = plt.figure()
        ax = plt.gca()
        ax.grid(False)
        ax.scatter(gene_dis1 ,gene_ex , c='red', alpha=0.5, edgecolors='red')
        plt.axvline(100000, color='black',linestyle=':')
        plt.axhline(mean_ex, color='black',linestyle=':')
        
        ax.set_xscale('log')

        plt.xlabel('SV breakpoint to TSS')
        plt.ylabel('normalized expression (FPKM-UQ)')
        format_list = [str(official_gene),str(disease)]
        plt.title('{} in {}'.format(*format_list))
        pdf.savefig()
