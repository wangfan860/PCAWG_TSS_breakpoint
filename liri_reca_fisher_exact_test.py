import pandas as pd
import scipy.stats as ss

distance_table1=pd.read_csv('/gpfs/data/lyang-lab/users/fan/distance84149x1214_0627.csv',index_col=0)
expr=pd.read_csv('/gpfs/data/lyang-lab/users/fan/correct_oldnorm_expression84149x1214_0630.csv', index_col=0)
cutoff_reca=pd.read_csv('/gpfs/data/lyang-lab/users/fan/fisher_liri_reca/reca_cutoff.csv',index_col=0)
cutoff_liri=pd.read_csv('/gpfs/data/lyang-lab/users/fan/fisher_liri_reca/liri_cutoff.csv',index_col=0)
##
disease_sv=distance_table1[distance_table1.dcc_project_code == 'LIRI']
disease_expr=expr[expr.dcc_project_code == 'LIRI']
gene_list=list(expr.columns)[2:]
table = pd.DataFrame()
for j in gene_list:
    print(j)
    mean_ex=cutoff_liri.loc[j,'liri_cutoff']
    smaller=disease_expr[disease_sv[j] < 100000][j].astype(float).dropna()
    s_high=smaller[smaller>= mean_ex]
    s_low=smaller[smaller< mean_ex]
    bigger=disease_expr[(disease_sv[j] >= 100000) | disease_sv[j].isnull()][j].astype(float).dropna()
    b_high=bigger[bigger>= mean_ex]
    b_low=bigger[bigger< mean_ex]
    results = ss.fisher_exact([[len(s_high), len(b_high)], [len(s_low), len(b_low)]])
    table = table.append({'gene':j, 'pvalue':results[1]}, ignore_index=True)
table.to_csv('/gpfs/data/lyang-lab/users/fan/fisher_liri_reca/fisher_liri.csv')
##
disease_sv=distance_table1[distance_table1.dcc_project_code == 'RECA']
disease_expr=expr[expr.dcc_project_code == 'RECA']
gene_list=list(expr.columns)[2:]
table = pd.DataFrame()
for j in gene_list:
    print(j)
    mean_ex=cutoff_reca.loc[j,'reca_cutoff']
    smaller=disease_expr[disease_sv[j] < 100000][j].astype(float).dropna()
    s_high=smaller[smaller>= mean_ex]
    s_low=smaller[smaller< mean_ex]
    bigger=disease_expr[(disease_sv[j] >= 100000) | disease_sv[j].isnull()][j].astype(float).dropna()
    b_high=bigger[bigger>= mean_ex]
    b_low=bigger[bigger< mean_ex]
    results = ss.fisher_exact([[len(s_high), len(b_high)], [len(s_low), len(b_low)]])
    table = table.append({'gene':j, 'pvalue':results[1]}, ignore_index=True)
table.to_csv('/gpfs/data/lyang-lab/users/fan/fisher_liri_reca/fisher_reca.csv')
