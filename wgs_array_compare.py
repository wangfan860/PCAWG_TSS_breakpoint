import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import glob

manifest=pd.read_csv('/Users/Fan/Desktop/wgs_vs_array/manifest.csv')
wgs_input=manifest.tumor_wgs_aliquot_id
array_input=manifest.file_name

#wgs breakpoint
for i in range(len(wgs_input)):
    wgs_file=wgs_input[i]
    wgs = pd.read_csv('~/Desktop/CNV/consensus.20170119.somatic.cna/original/{}'.format(wgs_file), sep='\t')
    chrom=[]
    bpt=[]
    change=[]
    for index, row in wgs.drop(wgs.index[len(wgs)-1]).iterrows():
        if row['end'] +1 == wgs.iloc[index +1]['start'] and row['total_cn']> wgs.iloc[index +1]['total_cn']:
            bpt.append(row['end'])
            chrom.append(row['chromosome'])
            change.append('deletion')
        if row['end'] +1 == wgs.iloc[index +1]['start'] and row['total_cn']< wgs.iloc[index +1]['total_cn']:
            bpt.append(row['end'])
            chrom.append(row['chromosome'])
            change.append('amplification')
    wgs_break= pd.DataFrame({'chrom':chrom,'breakpoint':bpt,'change':change})
    barcode=manifest.barcode[i]
    project=manifest.project_id[i]
    wgs_break.to_csv('~/Desktop/wgs_vs_array/wgs_{}_{}.csv'.format(project,barcode),index=False)
#array breakpoint
    array_file= array_input[i]
    array = pd.read_csv('/Users/Fan/Desktop/wgs_vs_array/{}'.format(array_file), sep='\t')
    chrom1=[]
    bpt1=[]
    change1=[]
    for index, row in array.drop(array.index[len(array)-1]).iterrows():
        if row['copy']> array.iloc[index +1]['copy']:
            bpt1.append(row['end'])
            chrom1.append(row['chromosome'])
            change1.append('deletion')
        if row['copy']< array.iloc[index +1]['copy']:
            bpt1.append(row['end'])
            chrom1.append(row['chromosome'])
            change1.append('amplification')
    array_break= pd.DataFrame({'chrom':chrom1,'breakpoint':bpt1,'change':change1})
    array_break.to_csv('~/Desktop/wgs_vs_array/array_{}_{}.csv'.format(project,barcode),index=False)

    min_abs_distance=[]
    for i in range(len(array_break)):
        select_chr=wgs_break[wgs_break.chrom.astype(str) == str(chrom1[i])]
        distance = select_chr.breakpoint - bpt1[i]
        min_distance = distance.abs().min()
        min_abs_distance.append(min_distance)

    difference= pd.DataFrame({'chrom':chrom1,'breakpoint':bpt1,'change':change1,'difference':min_abs_distance})
    difference.to_csv('~/Desktop/wgs_vs_array/difference_{}_{}.csv'.format(project,barcode),index=False)

#calculate hisogram of difference_import glob
input1=glob.glob('/Users/Fan/Desktop/wgs_vs_array/difference_*.csv')
appended_data = []
for infile in input1:
    data = pd.read_csv(infile)
    appended_data.append(data)
appended_data = pd.concat(appended_data, axis=1)
