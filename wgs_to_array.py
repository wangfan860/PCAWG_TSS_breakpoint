import pandas as pd


def Distance_to_tss(wgs_file,array_file):
    wgs = pd.read_csv(wgs_file, sep='\t')
    array = pd.read_csv(array_file,sep='\t',names = ["chromosome", "Start", "End", "probe",'segment','copy'])
    array['chromosome']=array.chromosome.str.replace('chr','')

    wgs_chrom=[]
    wgs_bpt=[]
    for index, row in wgs.drop(wgs.index[len(wgs)-1]).iterrows():
        if row['end'] +1 == wgs.iloc[index +1]['start'] and row['total_cn']!=wgs.iloc[index +1]['total_cn']:
            wgs_bpt.append(row['end'])
            wgs_chrom.append(row['chromosome'])
    gold_std = pd.DataFrame({'chrom':wgs_chrom,'breakpoint':wgs_bpt})

    array_chrom=[]
    array_bpt=[]
    for index, row in array.drop(array.index[len(array)-1]).iterrows():
        if row['copy']!=array.iloc[index +1]['copy']:
            array_bpt.append(row['End'])
            array_chrom.append(row['chromosome'])

    min_abs_distance=[]
    for i in range(len(array_chrom)):
        select_chr=gold_std[gold_std.chrom.astype(str) == str(array_chrom[i])]
        distance = select_chr.breakpoint - array_bpt[i]
        min_distance = distance.abs().min()
        min_abs_distance.append(min_distance)
    return min_abs_distance

wgs_file='/Users/Fan/Desktop/CNV/consensus.20170119.somatic.cna/original/b9d1a64e-d445-4174-a5b4-76dd6ea69419.consensus.20170119.somatic.cna.txt'
array_file='/Users/Fan/Documents/out_TCGA-C5-A1M9.csv'
result= Distance_to_tss(wgs_file,array_file)
result1 = [x for x in result if str(x) != 'nan']

wgs_file2='/Users/Fan/Desktop/CNV/consensus.20170119.somatic.cna/original/9ff21093-58d7-4b69-aade-c242a383ea56.consensus.20170119.somatic.cna.txt'
array_file2='/Users/Fan/Documents/out_TCGA-EK-A2R9.csv'
result2= Distance_to_tss(wgs_file2,array_file2)
result3 = [x for x in result2 if str(x) != 'nan']

wgs_file3='/Users/Fan/Desktop/CNV/consensus.20170119.somatic.cna/original/d8fbb398-d1da-4444-984a-22c8523625da.consensus.20170119.somatic.cna.txt'
array_file3='/Users/Fan/Documents/out_TCGA-A2-A3Y0.csv'
result4= Distance_to_tss(wgs_file3,array_file3)
result5 = [x for x in result4 if str(x) != 'nan']

wgs_file4='/Users/Fan/Desktop/CNV/consensus.20170119.somatic.cna/original/4838b5a9-968c-4178-bffb-3fafe1f6dc09.consensus.20170119.somatic.cna.txt'
array_file4='/Users/Fan/Documents/out_TCGA-DK-A3IL.csv'
result6= Distance_to_tss(wgs_file4,array_file4)
result7 = [x for x in result6 if str(x) != 'nan']

wgs_file5='/Users/Fan/Desktop/CNV/consensus.20170119.somatic.cna/original/9c70688d-6e43-4520-9262-eaae4e4d597d.consensus.20170119.somatic.cna.txt'
array_file5='/Users/Fan/Documents/out_TCGA-BH-A18R.csv'
result8= Distance_to_tss(wgs_file5,array_file5)
result9 = [x for x in result8 if str(x) != 'nan']

import matplotlib.pyplot as plt

plt.hist(final1, bins=10000)
ax = plt.gca()
ax.grid(False)
ax.set_xscale('log')
ax.set_yscale('log')
plt.ylim(0, 1000000)
plt.xlabel('distance between array breakpoint and wgs breakpoint')
plt.ylabel('frequency')
a= ['1bp','10bp','100bp','1kb','10kb','100kb','1mb','10mb','100mb']
#ax.set_xticklabels(a)
plt.title('Delta(array-wgs) histogram')
plt.show()
