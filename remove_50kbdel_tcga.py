#romove deletion smaller than 50kb
import pandas as pd
import glob

manifest=pd.read_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/8990cnv_with_expr.csv')
array_input=manifest.file_name

for i in range(len(array_input)):
    array_file= array_input[i]
    print(str(array_file))
    array = pd.read_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/CNV/{}'.format(array_file), sep='\t')
    array['copy']=2**(1+array.Segment_Mean)
    array.columns=['gdc_aliquot','chromosome','start','end','probe','seg','copy']
    chrom1=[]
    bpt1=[]
    change1=[]
    copy_before=[]
    copy_after=[]
    for index, row in array.drop(array.index[len(array)-1]).iterrows():
        if row['copy']> array.iloc[index +1]['copy']:
            bpt1.append(row['end'])
            chrom1.append(row['chromosome'])
            change1.append('deletion')
            copy_before.append(row['copy'])
            copy_after.append(array.iloc[index +1]['copy'])

        if row['copy']< array.iloc[index +1]['copy']:
            bpt1.append(row['end'])
            chrom1.append(row['chromosome'])
            change1.append('amplification')
            copy_before.append(row['copy'])
            copy_after.append(array.iloc[index +1]['copy'])
    array_break= pd.DataFrame({'chrom':chrom1,'breakpoint':bpt1,'change':change1,'copy_before':copy_before,'copy_after':copy_after})
    barcode=manifest.barcode[i]
    array_break.to_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/remove_50kb/array_{}.csv'.format(barcode),index=False)
#convert to two breakpoints format
input1=glob.glob('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/remove_50kb/array_*.csv')
for i in input1:
    print(i)
    file = pd.read_csv('{}'.format(i))
    if len(file) != 0:
        chrom=[]
        bpt1=[]
        bpt2=[]
        change=[]

        for index, row in file.drop(file.index[len(file)-1]).iterrows():
            if row['chrom'] == file.iloc[index +1]['chrom'] and row['copy_after']== file.iloc[index +1]['copy_before'] and row['change'] == 'deletion' and file.iloc[index +1]['change']=='amplification':
                bpt1.append(row['breakpoint'])
                bpt2.append(file.iloc[index +1]['breakpoint'])
                chrom.append(row['chrom'])
                change.append('deletion')

            if row['chrom'] == file.iloc[index +1]['chrom'] and row['copy_after']== file.iloc[index +1]['copy_before'] and row['change'] == 'amplification' and file.iloc[index +1]['change']=='deletion':
                bpt1.append(row['breakpoint'])
                bpt2.append(file.iloc[index +1]['breakpoint'])
                chrom.append(row['chrom'])
                change.append('amplification')

        file_break= pd.DataFrame({'chrom':chrom,'breakpoint1':bpt1,'breakpoint2':bpt2,'change':change})
        file_break['size']=(file_break.breakpoint2-file_break.breakpoint1).abs()
        remove_50kb= file_break.drop(file_break[(file_break['change']=='deletion')&(file_break['size']<50000)].index)
        remove_50kb.to_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/remove_50kb/{}_2breakpoints_format_removed_50kb.csv'.format(i),index=False)
