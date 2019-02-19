#romove deletion smaller than 50kb
import pandas as pd
import glob

manifest=pd.read_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/8990cnv_with_expr.csv')
array_input=manifest.file_name

for i in range(len(array_input)):
    array_file= array_input[i]
    print(str(array_file))
    array = pd.read_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/CNV/{}'.format(array_file), sep='\t')
    array.columns=['gdc_aliquot','chromosome','start','end','probe','seg']
    chrom1=[]
    bpt1=[]
    bpt2=[]
    gap_size=[]
    for index, row in array.drop(array.index[len(array)-1]).iterrows():
        if row['chromosome']== array.iloc[index +1]['chromosome']:
            if row['end']!= array.iloc[index +1]['start']:
                gap= int(array.iloc[index +1]['start']) - int(row['end'])
                bpt1.append(row['end'])
                bpt2.append(array.iloc[index +1]['start'])
                chrom1.append(row['chromosome'])
                gap_size.append(gap)

    array_break= pd.DataFrame({'chrom':chrom1,'breakpoint1':bpt1,'breakpoint2':bpt2,'gap_size':gap_size})
    barcode=manifest.barcode[i]
    array_break.to_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/gap_size_distribution/gap_size_{}.csv'.format(barcode),index=False)
#convert to two breakpoints format
input1=glob.glob('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/gap_size_distribution/gap_size_*.csv')
appended_data = []
for infile in input1:
    print(str(infile))
    data = pd.read_csv(infile)
    appended_data.append(data)
appended_data = pd.concat(appended_data, axis=0)
appended_data1 = appended_data[appended_data.chrom=='X']
appended_data2 = appended_data[appended_data.chrom!='X']
appended_data2 = appended_data2.apply(pd.to_numeric)
new=appended_data1.append(appended_data2)
new.to_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/gap_size_distribution/all_gap_size.csv',index=False)
