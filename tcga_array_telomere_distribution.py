import pandas as pd
import glob

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
        file_break.to_csv('{}_2breakpoints_format.txt'.format(i),index=False)
#add telemere distance
array_input=glob.glob('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/remove_50kb/*_2breakpoints_format.txt')
telomere=pd.read_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/telomere_distance_distribution/hg38_telomere.txt',sep='\t')

for i in range(len(array_input)):
    array_file= array_input[i]
    print(str(array_file))
    array = pd.read_csv(array_file)
    bpt1_telo=[]
    bpt2_telo=[]
    for index, row in array.iterrows():
        select = telomere[telomere.chrom == str(row['chrom'])]
        if (list(select.chromStart)[0]< row['breakpoint1']) & (list(select.chromEnd)[0]> row['breakpoint1']):
            bpt1_telo.append(0)
        if (list(select.chromStart)[1]< row['breakpoint1']) & (list(select.chromEnd)[1]> row['breakpoint1']):
            bpt1_telo.append(0)
        else:
            distance_1 = list(select.chromStart)[1] - int(row['breakpoint1'])
            distance_2 = int(row['breakpoint1']) - list(select.chromEnd)[0]
            total_distance = [distance_1, distance_2]
            abs_total=map(abs, total_distance)
            min_abs_distance = min(abs_total)
            bpt1_telo.append(min_abs_distance)

    for index, row in array.iterrows():
        select = telomere[telomere.chrom.astype(str) == str(row['chrom'])]
        if (list(select.chromStart)[0]< row['breakpoint2']) & (list(select.chromEnd)[0]> row['breakpoint2']):
            bpt2_telo.append(0)
        if (list(select.chromStart)[1]< row['breakpoint2']) & (list(select.chromEnd)[1]> row['breakpoint2']):
            bpt2_telo.append(0)
        else:
            distance_1 = list(select.chromStart)[1] - int(row['breakpoint2'])
            distance_2 = int(row['breakpoint2']) - list(select.chromEnd)[0]
            total_distance = [distance_1, distance_2]
            abs_total=map(abs, total_distance)
            min_abs_distance = min(abs_total)
            bpt2_telo.append(min_abs_distance)

    array['bpt1_to_telo']=bpt1_telo
    array['bpt2_to_telo']=bpt2_telo
    array.to_csv('{}_telomere.csv'.format(array_file),index=False)
#convert to two breakpoints format
input1=glob.glob('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/remove_50kb/*_telomere.csv')
appended_data = []
for infile in input1:
    print(str(infile))
    data = pd.read_csv(infile)
    appended_data.append(data)
appended_data = pd.concat(appended_data, axis=0)
appended_data1 = appended_data[appended_data.chrom=='X']
appended_data2 = appended_data[appended_data.chrom!='X']
appended_data2[["chrom"]]= appended_data2[["chrom"]].apply(pd.to_numeric)
new=appended_data1.append(appended_data2)
new.to_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/telomere_distance_distribution/all_telomere_distance.csv',index=False)

#wgs telemere distance
array_input=glob.glob('wgs_*_2breakpoints_format.csv')
telomere=pd.read_csv('hg19_telomere.txt',sep='\t')

for i in range(len(array_input)):
    array_file= array_input[i]
    print(str(array_file))
    array = pd.read_csv(array_file)
    bpt1_telo=[]
    bpt2_telo=[]
    for index, row in array.iterrows():
        select = telomere[telomere.chrom == str(row['chrom'])]
        if (list(select.chromStart)[0]< row['breakpoint1']) & (list(select.chromEnd)[0]> row['breakpoint1']):
            bpt1_telo.append(0)
        if (list(select.chromStart)[1]< row['breakpoint1']) & (list(select.chromEnd)[1]> row['breakpoint1']):
            bpt1_telo.append(0)
        else:
            distance_1 = list(select.chromStart)[1] - int(row['breakpoint1'])
            distance_2 = int(row['breakpoint1']) - list(select.chromEnd)[0]
            total_distance = [distance_1, distance_2]
            abs_total=map(abs, total_distance)
            min_abs_distance = min(abs_total)
            bpt1_telo.append(min_abs_distance)

    for index, row in array.iterrows():
        select = telomere[telomere.chrom.astype(str) == str(row['chrom'])]
        if (list(select.chromStart)[0]< row['breakpoint2']) & (list(select.chromEnd)[0]> row['breakpoint2']):
            bpt2_telo.append(0)
        if (list(select.chromStart)[1]< row['breakpoint2']) & (list(select.chromEnd)[1]> row['breakpoint2']):
            bpt2_telo.append(0)
        else:
            distance_1 = list(select.chromStart)[1] - int(row['breakpoint2'])
            distance_2 = int(row['breakpoint2']) - list(select.chromEnd)[0]
            total_distance = [distance_1, distance_2]
            abs_total=map(abs, total_distance)
            min_abs_distance = min(abs_total)
            bpt2_telo.append(min_abs_distance)

    array['bpt1_to_telo']=bpt1_telo
    array['bpt2_to_telo']=bpt2_telo
    array.to_csv('{}_telomere.csv'.format(array_file),index=False)

input1=glob.glob('wgs*2breakpoints_format.csv_telomere.csv')
appended_data = []
for infile in input1:
    print(str(infile))
    data = pd.read_csv(infile)
    appended_data.append(data)
appended_data = pd.concat(appended_data, axis=0)
appended_data1 = appended_data[(appended_data.chrom=='X') |(appended_data.chrom=='Y')]
appended_data2 = appended_data[(appended_data.chrom!='X')&(appended_data.chrom!='Y')]
appended_data2[["chrom"]]= appended_data2[["chrom"]].apply(pd.to_numeric)
new=appended_data1.append(appended_data2)
new.to_csv('~/Desktop/TCGA_RESULT/wgs_all_telomere_distance.csv',index=False)

#local array
array_input=glob.glob('diff*_2breakpoints_format.csv')
telomere=pd.read_csv('hg19_telomere.txt',sep='\t')

for i in range(len(array_input)):
    array_file= array_input[i]
    print(str(array_file))
    array = pd.read_csv(array_file)
    bpt1_telo=[]
    bpt2_telo=[]
    for index, row in array.iterrows():
        select = telomere[telomere.chrom == str(row['chrom'])]
        if (list(select.chromStart)[0]< row['breakpoint1']) & (list(select.chromEnd)[0]> row['breakpoint1']):
            bpt1_telo.append(0)
        if (list(select.chromStart)[1]< row['breakpoint1']) & (list(select.chromEnd)[1]> row['breakpoint1']):
            bpt1_telo.append(0)
        else:
            distance_1 = list(select.chromStart)[1] - int(row['breakpoint1'])
            distance_2 = int(row['breakpoint1']) - list(select.chromEnd)[0]
            total_distance = [distance_1, distance_2]
            abs_total=map(abs, total_distance)
            min_abs_distance = min(abs_total)
            bpt1_telo.append(min_abs_distance)

    for index, row in array.iterrows():
        select = telomere[telomere.chrom.astype(str) == str(row['chrom'])]
        if (list(select.chromStart)[0]< row['breakpoint2']) & (list(select.chromEnd)[0]> row['breakpoint2']):
            bpt2_telo.append(0)
        if (list(select.chromStart)[1]< row['breakpoint2']) & (list(select.chromEnd)[1]> row['breakpoint2']):
            bpt2_telo.append(0)
        else:
            distance_1 = list(select.chromStart)[1] - int(row['breakpoint2'])
            distance_2 = int(row['breakpoint2']) - list(select.chromEnd)[0]
            total_distance = [distance_1, distance_2]
            abs_total=map(abs, total_distance)
            min_abs_distance = min(abs_total)
            bpt2_telo.append(min_abs_distance)

    array['bpt1_to_telo']=bpt1_telo
    array['bpt2_to_telo']=bpt2_telo
    array.to_csv('{}_telomere.csv'.format(array_file),index=False)

input1=glob.glob('diff*2breakpoints_format.csv_telomere.csv')
appended_data = []
for infile in input1:
    print(str(infile))
    data = pd.read_csv(infile)
    appended_data.append(data)
appended_data = pd.concat(appended_data, axis=0)
appended_data1 = appended_data[appended_data.chrom=='X']
appended_data2 = appended_data[appended_data.chrom!='X']
appended_data2[["chrom"]]= appended_data2[["chrom"]].apply(pd.to_numeric)
new=appended_data1.append(appended_data2)
new.to_csv('~/Desktop/TCGA_RESULT/array_all_telomere_distance.csv',index=False)
