import pandas as pd
import glob
files_1=glob.glob('*cut_noheader.txt.1')
files_2=glob.glob('*cut_noheader.txt.2')
files_1.sort()
files_2.sort()
for i in range(len(files_1)):
    print i
    file1=pd.read_csv(files_1[i], header=None, sep='\t')
    file2=pd.read_csv(files_2[i], header=None, sep='\t')
    file1.columns=['chr','start','end','gene']
    file2.columns=['chr','start','end','probe','segment_mean']

    file1_s=file1.sort_values(['chr', 'start'], ascending=[True, True])
    file1_ss=file1_s.reset_index()
    file1_sss=file1_ss.drop(['index'], axis=1)

    file2_s=file2.sort_values(['chr', 'start'], ascending=[True, True])
    file2_ss=file2_s.reset_index()
    file2_sss=file2_ss.drop(['index'], axis=1)
    merge=pd.concat([file1_sss,file2_sss], axis=1)
    merge[['gene','segment_mean']].to_csv(files_1[i]+'_merge.csv')
 #for LUAD_driverless
import pandas as pd
import glob
files_1=glob.glob('*noheader.1')
files_2=glob.glob('*noheader.2')
files_1.sort()
files_2.sort()
for i in range(len(files_1)):
    print i
    file1=pd.read_csv(files_1[i], header=None, sep='\t')
    file2=pd.read_csv(files_2[i], header=None, sep='\t')
    file1.columns=['chr','start','end','gene']
    file2.columns=['chr','start','end','width','strand','id','mark','segment_mean']

    file1_s=file1.sort_values(['chr', 'start'], ascending=[True, True])
    file1_ss=file1_s.reset_index()
    file1_sss=file1_ss.drop(['index'], axis=1)

    file2_s=file2.sort_values(['chr', 'start'], ascending=[True, True])
    file2_ss=file2_s.reset_index()
    file2_sss=file2_ss.drop(['index'], axis=1)
    merge=pd.concat([file1_sss,file2_sss], axis=1)
    merge[['gene','segment_mean']].to_csv(files_1[i]+'_merge.csv')
