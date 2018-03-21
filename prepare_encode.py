import pandas as pd

if __name__ == '__main__':
    print('preparing encode df')
    #read encode csv
    encode_tss = pd.read_csv('encodehg19_tss_v2.csv')
    #remove decoy
    decoy = []
    for i in encode_tss['chrom']:
        if i.startswith('ch'):
            decoy.append(i)
    for i in decoy:
        encode_tss = encode_tss[encode_tss.chrom != i]
    #create TSS window
    encode_tss['TSS_back'] = encode_tss['TSS'] - 100000
    encode_tss['TSS_forth'] = encode_tss['TSS'] + 100000
    #construct a dataframe with 'aliquot_id' and 'ENSTid' as column name
    eid = list(encode_tss['ENST_id'])
    empty_df = pd.DataFrame(columns = eid)
    #rename duplicated column name with X_ENSTid and Y_ENSTid
    cols=pd.Series(empty_df.columns)
    for dup in empty_df.columns.get_duplicates():
        cols[empty_df.columns.get_loc(dup)]=[str(d_idx)+'_'+dup if d_idx!=0 else dup for d_idx in ['X', 'Y']]
    #replace original 'ENSTid' column with renamed ids
    new_id = list(cols)
    new_id_cols=pd.Series(new_id)
    encode_tss['ENSTid'] = new_id_cols
    encode_tss.to_csv('renamed_encode_tss.csv', index=False)
