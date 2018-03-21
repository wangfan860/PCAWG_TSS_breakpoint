"""
pcawg sv aggregate
"""

import glob
import os
import pandas as pd

def count_gene(bedpe, encode_tss, empty_dict):
    #read one bedpe file
    aliquot_id = os.path.basename(bedpe).split('.')[0]
    aliquot_id_data = pd.read_csv(bedpe)
    aliquot_id_dict = empty_dict
    #add aliquot_id to dict
    aliquot_id_dict['aliquot_id'] = aliquot_id
    #filter encode_tss df
    start1 = list(set(aliquot_id_data['start1']))
    start2 = list(set(aliquot_id_data['start2']))
    all_start = start1 + list(set(start2) - set(start1))
    possible_frame = [encode_tss.query('TSS_back <= {} <= TSS_forth'.format(i)) for i in all_start]
    if possible_frame:
        possible_df = pd.concat(possible_frame).drop_duplicates()
        strand1 = []
        for i in range(aliquot_id_data.shape[0]):
            if aliquot_id_data['strand1'].iloc[i] == '-':
                strand1.append(possible_df.query('chr == "{}" & TSS_back <= {} <= TSS_forth'.format(str(aliquot_id_data['chrom1'].iloc[i]), int(aliquot_id_data['start1'].iloc[i]))))
        if strand1:
            strand1_pd = pd.concat(strand1)
            for i in strand1_pd['ENSTid']:
                aliquot_id_dict[i] += 1
        strand2 = []
        for i in range(aliquot_id_data.shape[0]):
            if aliquot_id_data['strand2'].iloc[i] == '-':
                strand2.append(possible_df.query('chr == "{}" & TSS_back <= {} <= TSS_forth'.format(str(aliquot_id_data['chrom2'].iloc[i]), int(aliquot_id_data['start2'].iloc[i]))))
        if strand2:
            strand2_pd = pd.concat(strand2)
            for i in strand2_pd['ENSTid']:
                aliquot_id_dict[i] += 1
        dup = []
        for i in range(aliquot_id_data.shape[0]):
            if aliquot_id_data['strand1'].iloc[i] == '-' and aliquot_id_data['strand2'].iloc[i] == '-':
                dup.append(possible_df.query('chr == "{}" & TSS_back <= {} <= TSS_forth & chr == "{}" & TSS_back <= {} <= TSS_forth'.format(str(aliquot_id_data['chrom1'].iloc[i]), int(aliquot_id_data['start1'].iloc[i]), str(aliquot_id_data['chrom2'].iloc[i]), int(aliquot_id_data['start2'].iloc[i]))))
        if dup:
            dup_pd = pd.concat(dup)
            for i in dup_pd['ENSTid']:
                aliquot_id_dict[i] -= 1
    return aliquot_id_dict

if __name__ == '__main__':
    print('preparing')
    #read encode csv
    encode_tss = pd.read_csv('renamed_encode_tss.csv')
    #construct a dataframe with 'aliquot_id' and 'ENSTid' as column name
    eid = list(encode_tss['ENSTid'])
    eid.insert(0, 'aliquot_id')
    empty_df = pd.DataFrame(columns = eid)
    #get a empty dict upon the empty dataframe which will store counts
    empty_dict = empty_df.to_dict()
    #compute
    pcawg_input = glob.glob('part_all/*.part.bedpe')
    for i in pcawg_input:
        #with an empty template output dict, encode_tss input df
        print('processing {}'.format(os.path.basename(i).split('.')[0]))
        #set default value
        for key in empty_dict.keys():
            empty_dict[key] = 0
        aliquot_dict = count_gene(i, encode_tss=encode_tss, empty_dict=empty_dict)
        empty_df = empty_df.append(aliquot_dict, ignore_index=True)
    #write to output
    empty_df.to_csv('2784_aggregated_pcawg_encode.output')


    # Wrong, inefficient concept


    # aliquot_id_data = pd.read_csv(args.input)
    # for i in range(encode_tss.shape[0]):
    #     for j in range(aliquot_id_data.shape[0]):
    #         if str(encode_tss['chr'].iloc[i]) == str(aliquot_id_data['chrom1'].iloc[j]) and aliquot_id_data['strand1'].iloc[j] == '-' and int(encode_tss['TSS_back'].iloc[i]) <= int(aliquot_id_data['start1'].iloc[j]) <= int(encode_tss['TSS_forth'].iloc[i]):
    #             empty_dict[encode_tss['ENSTid'].iloc[i]] += 1
    #             continue
    #         elif str(encode_tss['chr'].iloc[i]) == str(aliquot_id_data['chrom2'].iloc[j]) and aliquot_id_data['strand2'].iloc[j] == '-' and int(encode_tss['TSS_back'].iloc[i]) <= int(aliquot_id_data['start2'].iloc[j]) <= int(encode_tss['TSS_forth'].iloc[i]):
    #             empty_dict[encode_tss['ENSTid'].iloc[i]] += 1
    # #set index for 2nd round
    # secround = [x for x in range(aliquot_id_data.shape[0]) if x not in list(set(comp))]
    # for i in range(encode_tss.shape[0]):
    #     for j in secround:
    #         if aliquot_id_data['strand2'].iloc[j] == '+':
    #             continue
    #         if str(encode_tss['chr'].iloc[i]) == str(aliquot_id_data['chrom2'].iloc[j]) and int(encode_tss['TSS_back'].iloc[i]) <= int(aliquot_id_data['start2'].iloc[j]) <= int(encode_tss['TSS_forth'].iloc[i]):
    #             empty_dict[encode_tss['ENSTid'].iloc[i]] += 1
    # aliquot_id_df = empty_df.append(empty_dict, ignore_index=True)
    # aliquot_id_df.to_csv('{}.output'.format(aliquot_id))
