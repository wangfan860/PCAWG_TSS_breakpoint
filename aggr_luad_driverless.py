"""
luad driverless sv aggregate
"""
import glob
import os
import pandas as pd
import argparse

def get_args():
    '''
    Loads the parser
    '''
    # Main parser
    parser = argparse.ArgumentParser(description="get closest distance to tss")
    # Args
    required = parser.add_argument_group("Required input parameters")
    # Metadata from input table
    required.add_argument('-g', '--gene_block', required=True, help='Gene block.')
    return parser.parse_args()

def Distance_to_tss(bedpe, chr, tss, gene_end):
    #read one bedpe file
    aliquot_id_data = pd.read_csv(bedpe)
    select1 = aliquot_id_data[aliquot_id_data.chrom1.astype(str) == str(chr)]
    select2 = aliquot_id_data[aliquot_id_data.chrom2.astype(str) == str(chr)]
    if tss + 100000 - gene_end < 0:
        left1 = select1.query('{} <= start1 <= {} & strand1=="-"'.format(tss-100000, tss+100000))
        distance_left1 = left1.start1 - int(tss)

        left2 = select2.query('{} <= start2 <= {} & strand2=="-"'.format(tss-100000, tss+100000))
        distance_left2 = left2.start2 - int(tss)

        total_distance = pd.concat([distance_left1, distance_left2]).drop_duplicates()
        min_abs_distance = total_distance.abs().min()
        return min_abs_distance

    elif tss + 100000 - gene_end >= 0:
        left1 = select1.query('{} <= start1 <= {} & strand1=="-"'.format(tss-100000, gene_end))
        distance_left1 = left1.start1 - int(tss)

        right1 = select1.query('{} < start1 <= {} & strand1=="+"'.format(gene_end, tss+100000))
        distance_right1 = right1.start1 - int(tss)

        left2 = select2.query('{} <= start2 <= {} & strand2=="-"'.format(tss-100000, gene_end))
        distance_left2 = left2.start2 - int(tss)

        right2 = select2.query('{} < start2 <= {} & strand2=="+"'.format(gene_end, tss+100000))
        distance_right2 = right2.start2 - int(tss)

        total_distance = pd.concat([distance_left1, distance_right1,distance_left2, distance_right2]).drop_duplicates()
        min_abs_distance = total_distance.abs().min()
        return min_abs_distance

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

if __name__ == '__main__':
    args = get_args()
    gb = int(args.gene_block)
    print('preparing')
    pcawg_input = glob.glob('/gpfs/data/lyang-lab/users/fan/LUAD_driverless/bedpe/*.bedpe')
    map1=pd.read_csv('/gpfs/data/lyang-lab/users/fan/LUAD_driverless/gencode_hg19.csv')
    total_gene_list=list(chunks(map1.gene_id.tolist(), 55))
    total_chrom_list=list(chunks(map1.chrom.tolist(), 55))
    total_tss_list=list(chunks(map1.TSS.tolist(), 55))
    total_gene_end_list=list(chunks(map1.gene_end.tolist(), 55))
    gene_list = total_gene_list[gb]
    chrom_list = total_chrom_list[gb]
    tss_list = total_tss_list[gb]
    gene_end_list = total_gene_end_list[gb]
    table_list = []
    for i in range(len(gene_list)):
        print(gene_list[i])
        table = pd.DataFrame(columns=['aliquot_id', gene_list[i]])
        for bedpe in pcawg_input:
            print('processing {}'.format(os.path.basename(bedpe).split('.')[0]))
            aliquot_id = os.path.basename(bedpe).split('-')[2]
            aliquot_dict = Distance_to_tss(bedpe, chrom_list[i], tss_list[i], gene_end_list[i])
            table = table.append({'aliquot_id':aliquot_id, gene_list[i]:aliquot_dict},ignore_index=True)
        table_list.append(table)

    df = table_list[0]
    for df_ in table_list[1:]:
        df = df.merge(df_, on='aliquot_id')
    df.to_csv('/gpfs/data/lyang-lab/users/fan/LUAD_driverless/Table_distance_to_tss.{}.csv'.format(gb), index=False)
