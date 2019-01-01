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

def Distance_to_tss(bedpe, chr, tss):
    aliquot_id_data = pd.read_csv(bedpe)
    select = aliquot_id_data[aliquot_id_data.chrom.astype(str) == str(chr)]
    distance_1 = select.breakpoint1 - int(tss)
    distance_2 = select.breakpoint2- int(tss)
    total_distance = pd.concat([distance_1, distance_2]).drop_duplicates()
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
    Input = glob.glob('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/remove_50kb/*2breakpoints_format_removed_50kb.csv')
    map1=pd.read_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/gencode_hg38_56354gene.csv')
    total_gene_list=list(chunks(map1.gene_id2.tolist(), 54))
    total_chrom_list=list(chunks(map1.chr.tolist(), 54))
    total_tss_list=list(chunks(map1.TSS.tolist(), 54))

    gene_list = total_gene_list[gb]
    chrom_list = total_chrom_list[gb]
    tss_list = total_tss_list[gb]

    table_list = []
    for i in range(len(gene_list)):
        print(gene_list[i])
        table = pd.DataFrame(columns=['barcode', gene_list[i]])
        for bedpe in Input:
            print('processing {}'.format(os.path.basename(bedpe).split('.')[0]))
            aliquot_id = os.path.basename(bedpe)
            aliquot_dict = Distance_to_tss(bedpe, chrom_list[i], tss_list[i])
            table = table.append({'barcode':aliquot_id, gene_list[i]:aliquot_dict},ignore_index=True)
        table_list.append(table)

    df = table_list[0]
    for df_ in table_list[1:]:
        df = df.merge(df_, on='barcode')
    df.to_csv('/gpfs/data/lyang-lab/users/fan/breakpoint_tcga/remove_50kb/TCGA_rm50kb_distance_to_tss.{}.csv'.format(gb), index=False)
