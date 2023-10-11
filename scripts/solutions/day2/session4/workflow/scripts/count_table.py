'''
Merge gene counts from all samples of an assembly into a single table.
'''


import os
import pandas as pd
import sys


# Constants
FIELDS = ['Geneid', 'Reads_quant']
STR_TO_REMOVE = '_genes_read_quantification.tsv'


# Functions
def import_clean(table):
    print(f'Importing and cleaning quantification data from <{table}>')
    reads = pd.read_csv(table, sep='\t', comment='#')
    reads.rename(columns={reads.columns[-1]: 'Reads_quant'}, inplace=True)
    print('Sorting <gene> table by Chromosome then Start position')
    # New columns are simpler and will be used to properly reorder the table
    print('\tCreating temporary columns')
    # We get the unique Chr ID with set
    reads['Chr_new'] = reads['Chr'].apply(lambda x: ''.join(set(x.split(';'))))
    # We select the start of the first exon
    reads['Start_new'] = reads['Start'].apply(lambda x: int(x.split(';')[0]))
    print('\tSorting table')
    reads.sort_values(['Chr_new', 'Start_new'], ascending=[True, True],
                      inplace=True)
    print('\tRemoving temporary columns')
    reads.drop(['Chr_new', 'Start_new'], axis='columns', inplace=True)
    final_table = reads[FIELDS].set_index('Geneid', drop=True)
    return final_table


# Main code execution
if __name__ == '__main__':

    with open(snakemake.log[0], 'w') as logfile:
    
        # Redirect everything from the script to Snakemake log
        sys.stderr = sys.stdout = logfile

        print('Getting data from snakemake')
        list_of_files = snakemake.input
        count_table = snakemake.output.table

        output_dir = os.path.dirname(count_table)
        os.makedirs(output_dir, exist_ok=True)

        print(f'Initializing global table with <{list_of_files[0]}>')
        total_table = import_clean(list_of_files[0])

        for file in list_of_files[1:]:
            print(f'\tAdding data from <{file}>')
            tmp_table = import_clean(file)
            total_table = pd.concat([total_table, tmp_table], axis=1)

        print('Renaming columns')
        column_titles = [os.path.basename(x).replace(STR_TO_REMOVE, '')
                         for x in list_of_files]
        total_table.columns = column_titles

        print(f'Saving final table in <{count_table}>')
        total_table.to_csv(count_table, sep='\t', header=True, index=True)
        print('Done')
