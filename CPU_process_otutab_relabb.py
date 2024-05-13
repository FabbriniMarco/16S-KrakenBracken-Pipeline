#!/usr/bin/python
import os
import pandas as pd
import concurrent.futures
import numpy as np

# Read the file using Pandas
otutab = pd.read_csv("otu_table_final.tsv", sep='\t')

# Convert 'fraction_total_reads' column to float if it's not already
otutab['fraction_total_reads'] = otutab['fraction_total_reads'].astype(float)

# Convert 'taxa' and 'taxid' columns to string
otutab['taxa'] = otutab['taxa'].astype(str)
otutab['taxid'] = otutab['taxid'].astype(str)

# Get unique samples
samples = otutab['ID'].unique()

# Create matrix_otutab DataFrame with appropriate index and columns
matrix_otutab = pd.DataFrame(0.0, index=otutab['taxid'].unique(), columns=samples)

# Set the interval for printing progress
print_interval = 1000

# Convert 'fraction_total_reads' column to NumPy array with explicit dtype
fraction_total_reads_cpu = otutab['fraction_total_reads'].values.astype(np.float64)

# Create a dictionary to map taxa to their corresponding index
taxid_index_map = {tax: i for i, tax in enumerate(matrix_otutab.index)}

# Use a list for column names
matrix_otutab_columns_list = list(matrix_otutab.columns)

# Create a 2D NumPy array for matrix_otutab on CPU
matrix_otutab_cpu = matrix_otutab.values.astype(np.float64)

def process_row(taxid, sample, value):
    if taxid in taxid_index_map:
        row_index = taxid_index_map[taxid]
        col_index = matrix_otutab_columns_list.index(sample)
        matrix_otutab_cpu[row_index, col_index] = value

# Create a ThreadPoolExecutor
with concurrent.futures.ThreadPoolExecutor() as executor:
    # Iterate through rows and columns using ThreadPoolExecutor
    futures = []
    for i, (taxid, sample, value) in enumerate(zip(otutab['taxid'], otutab['ID'], fraction_total_reads_cpu)):
        if i % print_interval == 0:
            print(f"{i}/{len(otutab)}")
        
        futures.append(executor.submit(process_row, taxid, sample, value))
    # Wait for all threads to finish
    concurrent.futures.wait(futures)

# Convert the result back to a Pandas DataFrame
matrix_otutab_result = pd.DataFrame(matrix_otutab_cpu, index=matrix_otutab.index, columns=matrix_otutab.columns)

matrix_otutab_result.index.name = 'taxid'
otutab_taxa = otutab[['taxid', 'taxa']]
otutab_taxa = otutab_taxa.drop_duplicates()

result = pd.merge(matrix_otutab_result, otutab_taxa, on='taxid', how='left')

# Save result to a CSV file
result.to_csv("final_otu_table_matrix_relabb.csv", index=False)

if not os.path.exists('taxa_summary'):
    os.makedirs('taxa_summary')

grouped_otutab = otutab[['taxonomy_lvl', 'taxa']].groupby('taxonomy_lvl')
taxlevel_to_path = {'P': 'taxa_summary/otu_table_L2',
                    'F': 'taxa_summary/otu_table_L5',
                    'G': 'taxa_summary/otu_table_L6'}

for lvl, group in grouped_otutab:
    subset_result = result[result['taxa'].isin(group['taxa'])]
    subset_result = subset_result.iloc[: , 1:]
    subset_result.set_index(subset_result.columns[-1], inplace=True)
    filename = taxlevel_to_path[lvl] + '.tsv'
    subset_result.to_csv(filename,sep='\t', index=True)

