#!/usr/bin/python
import pandas as pd
import concurrent.futures
import numpy as np

# Read the file using Pandas
otutab = pd.read_csv("otu_table_final.tsv", sep='\t')

# Convert 'new_est_reads' column to float if it's not already
otutab['new_est_reads'] = otutab['new_est_reads'].astype(float)

# Convert 'taxa' and 'taxid' columns to string
otutab['taxa'] = otutab['taxa'].astype(str)
otutab['taxid'] = otutab['taxid'].astype(str)

# Get unique samples
samples = otutab['ID'].unique()

# Create matrix_otutab DataFrame with appropriate index and columns
matrix_otutab = pd.DataFrame(0.0, index=otutab['taxid'].unique(), columns=samples)

# Set the interval for printing progress
print_interval = 1000

# Convert 'new_est_reads' column to NumPy array with explicit dtype
new_est_reads_cpu = otutab['new_est_reads'].values.astype(np.float64)

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
    for i, (taxid, sample, value) in enumerate(zip(otutab['taxid'], otutab['ID'], new_est_reads_cpu)):
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

# Merge the 'taxa' column from 'otutab' DataFrame to 'matrix_otutab_result' DataFrame based on 'taxid' column
result = pd.merge(matrix_otutab_result, otutab_taxa, on='taxid', how='left')

# Save matrix_otutab_result to a CSV file
result.to_csv("final_otu_table_matrix.csv", index=False)

