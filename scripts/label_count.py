"""
Purpose: Run after initial_check.py
Count the number of labels 5' of the refstartpos and the number of labels 3' of the refendpos per molecule spanning the region.
Gives us an idea of how many molecules really span vs those which barely do.
Helpful if some regions where we may need to go quite a bit 5' and/or 3' of the seg dups to truly obtain unambiguous alignments.


Notes:
Still need to figure out how to count molecules: convert start/end position to molecule position. 


Last updated: JW 5/27/2021
"""

from numpy.lib.shape_base import column_stack
import pandas as pd
import csv
import re
import os
import numpy as np


with open("config", 'r') as c:
    reader = csv.reader(c)
    for line in reader:
        if len(line) > 0:
            if "SAMPLE" in line[0]:
                sample = (' '.join([str(char) for char in line])).split('=')[1]
            elif "COMPLEX_CHR" in line[0]:
                complex_chr = (' '.join([str(char) for char in line])).split('=')[1]
                if complex_chr == "X":
                    complex_chr = 23
                elif complex_chr == "Y":
                    complex_chr = 24

header_lines = []
results_dir = 'results/{}_initial_genome_check'.format(sample) ## Where the output files from initial_check.py are

## Save header of mapping files
def read_header_comments(file):
    header_lines.clear()
    for row in file:
        if row[0] == '#':
            header_lines.append(row)
            yield row.split('#')[1].strip()


## Get header for dataframe and returns header
def get_header_line(file):
    with open(file, 'r', newline='') as f:
        decommented_lines = [line for line in csv.reader(read_header_comments(f))]
        header = decommented_lines[-2]
        header = [word for line in header for word in line.split("\t")]
        header[0] = header[0][2:] ## Remove the "h " from the first column name

        return header


def count_labels(xmap, qmap, contig):
    header = get_header_line(xmap)
    xmap_df = pd.read_csv(xmap, sep='\t', comment='#', names=header, index_col=None)
    header = get_header_line(qmap)
    qmap_df = pd.read_csv(qmap, sep='\t', comment='#', names=header, index_col=None)

    grouped = qmap_df.groupby('CMapId')
    count_list = []
    for name, group in grouped: 
        fiveprime_count = 0
        threeprime_count = 0
        for index, row in xmap_df.iterrows(): 
            if name == row['QryContigID']:
                queryID = row['QryContigID']
                queryStart = row['QryStartPos']
                queryEnd = row['QryEndPos']
                orientation = row['Orientation']

                for index, row in group.iterrows():
                    if row['CMapId'] == queryID and orientation == '+':
                        if row['Position'] < queryStart:
                            fiveprime_count += 1
                        elif row['Position'] > queryEnd:
                            threeprime_count += 1
                    else: ## if orientation is negative
                        if row['Position'] > queryStart:
                            fiveprime_count += 1
                        elif row['Position'] < queryEnd:
                            threeprime_count += 1
        # count_df['Molecule_ID'] = name
        # count_df["5prime_Labels"] = fiveprime_count
        # count_df["3prime_Labels"] = threeprime_count
        count_list.append([name, fiveprime_count, threeprime_count])
    count_df = pd.DataFrame(count_list, columns=['Molecule_ID', '5prime_Labels', '3prime_Labels'])

    count_df.to_csv("{}/contig{}_moleculeLabelCounts.txt".format(results_dir, contig), sep='\t', index=None)


def main():
    xmap_file = '{}/{}_fullContigs.xmap'.format(results_dir, sample)
    header = get_header_line(xmap_file)
    xmap_df = pd.read_csv(xmap_file, sep='\t', comment='#', names=header, index_col=None)
    contigs_list = xmap_df['QryContigID'].unique().tolist()

    for contig in contigs_list:
        mol_xmap_file = '{}/{}_fullContig{}_molecules.xmap'.format(results_dir, sample, contig)
        mol_qmap_file = '{}/{}_fullContig{}_molecules_q.cmap'.format(results_dir, sample, contig)
        count_labels(mol_xmap_file, mol_qmap_file, contig)


if __name__ == '__main__':
    main()