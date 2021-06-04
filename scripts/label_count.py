"""
Purpose: Run after initial_check.py
Count the number of labels 5' of the refstartpos and the number of labels 3' of the refendpos per molecule spanning the region.
Gives us an idea of how many molecules really span vs those which barely do.
Helpful if some regions where we may need to go quite a bit 5' and/or 3' of the seg dups to truly obtain unambiguous alignments.

Last updated: JW 6/4/2021
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
            elif "COMPLEX_START" in line[0]:
                complex_start = int(re.sub('[COMPLEX_START=]', '', line[0]))
            elif "COMPLEX_END" in line[0]:
                complex_end = int(re.sub('[COMPLEX_END=]', '', line[0]))


header_lines = []
results_dir = 'results/{}_initial_genome_check'.format(sample) ## Where the output files from initial_check.py are
xmap_file = '{}/{}_fullContigs.xmap'.format(results_dir, sample)

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


def closest(lst, K):
    return lst[min(range(len(lst)), key=lambda i: abs(lst[i] - K))]


## Get rows contig xmap that overlap with region start and end (5' and 3' and exclude nested)
def contig_xmap(contigs_list):
    header = get_header_line(xmap_file)
    xmap_df = pd.read_csv(xmap_file, sep='\t', comment='#', names=header, index_col=None)

    grouped = xmap_df.groupby('QryContigID')
    rows_list = []
    for name, group in grouped:
        for index, row in group.iterrows():
            if row['RefStartPos'] < complex_start:
                rows_list.append(row)
            elif row['RefEndPos'] > complex_end:
                rows_list.append(row)
    xmap_prime_df = pd.DataFrame(rows_list, columns=header)

    return (xmap_prime_df)


def count_labels(contig, mol_xmap, mol_qmap, mol_rmap, start, end, contig_orientation):
    print("Getting relevant molecule coordinates to contig coordinates for contig{}".format(contig))
    header = get_header_line(mol_xmap)
    mol_xmap_df = pd.read_csv(mol_xmap, sep='\t', comment='#', names=header, index_col=None)

    header = get_header_line(mol_qmap)
    mol_qmap_df = pd.read_csv(mol_qmap, sep='\t', comment='#', names=header, index_col=None)

    header = get_header_line(mol_rmap)
    merged_rmap_df = pd.read_csv(mol_rmap, sep='\t', comment='#', names=header, index_col=None)

    count_list = []
    for index, row in mol_xmap_df.iterrows():
        molecule_ID = row['QryContigID']
        orientation = row['Orientation']
        alignment = row['Alignment']
        sub_qmap =  mol_qmap_df.query("CMapId == @molecule_ID")
        querydict = pd.Series(sub_qmap.Position.values, index=sub_qmap.SiteID).to_dict()
        chrom_ref = merged_rmap_df.query("CMapId == @contig")
        refdict = pd.Series(chrom_ref.SiteID.values, index=chrom_ref.Position).to_dict()
        keylist = list(refdict.keys())

        ## Finds the nicking sites on the contig closest to region start and end
        refstart_nicksite = closest(keylist, start)
        refend_nicksite = closest(keylist, end) 
        alignmentlist = re.split('\(|\)', alignment)
        alignmentlist = list(filter(None, alignmentlist))
        alignment_df = pd.DataFrame(alignmentlist)
        alignment_df[['refsite', 'querysite']] = alignment_df[0].str.split(',', expand=True)
        alignment_df.drop(columns=0, inplace=True)
        alignment_df = alignment_df.set_index('refsite', drop=False)
        refdictstart = refdict[refstart_nicksite]
        refdictend = refdict[refend_nicksite]

        ## Looks at alignment string to see what matched the nick sites closest to the region and gives SiteID
        refsite_list = alignment_df['refsite'].to_list()
        if str(refdictstart) in refsite_list:
            molqstart = int(alignment_df.loc[str(refdictstart), 'querysite'])
        else:
            molqstart = 0
        if str(refdictend) in refsite_list:
            molqend = int(alignment_df.loc[str(refdictend), 'querysite'])
        else:
            molqend = 0
        
        site_list = alignment_df['querysite'].to_list()
        fiveprime_count = 0
        threeprime_count = 0
        if contig_orientation == '+' and orientation == '+':
            for site in site_list:
                if int(site) < molqstart:
                    fiveprime_count += 1
                elif int(site) > molqend:
                    threeprime_count += 1
        elif contig_orientation == '-' and orientation == '-':
            for site in site_list:
                if int(site) > molqstart:
                    fiveprime_count += 1
                elif int(site) < molqend:
                    threeprime_count += 1
        elif contig_orientation == '+' and orientation == '-':
            for site in site_list:
                if int(site) > molqstart:
                    threeprime_count += 1
                elif int(site) < molqend:
                    fiveprime_count += 1
        elif contig_orientation == '-' and orientation == '+':
            for site in site_list:
                if int(site) < molqstart:
                    threeprime_count += 1
                elif int(site) > molqend:
                    fiveprime_count += 1

        count_list.append([molecule_ID, fiveprime_count, threeprime_count, orientation, contig_orientation])

    count_df = pd.DataFrame(count_list, columns=['Molecule_ID', '5prime_Labels', '3prime_Labels', 'mol_orientation', 'contig_orientation'])
    count_df.to_csv("{}/contig{}_moleculeLabelCounts.txt".format(results_dir, contig), sep='\t', index=None)

    # cmd = 'rm {}/contig{}_startEndPos.txt'.format(results_dir, contig)
    # os.system(cmd)


def main():
    header = get_header_line(xmap_file)
    xmap_df = pd.read_csv(xmap_file, sep='\t', comment='#', names=header, index_col=None)
    contigs_list = xmap_df['QryContigID'].unique().tolist()

    xmap_prime_df = contig_xmap(contigs_list)

    for contig in contigs_list:
        mol_xmap_file = '{}/{}_fullContig{}_molecules.xmap'.format(results_dir, sample, contig)
        mol_qmap_file = '{}/{}_fullContig{}_molecules_q.cmap'.format(results_dir, sample, contig)
        mol_rmap_file = '{}/exp_refineFinal1_contig{}_r.cmap'.format(results_dir, contig)
        with open("{}/contig{}_startEndPos.txt".format(results_dir, contig), 'r') as c:
            reader = csv.reader(c)
            for line in reader:
                if len(line) > 0:
                    if "START" in line[0]:
                        start = int(float((' '.join([str(char) for char in line])).split('=')[1]))
                    elif "END" in line[0]:
                        end = int(float((' '.join([str(char) for char in line])).split('=')[1]))

        for index, row in xmap_prime_df.iterrows():
            if row['QryContigID'] == contig and row['RefStartPos'] < complex_start:
                contig_orientation = row['Orientation']
            if row['QryContigID'] == contig and row['RefEndPos'] > complex_end:
                contig_orientation = row['Orientation']

        count_labels(contig, mol_xmap_file, mol_qmap_file, mol_rmap_file, start, end, contig_orientation)


if __name__ == '__main__':
    main()