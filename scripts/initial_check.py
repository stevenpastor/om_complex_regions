"""
Purpose: Check if a region has 2 contiguous CMAPs from the de novo assembly. 
Get the molecules per CMAP and see if they support the 2 haplotypes. 
If they do, then the genome is done and does not need the pipeline.

Notes:
Script assumes assemblies are unzipped and named sample_output.
File paths in def main() may need to be changed.
Needs more testing.
Currently testing on assembly with no BNX file, so extracting molecules from
output/contigs/exp_refineFinal1/alignmol/merge/.

Last updated: JW 5/25/2021
"""

import pandas as pd
import csv
import re
import os
import numpy as np

## Get sample/chr/start/end from config file:
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
list_fullcontig_genomes = []
output_dir = 'results/{}_initial_genome_check'.format(sample) ## Where the output files will appear
molecule_dir = '{}_output/contigs/exp_refineFinal1/alignmol/merge'.format(sample) ## Where the molecule xmap, qmap, rmap are (no BNX file in this assembly)

## Save header of mapping files
def read_header_comments(file):
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


## Outputs dataframes to tab-delimited files
def output_file(df_to_output, output_filename):
    with open('header_file.txt', 'w') as filehandle:
        for listitem in header_lines:
            filehandle.write('%s' % listitem)

    df_to_output.to_csv('data_file.txt', sep='\t', header=None, index=None)
    cmd = 'cat header_file.txt data_file.txt > {} && rm header_file.txt && rm data_file.txt'.format(output_filename)
    os.system(cmd)
    header_lines.clear()


## Open xmap file and check for contigs that span the entire region
## Should we have a 5kb buffer for what is considered spanning all the way? 
def open_xmap_file(xmap_file, chrom, start, end, output):
    header = get_header_line(xmap_file)
    xmap_df = pd.read_csv(xmap_file, sep='\t', comment='#', names=header, index_col=None)
    xmap_df = xmap_df.query("RefContigID == @chrom")
    sub_df = xmap_df.query("RefStartPos <= @start and RefEndPos >= @end")

    ## We want to know which genomes have 2 (or more) full-length contigs. Save list of genomes. 
    if sub_df.shape[0] >= 2: 
        if sample not in list_fullcontig_genomes:
            list_fullcontig_genomes.append(sample)         

        output_file(sub_df, output)

        return sub_df


## Get contigs from cmap file
def open_cmap_file(cmap_file, IDsToExtract, output):
    header = get_header_line(cmap_file)
    cmap_df = pd.read_csv(cmap_file, sep='\t', comment='#', names=header, index_col=None)
    cmap_df = cmap_df.query(" CMapId == @IDsToExtract")

    output_file(cmap_df, output)


def closest(lst, K):
    return lst[min(range(len(lst)), key=lambda i: abs(lst[i] - K))]


## Get coordinates of molecules relative to the contig
## This takes a long time, maybe there is a quicker way to do this
def getRelevantContigCoordinates(contig):
    merged_xmap_file = 'results/{}_initial_genome_check/{}_fullContigs.xmap'.format(sample, sample)
    header = get_header_line(merged_xmap_file)
    merged_xmap_df = pd.read_csv(merged_xmap_file, sep='\t', comment='#', names=header, index_col=None)

    merged_qmap_file = 'results/{}_initial_genome_check/{}_fullContigs_q.cmap'.format(sample, sample)
    header = get_header_line(merged_qmap_file)
    merged_qmap_df = pd.read_csv(merged_qmap_file, sep='\t', comment='#', names=header, index_col=None)

    merged_rmap_file = '{}_output/contigs/exp_refineFinal1_sv/merged_smaps/exp_refineFinal1_merged_r.cmap'.format(sample, contig)
    header = get_header_line(merged_rmap_file)
    merged_rmap_df = pd.read_csv(merged_rmap_file, sep='\t', comment='#', names=header, index_col=None)
    header_lines.clear()

    ## Get relevant positions on contig
    alignment = merged_xmap_df.loc[0, 'Alignment']
    sub_qmap =  merged_qmap_df.query("CMapId == @contig")
    querydict = pd.Series(sub_qmap.Position.values, index=sub_qmap.SiteID).to_dict()
    chrom_ref = merged_rmap_df.query("CMapId == @complex_chr")
    refdict = pd.Series(chrom_ref.SiteID.values, index=chrom_ref.Position).to_dict()
    keylist = list(refdict.keys())

    ## Finds the nickings sites on the contig closest to region start and end
    refstart_nicksite = closest(keylist, complex_start) 
    refend_nicksite = closest(keylist, complex_end) 
    alignmentlist = re.split('\(|\)', alignment)
    alignmentlist = list(filter(None, alignmentlist))
    alignment_df = pd.DataFrame(alignmentlist)
    alignment_df[['refsite', 'querysite']] = alignment_df[0].str.split(',', expand=True)
    alignment_df.drop(columns=0, inplace=True)
    alignment_df = alignment_df.set_index('refsite', drop=False)
    reflist = alignment_df.index.values.tolist()
    refdictstart = refdict[refstart_nicksite]
    refdictend = refdict[refend_nicksite]
    while str(refdictstart) not in reflist:
        refdictstart = refdictstart - 1
    while str(refdictend) not in reflist:
        refdictend = refdictend + 1

    ## Looks at alignment string to see what matched the nick sites closest to the region and gives SiteID
    contigqstart = alignment_df.loc[str(refdictstart), 'querysite'] 
    contigqend = alignment_df.loc[str(refdictend), 'querysite']
    contigStartPos = min([querydict[int(contigqstart)], querydict[int(contigqend)]]) # convert site ID to contig position
    contigEndPos = max([querydict[int(contigqstart)], querydict[int(contigqend)]])

    return contigStartPos, contigEndPos


## If contigs span the entire region, extract those molecules 
## Extracting from alignmol map files because newer assemblies do not have bnx in the folder
## Would have to download bnx and finish this block of code if extracting from bnx file
def extract_molecules(contigs_list):
        for contig in contigs_list: 
            contigStartPos, contigEndPos = getRelevantContigCoordinates(contigs_list)

            ## Extract from xmap
            xmap_file = '{}/exp_refineFinal1_contig{}.xmap'.format(molecule_dir, contig)
            header = get_header_line(xmap_file)
            xmap_df = pd.read_csv(xmap_file, sep='\t', comment='#', names=header, index_col=0)
            molecule_df = xmap_df.query("RefStartPos < @contigEndPos and RefEndPos > @contigStartPos")
            molecules_list = list(molecule_df['QryContigID'])

            sub_xmap_df = xmap_df.query("QryContigID in @molecules_list")

            xmap_output = '{}/{}_fullContig{}_molecules.xmap'.format(output_dir, sample, contig)
            output_file(sub_xmap_df, xmap_output)

            ## Extract from qmap
            qmap_file = '{}/exp_refineFinal1_contig{}_q.cmap'.format(molecule_dir, contig)
            header = get_header_line(qmap_file)
            qmap_df = pd.read_csv(qmap_file, sep='\t', comment='#', names=header, index_col=0)
            sub_qmap_df = qmap_df.query("CMapId in @molecules_list")

            qmap_output = '{}/{}_fullContig{}_molecules_q.cmap'.format(output_dir, sample, contig)
            output_file(sub_qmap_df, qmap_output)

            ## Extract from BNX file instead
            # for molecule in molecules_list:
            #     bnxfile = '{}_output'.format(sample)


def main():
    if not os.path.exists('results/{}_initial_genome_check'.format(sample)):
        os.makedirs('results/{}_initial_genome_check'.format(sample))

    ## Files
    merged_xmap = '{}_output/contigs/exp_refineFinal1_sv/merged_smaps/exp_refineFinal1_merged.xmap'.format(sample)
    output_fullcontigs_xmap = '{}/{}_fullContigs.xmap'.format(output_dir, sample)
    merged_cmap = '{}_output/contigs/exp_refineFinal1_sv/merged_smaps/exp_refineFinal1_merged_q.cmap'.format(sample)
    output_fullcontigs_qcmap = '{}/{}_fullContigs_q.cmap'.format(output_dir, sample)


    fullcontigs_df = open_xmap_file(merged_xmap, complex_chr, complex_start, complex_end, output_fullcontigs_xmap)
    contigs_list = fullcontigs_df['QryContigID'].tolist()
    if len(contigs_list) >= 2:
        open_cmap_file(merged_cmap, contigs_list, output_fullcontigs_qcmap)
        extract_molecules(contigs_list)
    else: 
        print("Only {} full-length contigs for sample {}".format(len(contigs_list), sample))


    with open('results/{}_initial_genome_check/genomes_completeContigs_list.txt'.format(sample), 'w') as filehandle:
        for genomes in list_fullcontig_genomes:
            filehandle.write('%s\n' % genomes)


if __name__ == '__main__':
    main()