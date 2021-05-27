"""
Purpose: Check if a region has 2 contiguous CMAPs from the de novo assembly. 
Get the molecules per CMAP and see if they support the 2 haplotypes. 
If they do, then the genome is done and does not need the pipeline.

Notes:
Script assumes assemblies are unzipped and named sample_output inside data folder.
File paths may need to be changed depending on assembly version.
Needs more testing if the correct molecules are being extracted.
Does not extract from bnx file currently. 

Last updated: JW 5/27/2021
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
            elif "COMPRESSION_TYPE" in line[0]:
                compression_type = (' '.join([str(char) for char in line])).split('=')[1]
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

## Get paths within de novo assembly folder for different versions
if compression_type == "zip":
    molecule_dir = 'data/{}_output/contigs/exp_refineFinal1/alignmol/merge'.format(sample) ## Where the molecule xmap, qmap, rmap are (no BNX file in this assembly)
else: 
    molecule_dir = 'data/{}_output/molecules'.format(sample)


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


## Outputs dataframes to tab-delimited files
def output_file(df_to_output, output_filename):
    with open('header_file.txt', 'w') as filehandle:
        for listitem in header_lines:
            filehandle.write('%s' % listitem)

    df_to_output.to_csv('data_file.txt', sep='\t', header=None, index=None)
    cmd = 'cat header_file.txt data_file.txt > {} && rm header_file.txt && rm data_file.txt'.format(output_filename)
    os.system(cmd)


## Open xmap file and check for contigs that span the entire region
## Should we have a 5kb buffer for what is considered spanning all the way? 
def open_xmap_file(xmap_file, chrom, start, end):
    print("Extracting contigs from xmap")
    header = get_header_line(xmap_file)
    xmap_df = pd.read_csv(xmap_file, sep='\t', comment='#', names=header, index_col=None)
    xmap_df = xmap_df.query("RefContigID == @chrom")
    contig_df = xmap_df.query("RefStartPos <= @start & RefEndPos >= @end or RefStartPos >= @end & RefEndPos <= @start") 
    ## Check if they are split mapped
    if contig_df.shape[0] == 0:
        rows_loop = []
        grouped = xmap_df.groupby('QryContigID')
        for name, group in grouped: 
            for index, row in group.iterrows():
                ## greater than start and less than end
                if (row['RefStartPos'] >= start and row['RefStartPos'] <= end) or (row['RefEndPos'] >= start and row['RefEndPos'] <= end):
                    rows_loop.append(row)
        sub_df = pd.DataFrame(rows_loop)
       
        ## Check if there is a gap between maps
        ## Gap is currently set to 100kb
        gap = 100000
        grouped2 = sub_df.groupby("QryContigID")
        rows_loop2 = []
        for name, group in grouped:
            group = group.assign(shifted_start=group.RefStartPos.shift(-1)).fillna(0)
            group = group.assign(shifted_end=group.RefEndPos.shift(-1)).fillna(0)
            for index, row in group.iterrows():
                if row['shifted_start'] - row['RefEndPos'] <= gap:
                    rows_loop2.append(row)
        contig_df = pd.DataFrame(rows_loop)

    return contig_df


## Get contigs from cmap file
def open_cmap_file(cmap_file, IDsToExtract):
    print("Extracting contigs from cmap")
    header = get_header_line(cmap_file)
    cmap_df = pd.read_csv(cmap_file, sep='\t', comment='#', names=header, index_col=None)
    cmap_df = cmap_df.query("CMapId == @IDsToExtract")

    return cmap_df


def closest(lst, K):
    return lst[min(range(len(lst)), key=lambda i: abs(lst[i] - K))]


## Get coordinates of molecules relative to the contig
## This takes a long time, maybe there is a quicker way to do this
def getRelevantContigCoordinates(contig):
    print("Getting relevant contig coordinates to molecule coordinates")
    merged_xmap_file = '{}/{}_fullContigs.xmap'.format(output_dir, sample)
    header = get_header_line(merged_xmap_file)
    merged_xmap_df = pd.read_csv(merged_xmap_file, sep='\t', comment='#', names=header, index_col=None)

    merged_qmap_file = '{}/{}_fullContigs_q.cmap'.format(output_dir, sample)
    header = get_header_line(merged_qmap_file)
    merged_qmap_df = pd.read_csv(merged_qmap_file, sep='\t', comment='#', names=header, index_col=None)

    merged_rmap_file = 'data/{}_output/contigs/exp_refineFinal1_sv/merged_smaps/exp_refineFinal1_merged_r.cmap'.format(sample, contig)
    header = get_header_line(merged_rmap_file)
    merged_rmap_df = pd.read_csv(merged_rmap_file, sep='\t', comment='#', names=header, index_col=None)

    ## Get relevant positions on contig
    print("Getting relevant positions on contig")
    if merged_xmap_df.shape[0] == 1:
        alignment = merged_xmap_df.loc[0, 'Alignment']
    else: 
        alignment = ''.join(merged_xmap_df['Alignment'])
    sub_qmap =  merged_qmap_df.query("CMapId == @contig")
    querydict = pd.Series(sub_qmap.Position.values, index=sub_qmap.SiteID).to_dict()
    chrom_ref = merged_rmap_df.query("CMapId == @complex_chr")
    refdict = pd.Series(chrom_ref.SiteID.values, index=chrom_ref.Position).to_dict()
    keylist = list(refdict.keys())

    ## Finds the nicking sites on the contig closest to region start and end
    print("Finding the nicking sites on the contig closest to region start and end")
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
    print("Checkpoint 1")
    while str(refdictstart) not in reflist:
        refdictstart = refdictstart - 1
    print("Checkpoint 2")
    while str(refdictend) not in reflist:
        refdictend = refdictend + 1

    ## Looks at alignment string to see what matched the nick sites closest to the region and gives SiteID
    print("Looking at alignment string to see what matched the nick sites closest to the region and gives SiteID")
    contigqstart = alignment_df.loc[str(refdictstart), 'querysite'] 
    contigqend = alignment_df.loc[str(refdictend), 'querysite']
    contigStartPos = min([querydict[int(contigqstart)], querydict[int(contigqend)]]) # convert site ID to contig position
    contigEndPos = max([querydict[int(contigqstart)], querydict[int(contigqend)]])

    return contigStartPos, contigEndPos


## If contigs span the entire region, extract those molecules 
## Would have to download bnx and finish this block of code if extracting from bnx file
def extract_molecules(contigs_list):
    for contig in contigs_list: 
        map_filename = 'exp_refineFinal1_contig{}'.format(contig)
        contigStartPos, contigEndPos = getRelevantContigCoordinates(contigs_list)
        ## Copy over contig_rmap 
        rmap_file = '{}/{}_r.cmap'.format(molecule_dir, map_filename)
        cmd = 'cp {} {}'.format(rmap_file, output_dir)
        os.system(cmd)

        ## Extract from xmap
        print("Extracting molecules from molecule files for contig {}".format(contig))
        xmap_file = '{}/{}.xmap'.format(molecule_dir, map_filename)
        header = get_header_line(xmap_file)
        xmap_df = pd.read_csv(xmap_file, sep='\t', comment='#', names=header, index_col=None)
        molecule_df = xmap_df.query("RefStartPos < @contigEndPos and RefEndPos > @contigStartPos")
        molecules_list = list(molecule_df['QryContigID'])

        sub_xmap_df = xmap_df.query("QryContigID in @molecules_list")

        xmap_output = '{}/{}_fullContig{}_molecules.xmap'.format(output_dir, sample, contig)
        output_file(sub_xmap_df, xmap_output)

        ## Extract from qmap
        qmap_file = '{}/{}_q.cmap'.format(molecule_dir, map_filename)
        header = get_header_line(qmap_file)
        qmap_df = pd.read_csv(qmap_file, sep='\t', comment='#', names=header, index_col=None)
        sub_qmap_df = qmap_df.query("CMapId in @molecules_list")

        qmap_output = '{}/{}_fullContig{}_molecules_q.cmap'.format(output_dir, sample, contig)
        output_file(sub_qmap_df, qmap_output)

        rmap_file = '{}/{}_r.cmap'.format(molecule_dir, map_filename)   
        cmd = 'cp {} {}'.format(rmap_file, output_dir)
        os.system(cmd)

        ## Extract from BNX file instead
        # for molecule in molecules_list:
        #     bnxfile = '{}_output'.format(sample)


def main():
    if not os.path.exists('results/{}_initial_genome_check'.format(sample)):
        os.makedirs('results/{}_initial_genome_check'.format(sample))

    ## Files
    merged_xmap = 'data/{}_output/contigs/exp_refineFinal1_sv/merged_smaps/exp_refineFinal1_merged.xmap'.format(sample)
    output_fullcontigs_xmap = '{}/{}_fullContigs.xmap'.format(output_dir, sample)
    merged_cmap = 'data/{}_output/contigs/exp_refineFinal1_sv/merged_smaps/exp_refineFinal1_merged_q.cmap'.format(sample)
    output_fullcontigs_qcmap = '{}/{}_fullContigs_q.cmap'.format(output_dir, sample)

    fullcontigs_df = open_xmap_file(merged_xmap, complex_chr, complex_start, complex_end)
    ## We want to know which genomes have 2 (or more) full-length contigs. Save list of genomes. 
    ## Process rest of script if there are 2 full-length contigs. 
    if fullcontigs_df.shape[0] > 0: 
        if sample not in list_fullcontig_genomes:
            list_fullcontig_genomes.append(sample)         
        output_file(fullcontigs_df, output_fullcontigs_xmap)
        contigs_list = fullcontigs_df['QryContigID'].unique().tolist()
        cmap_df = open_cmap_file(merged_cmap, contigs_list)
        output_file(cmap_df, output_fullcontigs_qcmap)
        extract_molecules(contigs_list)
        print("{} full-length contigs for sample {}".format(fullcontigs_df.shape[0], sample))
    else: 
        print("No full-length contigs for sample {}".format(sample))

    with open('results/genomes_completeContigs_list.txt'.format(sample), 'a') as filehandle:
        for genomes in list_fullcontig_genomes:
            filehandle.write('%s\n' % genomes)


if __name__ == '__main__':
    main()