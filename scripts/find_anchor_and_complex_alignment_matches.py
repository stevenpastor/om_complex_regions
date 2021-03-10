"""
Filtering XMAP for molecules with >=complex_labels
"""

# standard library imports.
import csv
import re
import os


# Define sample and alignment files:
with open("config", 'r') as c:
    reader = csv.reader(c)
    for line in reader:
        if len(line) > 0:
            if "SAMPLE" in line[0]:
                sample = re.sub('[SAMPLE=]', '', line[0])
            elif "COMPLEX_LABELS" in line[0]:
                complex_labels = int(re.sub('[COMPLEX_LABELS=]', '', line[0]))
            elif "ANCHOR_LABELS" in line[0]:
                num_anchor_labels = int(re.sub('[ANCHOR_LABELS=]', '', line[0]))
                

five_anchor_file = "results/"+sample+"_extract_filter_anchor_labels/"+sample+"_five_query_labels_in_anchor"
five_anchored = "results/"+sample+"_separate_alignments/"+sample+"_five.xmap"
five_output_raw_file = "results/"+sample+"_extract_filter_anchor_labels/"+sample+"_five_query_labels_in_anchor_complex_filtered"

three_anchor_file = "results/"+sample+"_extract_filter_anchor_labels/"+sample+"_three_query_labels_in_anchor"
three_anchored = "results/"+sample+"_separate_alignments/"+sample+"_three.xmap"
three_output_raw_file = "results/"+sample+"_extract_filter_anchor_labels/"+sample+"_three_query_labels_in_anchor_complex_filtered"

spanning_anchor_file = "results/"+sample+"_extract_filter_anchor_labels/"+sample+"_spanning_query_labels_in_anchor"
spanning_anchored = "results/"+sample+"_separate_alignments/"+sample+"_completely_spanning.xmap"
spanning_output_raw_file = "results/"+sample+"_extract_filter_anchor_labels/"+sample+"_spanning_query_labels_in_anchor_complex_filtered"


def open_tsv_file(tsv_file):
    """
    open one TSV file as list of lists
    one line = one list
    """
    with open(tsv_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        tsv_file_lst_lsts = [line for line in reader if "#" not in line[0]]

    return tsv_file_lst_lsts

def fix_label_lists(tsv_file_lst_lsts):
    """
    complex files have first element of complex labels
    for an ID
    this list is imported as one string and needs to be
    list of integers, which is done in this function
    output is dict with ID as key and labels as value
    labels are ints
    """
    comp_dict = {}
    for ele in tsv_file_lst_lsts:
        comp_dict[ele[1]] = [int(i) for i in ele[0].split(',') ]
    
    return comp_dict

def find_complex_labels_above_threshold(anch_dict, group_anchored, complex_labels, output_raw_file):
    """
    In XMAP, find molecules with >=complex_labels mapping
    """

    complex_dict = {}
    for line in anch_dict:
        # go into XMAP (group_anchored) and find the ID matching anchor label file:
        for x in group_anchored:
            if line == x[1]:
                # following 2 will fix tuples in col14 and make into list
                # of only query labels: XMAP col14 format: (REF,QUERY)
                tmp = re.sub('[(]', '', x[13])
                tmp_two = re.sub('[)]', ',', tmp).split(',')[:-1][1::2]
                tmp_three = [int(i) for i in tmp_two]
                complex_labels_set = list(set(tmp_three) - set(anch_dict[line]))
                if x[1] not in complex_dict:
                    complex_dict[x[1]] = complex_labels_set
                else:
                    for n in complex_labels_set:
                        complex_dict[x[1]].append(n)

    myfile = open(output_raw_file, 'w')

    for i in complex_dict:
        if len(complex_dict[i]) >= complex_labels:
            myfile.write("%s\t%s\n" % (i, complex_dict[i]))

    myfile.close()

    return


def main():
    # 5':
    anch_file = open_tsv_file(five_anchor_file)
    anch_file_dict = fix_label_lists(anch_file)
    group_anchored = open_tsv_file(five_anchored)
    complex_above = find_complex_labels_above_threshold(anch_file_dict, group_anchored, complex_labels, five_output_raw_file)

    # 3':
    anch_file = open_tsv_file(three_anchor_file)
    anch_file_dict = fix_label_lists(anch_file)
    group_anchored = open_tsv_file(three_anchored)
    complex_above = find_complex_labels_above_threshold(anch_file_dict, group_anchored, complex_labels, three_output_raw_file)

    # Spanning:
    anch_file = open_tsv_file(spanning_anchor_file)
    anch_file_dict = fix_label_lists(anch_file)
    group_anchored = open_tsv_file(spanning_anchored)
    complex_above = find_complex_labels_above_threshold(anch_file_dict, group_anchored, complex_labels, spanning_output_raw_file)


main()

