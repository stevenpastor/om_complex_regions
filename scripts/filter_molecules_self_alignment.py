"""
# Goal:
# 1. Find matches where complex labels only sum >0 alignment scores
# 2. Of matches meeting criteria 1, further filter matches for those 
# with >=num_anchor_labels
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
                

self_align = "results/"+sample+"_merge_into_one_cmap_self_align/"+sample+"_merged_self_aligned.align"
idmap = "results/"+sample+"_merge_into_one_cmap_self_align/"+sample+"_merged.idmap"
five_anchor_file = "results/"+sample+"_extract_filter_anchor_labels/"+sample+"_five_query_labels_in_anchor"
three_anchor_file = "results/"+sample+"_extract_filter_anchor_labels/"+sample+"_three_query_labels_in_anchor"
spanning_anchor_file = "results/"+sample+"_extract_filter_anchor_labels/"+sample+"_spanning_query_labels_in_anchor"
output_file = 'results/'+sample+'_filter_molecules_self_alignment/'+sample+"_filtered_self_alignments"


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
    complex and anchor_files have first element of complex labels
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

def find_complex_score_above_threshold(five_anch_file_dict, three_anch_file_dict, spanning_anch_file_dict, self_aligned, mol_idmap, output_file, num_anchor_labels, complex_labels):
    """
    using self-align file, obtain list of indices 
    where complex labels' average scores are above 
    some threshold, defined by complex_score variable
    """

    complex_scores_to_write = []

    for ind,a in enumerate(self_aligned):
        anchor_labs = []
        complex_labs = []
        complex_scores_to_write.append([])

        # 5 vs 5
        if a[2] in list(five_anch_file_dict.keys()) and a[3] in list(five_anch_file_dict.keys()) and a[0] == ">0":
            # for each label in 3rd column ID:
            for label_ind,label in enumerate(self_aligned[ind+1]):
                # if this label is in the anchor list for this ID:
                if label != "" and int(label)+1 in five_anch_file_dict[a[2]]:
                    anchor_labs.append(float(self_aligned[ind+3][label_ind]))
                # if this label is NOT in the anchor list for this ID (i.e., complex label):
                if label != "" and int(label)+1 not in five_anch_file_dict[a[2]]:
                    complex_labs.append(float(self_aligned[ind+3][label_ind]))
            if len(anchor_labs) >= num_anchor_labels and len(complex_labs) > 0:
                complex_scores_to_write[ind].append(a[2])
                complex_scores_to_write[ind].append("five")
                complex_scores_to_write[ind].append(a[3])
                complex_scores_to_write[ind].append("five")
                complex_scores_to_write[ind].append(sum(complex_labs))

        # 3 vs 3
        elif a[2] in list(three_anch_file_dict.keys()) and a[3] in list(three_anch_file_dict.keys()) and a[0] == ">0":
            # for each label in 3rd column ID:
            for label_ind,label in enumerate(self_aligned[ind+1]):
                # if this label is in the anchor list for this ID:
                if label != "" and int(label)+1 in three_anch_file_dict[a[2]]:
                    anchor_labs.append(float(self_aligned[ind+3][label_ind]))
                # if this label is NOT in the anchor list for this ID (i.e., complex label):
                if label != "" and int(label)+1 not in three_anch_file_dict[a[2]]:
                    complex_labs.append(float(self_aligned[ind+3][label_ind]))
            if len(anchor_labs) >= num_anchor_labels and len(complex_labs) > 0:
                complex_scores_to_write[ind].append(a[2])
                complex_scores_to_write[ind].append("three")
                complex_scores_to_write[ind].append(a[3])
                complex_scores_to_write[ind].append("three")
                complex_scores_to_write[ind].append(sum(complex_labs))

        # 5 vs spanning or spanning vs 5
        elif a[2] in list(five_anch_file_dict.keys()) and a[3] in list(spanning_anch_file_dict.keys()) and a[0] == ">0":
            # for each label in 3rd column ID:
            for label_ind,label in enumerate(self_aligned[ind+1]):
                # if this label is in the anchor list for this ID:
                if label != "" and int(label)+1 in five_anch_file_dict[a[2]]:
                    anchor_labs.append(float(self_aligned[ind+3][label_ind]))
                # if this label is NOT in the anchor list for this ID (i.e., complex label):
                if label != "" and int(label)+1 not in five_anch_file_dict[a[2]]:
                    complex_labs.append(float(self_aligned[ind+3][label_ind]))
            if len(anchor_labs) >= num_anchor_labels and len(complex_labs) > 0:
                complex_scores_to_write[ind].append(a[2])
                complex_scores_to_write[ind].append("five")
                complex_scores_to_write[ind].append(a[3])
                complex_scores_to_write[ind].append("spanning")
                complex_scores_to_write[ind].append(sum(complex_labs))

        elif a[2] in list(spanning_anch_file_dict.keys()) and a[3] in list(five_anch_file_dict.keys()) and a[0] == ">0":
            # for each label in 3rd column ID:
            for label_ind,label in enumerate(self_aligned[ind+1]):
                # if this label is in the anchor list for this ID:
                if label != "" and int(label)+1 in spanning_anch_file_dict[a[2]]:
                    anchor_labs.append(float(self_aligned[ind+3][label_ind]))
                # if this label is NOT in the anchor list for this ID (i.e., complex label):
                if label != "" and int(label)+1 not in spanning_anch_file_dict[a[2]]:
                    complex_labs.append(float(self_aligned[ind+3][label_ind]))
            if len(anchor_labs) >= num_anchor_labels and len(complex_labs) > 0:
                complex_scores_to_write[ind].append(a[2])
                complex_scores_to_write[ind].append("spanning")
                complex_scores_to_write[ind].append(a[3])
                complex_scores_to_write[ind].append("five")
                complex_scores_to_write[ind].append(sum(complex_labs))

        # 3 vs spanning or spanning vs 3
        elif a[2] in list(three_anch_file_dict.keys()) and a[3] in list(spanning_anch_file_dict.keys()) and a[0] == ">0":
            # for each label in 3rd column ID:
            for label_ind,label in enumerate(self_aligned[ind+1]):
                # if this label is in the anchor list for this ID:
                if label != "" and int(label)+1 in three_anch_file_dict[a[2]]:
                    anchor_labs.append(float(self_aligned[ind+3][label_ind]))
                # if this label is NOT in the anchor list for this ID (i.e., complex label):
                if label != "" and int(label)+1 not in three_anch_file_dict[a[2]]:
                    complex_labs.append(float(self_aligned[ind+3][label_ind]))
            if len(anchor_labs) >= num_anchor_labels and len(complex_labs) > 0:
                complex_scores_to_write[ind].append(a[2])
                complex_scores_to_write[ind].append("three")
                complex_scores_to_write[ind].append(a[3])
                complex_scores_to_write[ind].append("spanning")
                complex_scores_to_write[ind].append(sum(complex_labs))

        elif a[2] in list(spanning_anch_file_dict.keys()) and a[3] in list(three_anch_file_dict.keys()) and a[0] == ">0":
            # for each label in 3rd column ID:
            for label_ind,label in enumerate(self_aligned[ind+1]):
                # if this label is in the anchor list for this ID:
                if label != "" and int(label)+1 in spanning_anch_file_dict[a[2]]:
                    anchor_labs.append(float(self_aligned[ind+3][label_ind]))
                # if this label is NOT in the anchor list for this ID (i.e., complex label):
                if label != "" and int(label)+1 not in spanning_anch_file_dict[a[2]]:
                    complex_labs.append(float(self_aligned[ind+3][label_ind]))
            if len(anchor_labs) >= num_anchor_labels and len(complex_labs) > 0:
                complex_scores_to_write[ind].append(a[2])
                complex_scores_to_write[ind].append("spanning")
                complex_scores_to_write[ind].append(a[3])
                complex_scores_to_write[ind].append("three")
                complex_scores_to_write[ind].append(sum(complex_labs))

        # spanning vs spanning
        elif a[2] in list(spanning_anch_file_dict.keys()) and a[3] in list(spanning_anch_file_dict.keys()) and a[0] == ">0":
            # for each label in 3rd column ID:
            for label_ind,label in enumerate(self_aligned[ind+1]):
                # if this label is in the anchor list for this ID:
                if label != "" and int(label)+1 in spanning_anch_file_dict[a[2]]:
                    anchor_labs.append(float(self_aligned[ind+3][label_ind]))
                # if this label is NOT in the anchor list for this ID (i.e., complex label):
                if label != "" and int(label)+1 not in spanning_anch_file_dict[a[2]]:
                    complex_labs.append(float(self_aligned[ind+3][label_ind]))
            if len(anchor_labs) >= num_anchor_labels and len(complex_labs) > 0:
                complex_scores_to_write[ind].append(a[2])
                complex_scores_to_write[ind].append("spanning")
                complex_scores_to_write[ind].append(a[3])
                complex_scores_to_write[ind].append("spanning")
                complex_scores_to_write[ind].append(sum(complex_labs))

        # 5 or 3 vs nested; nested vs spanning:
        elif a[0] == ">0":
            all_others = []
            for m in mol_idmap:
                if a[2] == m[0]:
                    all_others.append(a[2])
                    all_others.append(m[1])
                if a[3] == m[0]:
                    all_others.append(a[3])
                    all_others.append(m[1])
            all_others.append(float(a[4]))
            complex_scores_to_write[ind].append(all_others[0])
            complex_scores_to_write[ind].append(all_others[1])
            complex_scores_to_write[ind].append(all_others[2])
            complex_scores_to_write[ind].append(all_others[3])
            complex_scores_to_write[ind].append(all_others[4])

    # HERE: STEVEN: MAKE SURE THIS OUTPUTS LIST OF LISTS
    # if $3 [2] or $4 [3] is five or three and its opposing molecule is nested ($4 [3] or $3 [2], respectively), 
    # get the anchor label list for the 5 or 3 molecule ID and if none of those are in the alignment, keep the alignment 
    #complex_scores_to_write: 188218  three   993260  three   0.7699999999999999
    #self_aligned: self-align file raw, without header
    #five_anch_file_dict {key: [anchor_labels]}
    #three_anch_file_dict {key: [anchor_labels]}
    # for loop complex_scores_to_write and where [1] five and [3] nested or [1] nested and [3] five or [1] three and [3] nested or [1] nested and [3] three
    # get the five or three ID's anchor label list from five_anch_file_dict[key] or three_anch_file_dict[key]
    # if the labels in the ind+1 (if five/three in [2]) or ind+2 (if five/three in [3] - in self_aligned, not complex_scores_to_write!)
    # are NOT in the anchor (count and want <1), then keep the line
    # elif the [1] or [3] does not have five or three with a nested, keep the line and do not look up the anchor labels 
    # example: 188218  three   993260  three   0.7699999999999999: just keep the line as is
    # put into list of lists and replace complex_scores_to_write below with the name of the new variable from here.
    complex_scores_to_write_fixed = []

    for ind,a in enumerate(self_aligned):
        for cind,i in enumerate(complex_scores_to_write):
#            complex_scores_to_write_fixed.append([])
            tmp_l = []

            # 5 vs nested and self-align 3rd col is fi
            if len(i) > 0 and a[2] == i[0] and a[3] == i[2] and i[1] == "five" and i[3] == "nested" and a[0] == ">0":
                tmp_l = []
                for l in self_aligned[ind+1]:
                    if l != "" and int(l)+1 in five_anch_file_dict[i[0]]:
                        tmp_l.append(l)

            # 5 vs nested and self-align 4th col is first ID
            elif len(i) > 0 and a[2] == i[2] and a[3] == i[0] and i[1] == "five" and i[3] == "nested" and a[0] == ">0":
                tmp_l = []
                for l in self_aligned[ind+2]:
                    if l != "" and int(l)+1 in five_anch_file_dict[i[0]]:
                        tmp_l.append(l)

            # nested vs 5 and self-align 3rd col is first ID
            elif len(i) > 0 and a[2] == i[0] and a[3] == i[2] and i[1] == "nested" and i[3] == "five" and a[0] == ">0":
                tmp_l = []
                for l in self_aligned[ind+1]:
                    if l != "" and int(l)+1 in five_anch_file_dict[i[2]]:
                        tmp_l.append(l)

            # nested vs 5 and self-align 4th col is first ID
            elif len(i) > 0 and a[2] == i[2] and a[3] == i[0] and i[1] == "nested" and i[3] == "five" and a[0] == ">0":
                tmp_l = []
                for l in self_aligned[ind+2]:
                    if l != "" and int(l)+1 in five_anch_file_dict[i[2]]:
                        tmp_l.append(l)

            # 3 vs nested and self-align 3rd col is first ID
            elif len(i) > 0 and a[2] == i[0] and a[3] == i[2] and i[1] == "three" and i[3] == "nested" and a[0] == ">0":
                tmp_l = []
                for l in self_aligned[ind+1]:
                    if l != "" and int(l)+1 in three_anch_file_dict[i[0]]:
                        tmp_l.append(l)

            # 3 vs nested and self-align 4th col is first ID
            elif len(i) > 0 and a[2] == i[2] and a[3] == i[0] and i[1] == "three" and i[3] == "nested" and a[0] == ">0":
                tmp_l = []
                for l in self_aligned[ind+2]:
                    if l != "" and int(l)+1 in three_anch_file_dict[i[0]]:
                        tmp_l.append(l)

            # nested vs 3 and self-align 3rd col is first ID
            elif len(i) > 0 and a[2] == i[0] and a[3] == i[2] and i[1] == "nested" and i[3] == "three" and a[0] == ">0":
                tmp_l = []
                for l in self_aligned[ind+1]:
                    if l != "" and int(l)+1 in three_anch_file_dict[i[2]]:
                        tmp_l.append(l)

            # nested vs 3 and self-align 4th col is first ID
            elif len(i) > 0 and a[2] == i[2] and a[3] == i[0] and i[1] == "nested" and i[3] == "three" and a[0] == ">0":
                tmp_l = []
                for l in self_aligned[ind+2]:
                    if l != "" and int(l)+1 in three_anch_file_dict[i[2]]:
                        tmp_l.append(l)

            if tmp_l:
                # assumption: more than 2 labels in anchor mapped between nested and 5 or 3:
                if len(tmp_l) > 2:
#                    print(tmp_l, i)
                    complex_scores_to_write_fixed.append(i)

    first_set = set(map(tuple, complex_scores_to_write))
    secnd_set = set(map(tuple, complex_scores_to_write_fixed))
    complex_scores_to_write_final = list(first_set.symmetric_difference(secnd_set))

    # write out to file:
    with open(output_file, 'w') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator="\n")
        for i in complex_scores_to_write_final:
            if len(i) > 0:
                writer.writerow(i)

    return


def main():
    if not os.path.exists('results/'+sample+'_filter_molecules_self_alignment'):
        os.makedirs('results/'+sample+'_filter_molecules_self_alignment')

    five_anch_file = open_tsv_file(five_anchor_file)
    five_anch_file_dict = fix_label_lists(five_anch_file)
    
    three_anch_file = open_tsv_file(three_anchor_file)
    three_anch_file_dict = fix_label_lists(three_anch_file)

    spanning_anch_file = open_tsv_file(spanning_anchor_file)
    spanning_anch_file_dict = fix_label_lists(spanning_anch_file)

    self_aligned = open_tsv_file(self_align)

    mol_idmap = open_tsv_file(idmap)

    complex_above = find_complex_score_above_threshold(five_anch_file_dict, three_anch_file_dict, spanning_anch_file_dict, self_aligned, mol_idmap, output_file, num_anchor_labels, complex_labels)


main()

