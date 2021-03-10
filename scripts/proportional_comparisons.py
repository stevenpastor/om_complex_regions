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
import numpy as np
import pandas as pd


# Define sample and alignment files:
with open("config", 'r') as c:
    reader = csv.reader(c)
    for line in reader:
        if len(line) > 0:
            if "SAMPLE" in line[0]:
                sample = re.sub('[SAMPLE=]', '', line[0])

filtered_align = 'results/'+sample+'_filter_molecules_self_alignment/'+sample+"_filtered_self_alignments"


def open_align_to_df(align):
    """
    place alignment file into pandas df
    convert scores col into floats
    """

    align_df = pd.read_csv(align, header=None, sep='\t')
    align_df[4] = align_df[4].astype(float)

    return align_df

def filter_positive_scores(align_df):
    """
    filter for only positive scores in filtered_align
    """

    align_df_positive = align_df[align_df[4] > 0]

    return align_df_positive

def find_unique_ids(align_df_positive):
    """
    # find all unique IDs from 1st and 3rd columns of align_df_positive
    """

    first_uniques = align_df_positive[0].unique()
    third_uniques = align_df_positive[2].unique()
    unique_ids = np.unique(np.concatenate((first_uniques, third_uniques), 0))

    return unique_ids

def find_mean_sd_and_one_sd_below(unique_ids, align_df_positive):
    """
    iterate through unique IDs and find in 1st or 3rd col of align_df_positive 
    and get mean, std dev, and mean-std dev
    """

    filtered_align_list_positive = []

    for i in unique_ids:
        series_one = align_df_positive[align_df_positive[0] == i][4]
        series_two = align_df_positive[align_df_positive[2] == i][4]
        combined_scores = series_one.append(series_two)
        average_score = combined_scores.mean()
        sample_sd = combined_scores.std() # watch for na
        one_sd_below_mean = average_score - sample_sd

        # obtain mapped combos above one sd below mean:
        first_above = align_df_positive.loc[(align_df_positive[0] == i) & (align_df_positive[4] > one_sd_below_mean)]
        second_above = align_df_positive.loc[(align_df_positive[2] == i) & (align_df_positive[4] > one_sd_below_mean)]
        combined_df = first_above.append(second_above)
        filtered_align_list_positive.append(combined_df.values.tolist())

    fix_me = []
    for i in filtered_align_list_positive:
        for j in i:
            fix_me.append(j)

    filtered_align_df_positive = pd.DataFrame(fix_me, columns = [0,1,2,3,4]) 
    filtered_align_df_positive = filtered_align_df_positive.drop_duplicates()

    return filtered_align_df_positive

def mapping_ids_dicts(unique_ids, filtered_align_df_positive):
    """
    makes dict from ids mapping to one another
    meeting criteria set above
    """

    ids_dict = {}
    for i in unique_ids:
        ids_dict[i] = []
        ids_dict[i].append(filtered_align_df_positive[filtered_align_df_positive[0] == i][2].tolist())
        ids_dict[i].append(filtered_align_df_positive[filtered_align_df_positive[2] == i][0].tolist())

    ids_dict_fixed = {}
    for i in ids_dict:
        for j in ids_dict[i]:
            if len(j) > 0 and i in ids_dict_fixed:
                for z in j:
                    ids_dict_fixed[i].append(z)
            elif len(j) > 0 and i not in ids_dict_fixed:
                ids_dict_fixed[i] = j

    return ids_dict_fixed

def combinations_comparison(ids_dict_fixed, align_df):
    """
    """

    import itertools

    combination_ids_proportions = []

    # list of tuples, where each tuple has: (id1, id2)
    result_list = list(itertools.combinations(ids_dict_fixed.keys(), 2))
    for ind,i in enumerate(result_list):
        combination_ids_proportions.append([])
        # proportion with regard to the total of both ids in tuples
        # i.e., shared_mapped_ids_between_id1_id2/unique(total_id1_mapped_ids + total_id2_mapped_ids):
        combination_ids_proportions[ind].append(i[0])
        combination_ids_proportions[ind].append(i[1])
        combination_ids_proportions[ind].append((len(list(set(ids_dict_fixed[i[0]]).intersection(ids_dict_fixed[i[1]]))) / len(list(set(ids_dict_fixed[i[0]] + ids_dict_fixed[i[1]]))))*100)

    combination_ids_proportions_annotated = []

    for ind,i in enumerate(combination_ids_proportions):
        combination_ids_proportions_annotated.append([])
        if len(align_df[align_df[0] == i[0]][1].tolist()) > 0 and len(align_df[align_df[0] == i[1]][1].tolist()) > 0:
            combination_ids_proportions_annotated[ind].append(i[0])
            combination_ids_proportions_annotated[ind].append(align_df[align_df[0] == i[0]][1].tolist()[0])
            combination_ids_proportions_annotated[ind].append(i[1])
            combination_ids_proportions_annotated[ind].append(align_df[align_df[0] == i[1]][1].tolist()[0])
            combination_ids_proportions_annotated[ind].append(i[2])
        elif len(align_df[align_df[0] == i[0]][1].tolist()) == 0 and len(align_df[align_df[0] == i[1]][1].tolist()) > 0:
            combination_ids_proportions_annotated[ind].append(i[0])
            combination_ids_proportions_annotated[ind].append(align_df[align_df[2] == i[0]][3].tolist()[0])
            combination_ids_proportions_annotated[ind].append(i[1])
            combination_ids_proportions_annotated[ind].append(align_df[align_df[0] == i[1]][1].tolist()[0])
            combination_ids_proportions_annotated[ind].append(i[2])
        elif len(align_df[align_df[0] == i[0]][1].tolist()) > 0 and len(align_df[align_df[0] == i[1]][1].tolist()) == 0:
            combination_ids_proportions_annotated[ind].append(i[0])
            combination_ids_proportions_annotated[ind].append(align_df[align_df[0] == i[0]][1].tolist()[0])
            combination_ids_proportions_annotated[ind].append(i[1])
            combination_ids_proportions_annotated[ind].append(align_df[align_df[2] == i[1]][3].tolist()[0])
            combination_ids_proportions_annotated[ind].append(i[2])
        elif len(align_df[align_df[0] == i[0]][1].tolist()) == 0 and len(align_df[align_df[0] == i[1]][1].tolist()) == 0:
            combination_ids_proportions_annotated[ind].append(i[0])
            combination_ids_proportions_annotated[ind].append(align_df[align_df[2] == i[0]][3].tolist()[0])
            combination_ids_proportions_annotated[ind].append(i[1])
            combination_ids_proportions_annotated[ind].append(align_df[align_df[2] == i[1]][3].tolist()[0])
            combination_ids_proportions_annotated[ind].append(i[2])

    return combination_ids_proportions_annotated

def main():
    if not os.path.exists('results/'+sample+'_proportional_comparisons'):
        os.makedirs('results/'+sample+'_proportional_comparisons')

    align_df = open_align_to_df(filtered_align)
    align_df_positive = filter_positive_scores(align_df)
    unique_ids = find_unique_ids(align_df_positive)

    # keeping self-align complex scores with anchor >=X (only for groups with anchors) 1 sd above positive complex score means:
    point_statistics = find_mean_sd_and_one_sd_below(unique_ids, align_df_positive)
    point_statistics.to_csv('results/'+sample+'_proportional_comparisons/'+sample+"_removed_one_sd_below.tsv", sep='\t', index=False, header=None)

#    # DEPRECATED?
#    unique_ids = find_unique_ids(point_statistics)
#    mapping_dicts = mapping_ids_dicts(unique_ids, point_statistics)
#
#    # proportions of shared mapped ids (includes those below 1 sd of mean positive complex scores):
#    combos_compared = combinations_comparison(mapping_dicts, align_df)
#    with open(sample+"_combination_proportions_self_alignments.tsv", "w", newline="") as f:
#        writer = csv.writer(f, delimiter='\t', lineterminator="\n")
#        writer.writerows(combos_compared)

main()


