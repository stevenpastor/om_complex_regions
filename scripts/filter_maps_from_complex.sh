#!/bin/sh


#################
CONFIG="config"

SAMPLE=$(grep -w "^SAMPLE" $CONFIG | sed 's/SAMPLE=//g')
FIVE_FILTER="results/"$SAMPLE"_extract_filter_anchor_labels/"$SAMPLE"_five_query_labels_in_anchor_complex_filtered"
THREE_FILTER="results/"$SAMPLE"_extract_filter_anchor_labels/"$SAMPLE"_three_query_labels_in_anchor_complex_filtered"
SPANNING_FILTER="results/"$SAMPLE"_extract_filter_anchor_labels/"$SAMPLE"_spanning_query_labels_in_anchor_complex_filtered"
FIVE_XMAP="results/"$SAMPLE"_separate_alignments/"$SAMPLE"_five.xmap"
THREE_XMAP="results/"$SAMPLE"_separate_alignments/"$SAMPLE"_three.xmap"
SPANNING_XMAP="results/"$SAMPLE"_separate_alignments/"$SAMPLE"_completely_spanning.xmap"
FIVE_QCMAP="results/"$SAMPLE"_separate_alignments/"$SAMPLE"_five_q.cmap"
THREE_QCMAP="results/"$SAMPLE"_separate_alignments/"$SAMPLE"_three_q.cmap"
SPANNING_QCMAP="results/"$SAMPLE"_separate_alignments/"$SAMPLE"_completely_spanning_q.cmap"
#################

mkdir -p "results/"$SAMPLE"_filtered_maps"

while read line1 line2
do
    grep -v "#" $FIVE_XMAP | awk -v id="$line1" '$2 == id'
done <$FIVE_FILTER > "results/"$SAMPLE"_filtered_maps/"$SAMPLE"_five.xmap"

while read line1 line2
do
    grep -v "#" $THREE_XMAP | awk -v id="$line1" '$2 == id'
done <$THREE_FILTER > "results/"$SAMPLE"_filtered_maps/"$SAMPLE"_three.xmap"

while read line1 line2
do
    grep -v "#" $SPANNING_XMAP | awk -v id="$line1" '$2 == id'
done <$SPANNING_FILTER > "results/"$SAMPLE"_filtered_maps/"$SAMPLE"_completely_spanning.xmap"

while read line1 line2
do
    grep -v "#" $FIVE_QCMAP | awk -v id="$line1" '$1 == id'
done <$FIVE_FILTER > "results/"$SAMPLE"_filtered_maps/"$SAMPLE"_five_q.cmap"

while read line1 line2
do
    grep -v "#" $THREE_QCMAP | awk -v id="$line1" '$1 == id'
done <$THREE_FILTER > "results/"$SAMPLE"_filtered_maps/"$SAMPLE"_three_q.cmap"

while read line1 line2
do
    grep -v "#" $SPANNING_QCMAP | awk -v id="$line1" '$1 == id'
done <$SPANNING_FILTER > "results/"$SAMPLE"_filtered_maps/"$SAMPLE"_completely_spanning_q.cmap"

