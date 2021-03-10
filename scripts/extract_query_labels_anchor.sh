#!/bin/sh


# Define sample and alignment files:
CONFIG="config"

SAMPLE=$(grep -w "^SAMPLE" $CONFIG | sed 's/SAMPLE=//g')
COMPLEX_CHR=$(grep -w "^COMPLEX_CHR" $CONFIG | sed 's/COMPLEX_CHR=//g')
COMPLEX_START=$(grep -w "^COMPLEX_START" $CONFIG | sed 's/COMPLEX_START=//g')
COMPLEX_END=$(grep -w "^COMPLEX_END" $CONFIG | sed 's/COMPLEX_END=//g')
FIVE_PADDING=$(grep -w "^FIVE_PADDING" $CONFIG | sed 's/FIVE_PADDING=//g')
THREE_PADDING=$(grep -w "^THREE_PADDING" $CONFIG | sed 's/THREE_PADDING=//g')
XMAP_FIVE=results/"$SAMPLE"_separate_alignments/"$SAMPLE"_five.xmap
XMAP_THREE=results/"$SAMPLE"_separate_alignments/"$SAMPLE"_three.xmap
QCMAP_FIVE=results/"$SAMPLE"_separate_alignments/"$SAMPLE"_five_q.cmap
QCMAP_THREE=results/"$SAMPLE"_separate_alignments/"$SAMPLE"_three_q.cmap
XMAP_SPANNING=results/"$SAMPLE"_separate_alignments/"$SAMPLE"_completely_spanning.xmap
QCMAP_SPANNING=results/"$SAMPLE"_separate_alignments/"$SAMPLE"_completely_spanning_q.cmap

# add padding to start and end:
COMPLEX_START=$(echo $COMPLEX_START | awk -v p="$FIVE_PADDING" '{print $1-p}')
COMPLEX_END=$(echo $COMPLEX_END | awk -v p="$THREE_PADDING" '{print $1+p}')

# define SD(s)/complexity
chr=1
complex_ambiguous_region_five=$FIVE_PADDING
complex_ambiguous_region_three=$(echo $COMPLEX_START $COMPLEX_END $THREE_PADDING | tr ' ' '\t' | awk '{print ($2-$1)-$3}')

RCMAP=results/"$SAMPLE"_complex_alignments/"$SAMPLE"_chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END"_r.cmap

# find ref_labels corresponding to complex region coordinates:
first_rcmap_label=$(awk -v c="$chr" \
-v five="$complex_ambiguous_region_five" \
-v three="$complex_ambiguous_region_three" \
'$1 == c && $6 >= five && $6 <= three' "$RCMAP" | cut -f4 | head -1)

last_rcmap_label=$(awk -v c="$chr" \
-v five="$complex_ambiguous_region_five" \
-v three="$complex_ambiguous_region_three" \
'$1 == c && $6 >= five && $6 <= three' "$RCMAP" | cut -f4 | tail -1)

mkdir -p "$SAMPLE"_tmp_extract_anchors

# obtain query labels for 5' anchored molecules
while read line
do

    # get tuples first query label only:
    first_q=$(echo $line | grep -v "#" | tr ' ' '\t' | cut -f14 | \
    tr ')' '\n' | sed 's/(//g' | awk 'NF' | tr ',' '\t' | \
    awk -v fir="$first_rcmap_label" \
    '$1 < fir' | \
    cut -f2 | head -1);

    # repeat for last one:
    last_q=$(echo $line | grep -v "#" | tr ' ' '\t' | cut -f14 | \
    tr ')' '\n' | sed 's/(//g' | awk 'NF' | tr ',' '\t' | \
    awk -v fir="$first_rcmap_label" \
    '$1 < fir' | \
    cut -f2 | tail -1);

    # if first_q < last_q, print integers from 1 to last_q
    if [ $first_q -lt $last_q ];
    then
        for i in $( seq 1 $last_q ); do echo $i | tr '\n' ','; done | sed 's/,$/\t/g';
    # else, go into qcmap and find last label, print integers from last_q to this last label 
    else
        tmp_id=$(echo $line | grep -v "#" | tr ' ' '\t' | cut -f2);
        last_l=$(awk -v idfind="$tmp_id" '$1 == idfind' $QCMAP_FIVE | head -1 | cut -f3);
        for i in $( seq $last_q $last_l ); do echo $i; done | sort -k1nr | tr '\n' ',' | sed 's/,$/\t/g';
    fi

    # get molecule ID:
    echo $line | grep -v "#" | tr ' ' '\t' | cut -f2;
done <$XMAP_FIVE | sort -k2n | uniq | awk -F'\t' 'NF > 1' | awk '$1 != "1"' > "$SAMPLE"_tmp_extract_anchors/"$SAMPLE"_five_query_labels_in_anchor

# obtain query labels for 3' anchored molecules
while read line
do

    # get tuples first query label only:
    first_q=$(echo $line | grep -v "#" | tr ' ' '\t' | cut -f14 | \
    tr ')' '\n' | sed 's/(//g' | awk 'NF' | tr ',' '\t' | \
    awk -v thr="$last_rcmap_label" \
    '$1 > thr' | \
    cut -f2 | head -1);

    # repeat for last one:
    last_q=$(echo $line | grep -v "#" | tr ' ' '\t' | cut -f14 | \
    tr ')' '\n' | sed 's/(//g' | awk 'NF' | tr ',' '\t' | \
    awk -v thr="$last_rcmap_label" \
    '$1 > thr' | \
    cut -f2 | tail -1);

    # if first_q > last_q, print integers from 1 to first_q then rev:
    if [ $first_q -gt $last_q ];
    then
        for i in $( seq 1 $first_q ); do echo $i; done | sort -k1nr | tr '\n' ',' | sed 's/,$/\t/g';
    # else, go into qcmap and find last label, print integers from last_q to this last label 
    else
        tmp_id=$(echo $line | grep -v "#" | tr ' ' '\t' | cut -f2);
        last_l=$(awk -v idfind="$tmp_id" '$1 == idfind' $QCMAP_THREE | head -1 | cut -f3);
        for i in $( seq $first_q $last_l ); do echo $i; done | tr '\n' ',' | sed 's/,$/\t/g';
    fi

    # get molecule ID:
    echo $line | grep -v "#" | tr ' ' '\t' | cut -f2;
done <$XMAP_THREE | sort -k2n | uniq | awk -F'\t' 'NF > 1' | awk '$1 != "1"' > "$SAMPLE"_tmp_extract_anchors/"$SAMPLE"_three_query_labels_in_anchor

# obtain query labels for Spanning molecules
# 5' first
while read line
do

    # get tuples first query label only:
    first_q=$(echo $line | grep -v "#" | tr ' ' '\t' | cut -f14 | \
    tr ')' '\n' | sed 's/(//g' | awk 'NF' | tr ',' '\t' | \
    awk -v fir="$first_rcmap_label" \
    '$1 < fir' | \
    cut -f2 | head -1);

    # repeat for last one:
    last_q=$(echo $line | grep -v "#" | tr ' ' '\t' | cut -f14 | \
    tr ')' '\n' | sed 's/(//g' | awk 'NF' | tr ',' '\t' | \
    awk -v fir="$first_rcmap_label" \
    '$1 < fir' | \
    cut -f2 | tail -1);

    # if first_q < last_q, print integers from 1 to last_q
    if [ $first_q -lt $last_q ];
    then
        for i in $( seq 1 $last_q ); do echo $i | tr '\n' ','; done | sed 's/,$/\t/g';
    # else, go into qcmap and find last label, print integers from last_q to this last label 
    else
        tmp_id=$(echo $line | grep -v "#" | tr ' ' '\t' | cut -f2);
        last_l=$(awk -v idfind="$tmp_id" '$1 == idfind' $QCMAP_SPANNING | head -1 | cut -f3);
        for i in $( seq $last_q $last_l ); do echo $i; done | sort -k1nr | tr '\n' ',' | sed 's/,$/\t/g';
    fi

    # get molecule ID:
    echo $line | grep -v "#" | tr ' ' '\t' | cut -f2;
done <$XMAP_SPANNING | sort -k2n | uniq | awk -F'\t' 'NF > 1' | awk '$1 != "1"' > "$SAMPLE"_tmp_extract_anchors/"$SAMPLE"_spanning_five_query_labels_in_anchor

# 3' next
while read line
do

    # get tuples first query label only:
    first_q=$(echo $line | grep -v "#" | tr ' ' '\t' | cut -f14 | \
    tr ')' '\n' | sed 's/(//g' | awk 'NF' | tr ',' '\t' | \
    awk -v thr="$last_rcmap_label" \
    '$1 > thr' | \
    cut -f2 | head -1);

    # repeat for last one:
    last_q=$(echo $line | grep -v "#" | tr ' ' '\t' | cut -f14 | \
    tr ')' '\n' | sed 's/(//g' | awk 'NF' | tr ',' '\t' | \
    awk -v thr="$last_rcmap_label" \
    '$1 > thr' | \
    cut -f2 | tail -1);

    # if first_q > last_q, print integers from 1 to first_q then rev:
    if [ $first_q -gt $last_q ];
    then
        for i in $( seq 1 $first_q ); do echo $i; done | sort -k1nr | tr '\n' ',' | sed 's/,$/\t/g';
    # else, go into qcmap and find last label, print integers from last_q to this last label 
    else
        tmp_id=$(echo $line | grep -v "#" | tr ' ' '\t' | cut -f2);
        last_l=$(awk -v idfind="$tmp_id" '$1 == idfind' $QCMAP_SPANNING | head -1 | cut -f3);
        for i in $( seq $first_q $last_l ); do echo $i; done | tr '\n' ',' | sed 's/,$/\t/g';
    fi

    # get molecule ID:
    echo $line | grep -v "#" | tr ' ' '\t' | cut -f2;
done <$XMAP_SPANNING | sort -k2n | uniq | awk -F'\t' 'NF > 1' | awk '$1 != "1"' > "$SAMPLE"_tmp_extract_anchors/"$SAMPLE"_spanning_three_query_labels_in_anchor

# combine 5' and 3':
while read line1 line2
do
    awk -v lone="$line1" -v id="$line2" '$2 == id {print lone","$1"\t"$2}' "$SAMPLE"_tmp_extract_anchors/"$SAMPLE"_spanning_three_query_labels_in_anchor
done <"$SAMPLE"_tmp_extract_anchors/"$SAMPLE"_spanning_five_query_labels_in_anchor > "$SAMPLE"_tmp_extract_anchors/"$SAMPLE"_spanning_query_labels_in_anchor


# clean-up:
mkdir -p results/"$SAMPLE"_extract_filter_anchor_labels
mv -f "$SAMPLE"_tmp_extract_anchors/"$SAMPLE"_five_query_labels_in_anchor results/"$SAMPLE"_extract_filter_anchor_labels
mv -f "$SAMPLE"_tmp_extract_anchors/"$SAMPLE"_three_query_labels_in_anchor results/"$SAMPLE"_extract_filter_anchor_labels
mv -f "$SAMPLE"_tmp_extract_anchors/"$SAMPLE"_spanning_query_labels_in_anchor results/"$SAMPLE"_extract_filter_anchor_labels
rm -rf "$SAMPLE"_tmp_extract_anchors


# filtration
ANCHOR_FIVE="$SAMPLE"_five_query_labels_in_anchor
ANCHOR_THREE="$SAMPLE"_three_query_labels_in_anchor
ANCHOR_SPANNING="$SAMPLE"_spanning_query_labels_in_anchor
# number of anchor labels desire a molecule to have mapped to reference
# >= than this number:
#num_anchor_labels=8
num_anchor_labels=$(grep -w "^ANCHOR_LABELS" $CONFIG | sed 's/ANCHOR_LABELS=//g')


# 5'
# obtain number of anchor labels:
cut -f1 results/"$SAMPLE"_extract_filter_anchor_labels/$ANCHOR_FIVE | awk -F',' '{print NF}' > results/"$SAMPLE"_extract_filter_anchor_labels/"$SAMPLE".tmp

# paste to original input:
paste results/"$SAMPLE"_extract_filter_anchor_labels/$ANCHOR_FIVE results/"$SAMPLE"_extract_filter_anchor_labels/"$SAMPLE".tmp > results/"$SAMPLE"_extract_filter_anchor_labels/"$SAMPLE".tmp.anchor

# filter for $num_anchor_labels:
awk -v num="$num_anchor_labels" '$3 >= num' results/"$SAMPLE"_extract_filter_anchor_labels/"$SAMPLE".tmp.anchor | cut -f1-2 > results/"$SAMPLE"_extract_filter_anchor_labels/"$SAMPLE".anchor.filt
rm -f results/"$SAMPLE"_extract_filter_anchor_labels/$ANCHOR_FIVE
mv -f results/"$SAMPLE"_extract_filter_anchor_labels/"$SAMPLE".anchor.filt results/"$SAMPLE"_extract_filter_anchor_labels/$ANCHOR_FIVE

# clean-up:
rm -f results/"$SAMPLE"_extract_filter_anchor_labels/"$SAMPLE".tmp
rm -f results/"$SAMPLE"_extract_filter_anchor_labels/"$SAMPLE".tmp.anchor


# 3'
# obtain number of anchor labels:
cut -f1 results/"$SAMPLE"_extract_filter_anchor_labels/$ANCHOR_THREE | awk -F',' '{print NF}' > results/"$SAMPLE"_extract_filter_anchor_labels/"$SAMPLE".tmp

# paste to original input:
paste results/"$SAMPLE"_extract_filter_anchor_labels/$ANCHOR_THREE results/"$SAMPLE"_extract_filter_anchor_labels/"$SAMPLE".tmp > results/"$SAMPLE"_extract_filter_anchor_labels/"$SAMPLE".tmp.anchor

# filter for $num_anchor_labels:
awk -v num="$num_anchor_labels" '$3 >= num' results/"$SAMPLE"_extract_filter_anchor_labels/"$SAMPLE".tmp.anchor | cut -f1-2 > results/"$SAMPLE"_extract_filter_anchor_labels/"$SAMPLE".anchor.filt
rm -f results/"$SAMPLE"_extract_filter_anchor_labels/$ANCHOR_THREE
mv -f results/"$SAMPLE"_extract_filter_anchor_labels/"$SAMPLE".anchor.filt results/"$SAMPLE"_extract_filter_anchor_labels/$ANCHOR_THREE

# clean-up:
rm -f results/"$SAMPLE"_extract_filter_anchor_labels/"$SAMPLE".tmp
rm -f results/"$SAMPLE"_extract_filter_anchor_labels/"$SAMPLE".tmp.anchor


# Spanning
# obtain number of anchor labels:
cut -f1 results/"$SAMPLE"_extract_filter_anchor_labels/$ANCHOR_SPANNING | awk -F',' '{print NF}' > results/"$SAMPLE"_extract_filter_anchor_labels/"$SAMPLE".tmp

# paste to original input:
paste results/"$SAMPLE"_extract_filter_anchor_labels/$ANCHOR_SPANNING results/"$SAMPLE"_extract_filter_anchor_labels/"$SAMPLE".tmp > results/"$SAMPLE"_extract_filter_anchor_labels/"$SAMPLE".tmp.anchor

# filter for $num_anchor_labels:
awk -v num="$num_anchor_labels" '$3 >= num' results/"$SAMPLE"_extract_filter_anchor_labels/"$SAMPLE".tmp.anchor | cut -f1-2 > results/"$SAMPLE"_extract_filter_anchor_labels/"$SAMPLE".anchor.filt
rm -f results/"$SAMPLE"_extract_filter_anchor_labels/$ANCHOR_SPANNING
mv -f results/"$SAMPLE"_extract_filter_anchor_labels/"$SAMPLE".anchor.filt results/"$SAMPLE"_extract_filter_anchor_labels/$ANCHOR_SPANNING

# clean-up:
rm -f results/"$SAMPLE"_extract_filter_anchor_labels/"$SAMPLE".tmp
rm -f results/"$SAMPLE"_extract_filter_anchor_labels/"$SAMPLE".tmp.anchor



