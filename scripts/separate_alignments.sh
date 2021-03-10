#!/bin/sh

# requirements: 
# RefAligner


#################
# Aligns molecules to a complex region
# requires RefAligner, config, fa2cmap_multi_color.pl (perl)

CONFIG="config"

SAMPLE=$(grep -w "^SAMPLE" $CONFIG | sed 's/SAMPLE=//g')
BNX=$(grep -w "^BNX" $CONFIG | sed 's/BNX=//g')
COMPLEX_GENOME=$(grep -w "^COMPLEX_GENOME" $CONFIG | sed 's/COMPLEX_GENOME=//g')
COMPLEX_CHR=$(grep -w "^COMPLEX_CHR" $CONFIG | sed 's/COMPLEX_CHR=//g')
COMPLEX_START=$(grep -w "^COMPLEX_START" $CONFIG | sed 's/COMPLEX_START=//g')
COMPLEX_END=$(grep -w "^COMPLEX_END" $CONFIG | sed 's/COMPLEX_END=//g')
FIVE_PADDING=$(grep -w "^FIVE_PADDING" $CONFIG | sed 's/FIVE_PADDING=//g')
THREE_PADDING=$(grep -w "^THREE_PADDING" $CONFIG | sed 's/THREE_PADDING=//g')
#################



# add padding to start and end:
COMPLEX_START=$(echo $COMPLEX_START | awk -v p="$FIVE_PADDING" '{print $1-p}')
COMPLEX_END=$(echo $COMPLEX_END | awk -v p="$THREE_PADDING" '{print $1+p}')

# separate alignments
XMAP="results/"$SAMPLE"_complex_alignments/"$SAMPLE"_chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END".xmap"
QCMAP="results/"$SAMPLE"_complex_alignments/"$SAMPLE"_chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END"_q.cmap"

# define SD(s)/complexity
chr=1
complex_ambiguous_region_five=$FIVE_PADDING
complex_ambiguous_region_three=$(echo $COMPLEX_START $COMPLEX_END $THREE_PADDING | tr ' ' '\t' | awk '{print ($2-$1)-$3}')

# distance outside complexity
# keep in mind different regions have different amounts of 
# labels:
#flanking_complex=100000
flanking_complex=$(grep -w "^FLANKING_COMPLEX" $CONFIG | sed 's/FLANKING_COMPLEX=//g')

mkdir -p $SAMPLE"_separate_alignments_tmp_results"

# find 5'-anchored molecules (normal and split-mapped):
# first, find ids 5' of complexity (even way upstream):
grep -v "#" "$XMAP" | \
awk -v st="$complex_ambiguous_region_five" \
-v fl="$flanking_complex" \
'$6 <= st && st-$6 >= fl' | cut -f2 | sort -k1n | uniq \
> "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_five_ids

# since considering split-mapped molecules, have to 
# group ids and loop through each's ref coords to 
# determine if map beyond anchor (into complexity):
while read line
do
    awk -v id="$line" '$2 == id' "$XMAP" | cut -f1-7 | \
    while read match
    do
        echo $match | tr ' ' '\t' | awk -v st="$complex_ambiguous_region_five" '$7 > st'
    done <"${1:-/dev/stdin}"
done <"$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_five_ids | cut -f2 | sort > "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_five_ids_anchor_complex

# go back and obtain all lines with these ids:
while read line
do
    awk -v id="$line" '$2 == id' "$XMAP"
done <"$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_five_ids_anchor_complex | sort -k1n | uniq | sort -k1n > "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_five.xmap


# find 3'-anchored molecules (normal and split-mapped):
grep -v "#" "$XMAP" | \
awk -v en="$complex_ambiguous_region_three" \
-v fl="$flanking_complex" \
'$7 >= en && $7-en >= fl' | cut -f2 | sort -k1n | uniq \
> "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_three_ids

# see above 5' for logic behind this:
while read line
do
    awk -v id="$line" '$2 == id' "$XMAP" | cut -f1-7 | \
    while read match
    do
        echo $match | tr ' ' '\t' | awk -v en="$complex_ambiguous_region_three" '$6 < en'
    done <"${1:-/dev/stdin}"
done <"$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_three_ids | cut -f2 | sort > "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_three_ids_anchor_complex

# go back and obtain all lines with these ids:
while read line
do
    awk -v id="$line" '$2 == id' "$XMAP"
done <"$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_three_ids_anchor_complex | sort -k1n | uniq | sort -k1n > "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_three.xmap


# find completely spanning molecules:
cut -f2 "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_five.xmap | sort | uniq | sort > "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_five_ids_tojoin
cut -f2 "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_three.xmap | sort | uniq | sort > "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_three_ids_tojoin
comm -12 "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_five_ids_tojoin "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_three_ids_tojoin > "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_completely_spanning_ids
while read line
do
    awk -v id="$line" '$2 == id' "$XMAP"
done <"$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_completely_spanning_ids | sort -k1n | uniq | sort -k1n > "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_completely_spanning.xmap


# separate 5' only from completely spanning:
grep -v -F -x -f "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_completely_spanning_ids "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_five_ids_anchor_complex > "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_five_ids.tmp
while read line
do
    awk -v id="$line" '$2 == id' "$XMAP"
done <"$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_five_ids.tmp | sort -k1n | uniq | sort -k1n > "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_five.xmap

# separate 3' only from completely spanning:
grep -v -F -x -f "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_completely_spanning_ids "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_three_ids_anchor_complex > "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_three_ids.tmp
while read line
do
    awk -v id="$line" '$2 == id' "$XMAP"
done <"$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_three_ids.tmp | sort -k1n | uniq | sort -k1n > "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_three.xmap


# find nested molecules:
grep -v "#" "$XMAP" | \
awk -v st="$complex_ambiguous_region_five" \
-v en="$complex_ambiguous_region_three" \
'$6 >= st && $7 <= en' | cut -f2 | sort -k1n | uniq \
> "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_nested_ids

# separate from completely spanning:
grep -v -F -x -f "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_completely_spanning_ids "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_nested_ids > "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_nested_ids.tmp
rm -f "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_nested_ids
mv -f "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_nested_ids.tmp "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_nested_ids

while read line
do
    keep=0;

    awk -v id="$line" '$2 == id' "$XMAP" | cut -f2,6,7 | \
    while read match
    do
        # echo $match | tr ' ' '\t' | awk -v en="$complex_ambiguous_region_three" '$6 < en'
        echo $match | tr ' ' '\t' | awk -v st="$complex_ambiguous_region_five" \
        -v en="$complex_ambiguous_region_three" -v k="$keep" '{if ($2 < st || $3 > en) {print $1"\t"k+1} else {print $1"\t"k}}' | tr '\n' '\t';
    done <"${1:-/dev/stdin}"
    echo "";
done <"$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_nested_ids | awk '{$1=$1};1' | tr ' ' '\t' > "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_nested_ids_tmp

# currently handles 3x split-mappers (add more if want):
awk -F'\t' '{ if (NF==2 && $2=="0") {print $1} else if (NF==4 && $2=="0" && $4=="0") {print $1} else if (NF==6 && $2=="0" && $4=="0" && $6 == "0") {print $1} else if (NF==8 && $2=="0" && $4=="0" && $6 == "0" && $8 == "0") {print $1} }' "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_nested_ids_tmp > "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_nested_ids_fixed

while read line
do
    awk -v id="$line" '$2 == id' "$XMAP"
done <"$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_nested_ids_fixed | sort -k1n | uniq | sort -k1n > "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_nested.xmap

mkdir -p results/"$SAMPLE"_separate_alignments

# obtain qcmaps based on group ids:
grep "#" $QCMAP > results/"$SAMPLE"_separate_alignments/"$SAMPLE"_five_q.cmap
grep "#" $QCMAP > results/"$SAMPLE"_separate_alignments/"$SAMPLE"_three_q.cmap
grep "#" $QCMAP > results/"$SAMPLE"_separate_alignments/"$SAMPLE"_nested_q.cmap
grep "#" $QCMAP > results/"$SAMPLE"_separate_alignments/"$SAMPLE"_completely_spanning_q.cmap

grep -v "#" "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_five.xmap | cut -f2 | sort | uniq > "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE".uniq.five.ids
grep -v "#" "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_three.xmap | cut -f2 | sort | uniq > "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE".uniq.three.ids
grep -v "#" "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_nested.xmap | cut -f2 | sort | uniq > "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE".uniq.nested.ids
grep -v "#" "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_completely_spanning.xmap | cut -f2 | sort | uniq > "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE".uniq.spanning.ids

while read line
do
    grep -v "#" $QCMAP | awk -v id="$line" '$1 == id'
done <"$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE".uniq.five.ids >> results/"$SAMPLE"_separate_alignments/"$SAMPLE"_five_q.cmap

while read line
do
    grep -v "#" $QCMAP | awk -v id="$line" '$1 == id'
done <"$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE".uniq.three.ids >> results/"$SAMPLE"_separate_alignments/"$SAMPLE"_three_q.cmap

while read line
do
    grep -v "#" $QCMAP | awk -v id="$line" '$1 == id'
done <"$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE".uniq.nested.ids >> results/"$SAMPLE"_separate_alignments/"$SAMPLE"_nested_q.cmap

while read line
do
    grep -v "#" $QCMAP | awk -v id="$line" '$1 == id'
done <"$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE".uniq.spanning.ids >> results/"$SAMPLE"_separate_alignments/"$SAMPLE"_completely_spanning_q.cmap

mv -f "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_five.xmap results/"$SAMPLE"_separate_alignments/
mv -f "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_three.xmap results/"$SAMPLE"_separate_alignments/
mv -f "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_nested.xmap results/"$SAMPLE"_separate_alignments/
mv -f "$SAMPLE"_separate_alignments_tmp_results/"$SAMPLE"_completely_spanning.xmap results/"$SAMPLE"_separate_alignments/

# clean-up:
rm -rf "$SAMPLE"_separate_alignments_tmp_results

