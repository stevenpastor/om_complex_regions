#!/bin/sh

# requirements: 
# multi-threaded environment


#################
CONFIG="config"

SAMPLE=$(grep -w "^SAMPLE" $CONFIG | sed 's/SAMPLE=//g')
COMPLEX_CHR=$(grep -w "^COMPLEX_CHR" $CONFIG | sed 's/COMPLEX_CHR=//g')
COMPLEX_START=$(grep -w "^COMPLEX_START" $CONFIG | sed 's/COMPLEX_START=//g')
COMPLEX_END=$(grep -w "^COMPLEX_END" $CONFIG | sed 's/COMPLEX_END=//g')
FIVE_PADDING=$(grep -w "^FIVE_PADDING" $CONFIG | sed 's/FIVE_PADDING=//g')
THREE_PADDING=$(grep -w "^THREE_PADDING" $CONFIG | sed 's/THREE_PADDING=//g')
# add padding to start and end:
COMPLEX_START=$(echo $COMPLEX_START | awk -v p="$FIVE_PADDING" '{print $1-p}')
COMPLEX_END=$(echo $COMPLEX_END | awk -v p="$THREE_PADDING" '{print $1+p}')

K_MEANS="$SAMPLE"_kmeans
#################


# use k means file (typically 2 lines) to obtain final molecules for haplotypes

num=1
while read line
do

    # start clean:
    rm -f "$SAMPLE"_chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END"_"$num".xmap

    # place line into new file with new delimiter:
    echo $line | tr ' ' '\n' | tr '_' '\t' > fixed_trav

    # use new file for traversal:
    while read h1 h2 h3
    do
        grep -v "#" results/416477_complex_alignments/416477_chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END".xmap | awk -v id="$h1" '$2 == id' >> "$SAMPLE"_chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END"_"$num".xmap
        grep -v "#" results/416477_complex_alignments/416477_chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END".xmap | awk -v id="$h2" '$2 == id' >> "$SAMPLE"_chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END"_"$num".xmap
        grep -v "#" results/416477_complex_alignments/416477_chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END".xmap | awk -v id="$h3" '$2 == id' >> "$SAMPLE"_chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END"_"$num".xmap
    done <fixed_trav

    # add header and obtain unique lines only (redundancy from above due to connections):
    grep "#" results/416477_complex_alignments/416477_chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END".xmap > "$SAMPLE".tmp
    sort "$SAMPLE"_chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END"_"$num".xmap | uniq | sort >> "$SAMPLE".tmp
    rm -f "$SAMPLE"_chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END"_"$num".xmap
    mv -f "$SAMPLE".tmp "$SAMPLE"_chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END"_"$num".xmap

    # increment for >1 haplotype:
    num=$((num+1)) 

done <$K_MEANS

# copy Q/RCMAP for use in Access (visualize):
cp -f results/416477_complex_alignments/416477_chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END"_q.cmap .
cp -f results/416477_complex_alignments/416477_chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END"_r.cmap .


# clean-up:
rm -f fixed_trav
