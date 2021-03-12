#!/bin/sh

#conda activate complex

# Define sample and alignment files:
CONFIG="config"

SAMPLE=$(grep -w "^SAMPLE" $CONFIG | sed 's/SAMPLE=//g')
MATRIX=results/"$SAMPLE"_five_n_three_all_vs_all_scores/"$SAMPLE"_tmp3
OUT="$SAMPLE"_kmeans
K_THRESHOLD=$(grep -w "^K_THRESHOLD" $CONFIG | sed 's/K_THRESHOLD=//g')

# fix blank tab at end:
sed 's/\t$//g' $MATRIX > "$MATRIX".fixed

python scripts/auto_clusters.py -i "$MATRIX".fixed -o ./$OUT -t $K_THRESHOLD

#mkdir -p results/"$SAMPLE"_kmeans
#mv -f $OUT results/"$SAMPLE"_kmeans

