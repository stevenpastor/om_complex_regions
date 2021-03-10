#!/bin/sh


#################
CONFIG="config"

SAMPLE=$(grep -w "^SAMPLE" $CONFIG | sed 's/SAMPLE=//g')
FIVE_QCMAP=results/"$SAMPLE"_filtered_maps/"$SAMPLE"_five_q.cmap
THREE_QCMAP=results/"$SAMPLE"_filtered_maps/"$SAMPLE"_three_q.cmap
NESTED_QCMAP=results/"$SAMPLE"_fix_nested/"$SAMPLE"_nested_q.cmap
SPANNING_QCMAP=results/"$SAMPLE"_separate_alignments/"$SAMPLE"_completely_spanning_q.cmap
#################


# merge:
./scripts/RefAligner -merge -i $FIVE_QCMAP -i $THREE_QCMAP -i $NESTED_QCMAP -i $SPANNING_QCMAP -o "$SAMPLE"_merged

# self-align:
./scripts/RefAligner \
-i "$SAMPLE"_merged.cmap \
-o "$SAMPLE"_merged_self_aligned \
-usecolor 1 -FP 0.44291 -FN 0.08797 -sf 0.06955 -sd 0.02311 -sr 0.01412 -res 3.033 -resSD 0.75 \
-bpp 498.07 -usecolor 1 \
-mres 0.9 -S 0.1 -A 8 -outlier 0.001 -endoutlier 0 \
-RepeatMask 5 0.01 -RepeatRec 0.7 0.6 1.4 -PVres 2 \
-alignscore -hashgen 5 4 2.2 1.2 0.05 3.0 1 1 -hash -MinSD 0.0 \
-minlen 50 -minsites 7 -maxsites 750 -maxSiteDensity 30 -maxContigSiteDensity 40 \
-MinFN 0.1 -MinFP 2.0 -MinSD 0.0 -MinSR 0.02 -MaxSF 0.15 -MaxSD 0.12 -MaxSR 0.03 

# find molecules belonging to certain groups and make idmap to inform next filter script:
grep -v "#" $FIVE_QCMAP | cut -f1 | sort | uniq | sort -k1n | awk '{print $1"\tfive"}' > "$SAMPLE"_merged.idmap
grep -v "#" $THREE_QCMAP | cut -f1 | sort | uniq | sort -k1n | awk '{print $1"\tthree"}' >> "$SAMPLE"_merged.idmap
grep -v "#" $NESTED_QCMAP | cut -f1 | sort | uniq | sort -k1n | awk '{print $1"\tnested"}' >> "$SAMPLE"_merged.idmap
grep -v "#" $SPANNING_QCMAP | cut -f1 | sort | uniq | sort -k1n | awk '{print $1"\tspanning"}' >> "$SAMPLE"_merged.idmap

# clean-up:
mkdir -p results/"$SAMPLE"_merge_into_one_cmap_self_align
mv -f "$SAMPLE"_merged* results/"$SAMPLE"_merge_into_one_cmap_self_align

