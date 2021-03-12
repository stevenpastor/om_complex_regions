#!/bin/sh


# requirements: 
# multi-threaded environment


#################
CONFIG="config"

SAMPLE=$(grep -w "^SAMPLE" $CONFIG | sed 's/SAMPLE=//g')
FIVE_N_THREE=results/"$SAMPLE"_long_haps/"$SAMPLE"_five_nested_three_final
COMPLEX_SCORES=results/"$SAMPLE"_proportional_comparisons/"$SAMPLE"_removed_one_sd_below.tsv

COMPLEX_START=$(grep -w "^COMPLEX_START" $CONFIG | sed 's/COMPLEX_START=//g')
COMPLEX_END=$(grep -w "^COMPLEX_END" $CONFIG | sed 's/COMPLEX_END=//g')
FIVE_PADDING=$(grep -w "^FIVE_PADDING" $CONFIG | sed 's/FIVE_PADDING=//g')
THREE_PADDING=$(grep -w "^THREE_PADDING" $CONFIG | sed 's/THREE_PADDING=//g')

# add padding to start and end:
COMPLEX_START=$(echo $COMPLEX_START | awk -v p="$FIVE_PADDING" '{print $1-p}')
COMPLEX_END=$(echo $COMPLEX_END | awk -v p="$THREE_PADDING" '{print $1+p}')

# define SD(s)/complexity
chr=1
complex_ambiguous_region_five=$FIVE_PADDING
complex_ambiguous_region_three=$(echo $COMPLEX_START $COMPLEX_END $THREE_PADDING | tr ' ' '\t' | awk '{print ($2-$1)-$3}')
#################



#cut -f1 $FIVE_N_THREE | sort | uniq | sort -k1n > ids
#
#mkdir -p results/"$SAMPLE"_post_five_n_three
#
#while read i
#do
#grep -w "$i" $FIVE_N_THREE > $i
#
#while read line1 line2 line3
#do 
#grep "#" results/"$SAMPLE"_complex_alignments/"$SAMPLE"_chr22_17800000_19400000_q.cmap > $line1.$line2.$line3.cmap

#grep -v "#" results/"$SAMPLE"_complex_alignments/"$SAMPLE"_chr22_17800000_19400000_q.cmap | awk -v id="$line1" '$1 == id' >> $line1.$line2.$line3.cmap; 
#grep -v "#" results/"$SAMPLE"_complex_alignments/"$SAMPLE"_chr22_17800000_19400000_q.cmap | awk -v id="$line2" '$1 == id' >> $line1.$line2.$line3.cmap; 
#grep -v "#" results/"$SAMPLE"_complex_alignments/"$SAMPLE"_chr22_17800000_19400000_q.cmap | awk -v id="$line3" '$1 == id' >> $line1.$line2.$line3.cmap; 
#
#./scripts/RefAligner -pairmerge -pairmergeRepeat -i $line1.$line2.$line3.cmap -o $line1.$line2.$line3.merged -f; done <$i
#
#for d in *merged_contig*cmap; do echo $d; done > input
#
#./scripts/RefAligner -merge -if input -o $i.merged.all
#
#./scripts/RefAligner -i $i.merged.all.cmap -ref results/"$SAMPLE"_complex_alignments/"$SAMPLE"_chr22_17800000_19400000_r.cmap -o "$i".merged.vs.ref -M 3 3 -FP 0.918057 -FN 0.099062 -sf 0.233588 -sd 0.090609 -S 0 -minlen 160 -minsites 10 -T 1e-11 -res 3.5 -resSD 0.7 -Mfast 0 -biaswt 0 -A 5 -BestRef 0 -nosplit 0 -outlier 1e-7 -endoutlier 1e-7 -f
#
#rm -f $i
#
#mv -f "$i".merged.all.cmap results/"$SAMPLE"_post_five_n_three
#mv -f "$i".merged.vs.ref*  results/"$SAMPLE"_post_five_n_three
#
#rm -f *cmap
#rm -f input
#rm -f *idmap
#
#done <ids
#
#
#for x in results/"$SAMPLE"_post_five_n_three/*xmap
#do
#
#xsam=$(echo $x | cut -d'/' -f3 | sed 's/\.merged\.vs\.ref\.xmap//g')
#
## filter for molecules mapped to both 5' and 3' anchors:
## first, find ids 5' of complexity (even way upstream):
#grep -v "#" "$x" | \
#awk -v st="$complex_ambiguous_region_five" \
#-v fl="$flanking_complex" \
#'$6 <= st && st-$6 >= fl' | cut -f2 | sort -k1n | uniq \
#> results/"$SAMPLE"_post_five_n_three/"$SAMPLE"_five_ids
#
## since considering split-mapped molecules, have to 
## group ids and loop through each's ref coords to 
## determine if map beyond anchor (into complexity):
#while read line
#do
#    awk -v id="$line" '$2 == id' "$x" | cut -f1-7 | \
#    while read match
#    do
#        echo $match | tr ' ' '\t' | awk -v st="$complex_ambiguous_region_five" '$7 > st'
#    done <"${1:-/dev/stdin}"
#done <results/"$SAMPLE"_post_five_n_three/"$SAMPLE"_five_ids | cut -f2 | sort > results/"$SAMPLE"_post_five_n_three/"$SAMPLE"_five_ids_anchor_complex
#
## go back and obtain all lines with these ids:
#while read line
#do
#    awk -v id="$line" '$2 == id' "$x"
#done <results/"$SAMPLE"_post_five_n_three/"$SAMPLE"_five_ids_anchor_complex | sort -k1n | uniq | sort -k1n > results/"$SAMPLE"_post_five_n_three/"$SAMPLE"_five.xmap
#
## find 3'-anchored molecules (normal and split-mapped):
#grep -v "#" "$x" | \
#awk -v en="$complex_ambiguous_region_three" \
#-v fl="$flanking_complex" \
#'$7 >= en && $7-en >= fl' | cut -f2 | sort -k1n | uniq \
#> results/"$SAMPLE"_post_five_n_three/"$SAMPLE"_three_ids
#
## see above 5' for logic behind this:
#while read line
#do
#    awk -v id="$line" '$2 == id' "$x" | cut -f1-7 | \
#    while read match
#    do
#        echo $match | tr ' ' '\t' | awk -v en="$complex_ambiguous_region_three" '$6 < en'
#    done <"${1:-/dev/stdin}"
#done <results/"$SAMPLE"_post_five_n_three/"$SAMPLE"_three_ids | cut -f2 | sort > results/"$SAMPLE"_post_five_n_three/"$SAMPLE"_three_ids_anchor_complex
#
## go back and obtain all lines with these ids:
#while read line
#do
#    awk -v id="$line" '$2 == id' "$x"
#done <results/"$SAMPLE"_post_five_n_three/"$SAMPLE"_three_ids_anchor_complex | sort -k1n | uniq | sort -k1n > results/"$SAMPLE"_post_five_n_three/"$SAMPLE"_three.xmap
#
## find completely spanning molecules:
#cut -f2 results/"$SAMPLE"_post_five_n_three/"$SAMPLE"_five.xmap | sort | uniq | sort > results/"$SAMPLE"_post_five_n_three/"$SAMPLE"_five_ids_tojoin
#cut -f2 results/"$SAMPLE"_post_five_n_three/"$SAMPLE"_three.xmap | sort | uniq | sort > results/"$SAMPLE"_post_five_n_three/"$SAMPLE"_three_ids_tojoin
#comm -12 results/"$SAMPLE"_post_five_n_three/"$SAMPLE"_five_ids_tojoin results/"$SAMPLE"_post_five_n_three/"$SAMPLE"_three_ids_tojoin > results/"$SAMPLE"_post_five_n_three/"$SAMPLE"_completely_spanning_ids
#while read line
#do
#    awk -v id="$line" '$2 == id' "$x"
#done <results/"$SAMPLE"_post_five_n_three/"$SAMPLE"_completely_spanning_ids | sort -k1n | uniq | sort -k1n > results/"$SAMPLE"_post_five_n_three/"$SAMPLE"_"$xsam"_completely_spanning.xmap
#
#done


# QCMAPs:
# obtain qcmaps based on group ids:
#for q in results/"$SAMPLE"_post_five_n_three/*.merged.vs.ref_q.cmap
#do
#
#rm -f results/"$SAMPLE"_post_five_n_three/"$SAMPLE"_"$qsam"_completely_spanning_q.cmap
#
#qsam=$(echo $q | cut -d'/' -f3 | sed 's/\.merged\.vs\.ref_q\.cmap//g')
#
## obtain qcmaps based on group ids:
#grep "#" $q > results/"$SAMPLE"_post_five_n_three/"$SAMPLE"_"$qsam"_completely_spanning_q.cmap
#
#grep -v "#" results/"$SAMPLE"_post_five_n_three/"$SAMPLE"_"$qsam"_completely_spanning.xmap | cut -f2 | sort | uniq > results/"$SAMPLE"_post_five_n_three/"$SAMPLE".uniq.spanning.ids
#
#while read line
#do
#    grep -v "#" $q | awk -v id="$line" '$1 == id'
#done <results/"$SAMPLE"_post_five_n_three/"$SAMPLE".uniq.spanning.ids >> results/"$SAMPLE"_post_five_n_three/"$SAMPLE"_"$qsam"_completely_spanning_q.cmap
#
#done


## merge into as few as possible:
#for m in results/"$SAMPLE"_post_five_n_three/*_completely_spanning_q.cmap
#do
#    sam=$(echo $m | sed 's/_completely_spanning_q\.cmap//g')
#    ./scripts/RefAligner -pairmerge -pairmergeRepeat -i $m -o "$sam".mostmerged
#
##    for c in "$sam".mostmerged_contig*.cmap
##    do
##        ./scripts/RefAligner -i $c -ref results/"$SAMPLE"_complex_alignments/"$SAMPLE"_chr22_17800000_19400000_r.cmap -o "$c".vsref -M 3 3 -FP 0.918057 -FN 0.099062 -sf 0.233588 -sd 0.090609 -S 0 -minlen 160 -minsites 10 -T 1e-11 -res 3.5 -resSD 0.7 -Mfast 0 -biaswt 0 -A 5 -BestRef 0 -nosplit 0 -outlier 1e-7 -endoutlier 1e-7 -f
##    done
#
#done


########## HERE:
# Take each one and align to all others except self (iterate and make new query and ref CMAPs each time); again, filter for only alignments where both 5’ and 3’ anchors mapped (>=X again) and can make a matrix from this and k means:
for mm in results/"$SAMPLE"_post_five_n_three/*mostmerged_contig*.cmap
do
    echo $mm
done > input

awk '{ a[$0] } END { for (i in a){ for (j in a){print (i "\t" j) } } }' input | awk '$1 != $2' > all_combos

cut -f1 all_combos | sort | uniq | sort -k1n | while read line; do echo $line | tr '\n' '\t'; awk -v id="$line" '$1 == id {print $2}' all_combos | tr '\n' '\t'; echo ""; done <"${1:-/dev/stdin}" | sed 's/,$//g' > all_combos_unique

# all_combos_unique has 13 rows (one per unique merged cmap) and cols2-end are the other merged cmaps
# per row, can make one cmap of col1, and another merged cmap of cols2-end
# refaligner these 2 cmaps with no splitmapping and bestref turned on
# obtain the highest confidence score in resulting XMAP and get the IDs (query and ref, cols2,3 in XMAP)
# do they make sense? can I do something with this?

##########


# clean-up:
rm -f ids
rm -f input
rm -f all_combos
rm -f all_combos_unique

