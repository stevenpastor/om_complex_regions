#!/bin/sh

# requires 
# RefAligner
# multi-threaded environment


#################
# nested molecules may map to more than one reference locus
# take them and ensure they uniquely map to defined 
# complex region and not another one outside your
# initial definition
# example, 160kb in LCR22A and LCR22D
# removes molecules mapping to either one
# as unhelpful
#################


# Define sample and alignment files:
CONFIG="config"

SAMPLE=$(grep -w "^SAMPLE" $CONFIG | sed 's/SAMPLE=//g')
COMPLEX_CHR=$(grep -w "^COMPLEX_CHR" $CONFIG | sed 's/COMPLEX_CHR=//g')
COMPLEX_START=$(grep -w "^COMPLEX_START" $CONFIG | sed 's/COMPLEX_START=//g')
COMPLEX_END=$(grep -w "^COMPLEX_END" $CONFIG | sed 's/COMPLEX_END=//g')
FIVE_PADDING=$(grep -w "^FIVE_PADDING" $CONFIG | sed 's/FIVE_PADDING=//g')
THREE_PADDING=$(grep -w "^THREE_PADDING" $CONFIG | sed 's/THREE_PADDING=//g')
WG_REFERENCE=$(grep -w "^WG_REFERENCE" $CONFIG | sed 's/WG_REFERENCE=//g')
XMAP_NESTED="$SAMPLE"_nested.xmap
QCMAP_NESTED="$SAMPLE"_nested_q.cmap
#XMAP_NESTED=results/"$SAMPLE"_separate_alignments/"$SAMPLE"_nested.xmap
#QCMAP_NESTED=results/"$SAMPLE"_separate_alignments/"$SAMPLE"_nested_q.cmap

# add padding to start and end:
COMPLEX_START=$(echo $COMPLEX_START | awk -v p="$FIVE_PADDING" '{print $1-p}')
COMPLEX_END=$(echo $COMPLEX_END | awk -v p="$THREE_PADDING" '{print $1+p}')

mkdir -p "$SAMPLE"_fix_nested_tmp_results
gzip -d $WG_REFERENCE
WG_REFERENCE=$(echo $WG_REFERENCE | sed 's/\.gz//g')

./scripts/RefAligner \
-i results/"$SAMPLE"_separate_alignments/$QCMAP_NESTED \
-ref $WG_REFERENCE \
-o "$SAMPLE"_fix_nested_tmp_results/$SAMPLE"_nested_WG" \
-M 3 3 -FP 0.918057 -FN 0.099062 -sf 0.233588 -sd 0.090609 -S 0 \
-minlen 160 -minsites 10 -T 1e-11 -res 3.5 -resSD 0.7 -Mfast 0 -biaswt 0 \
-A 5 -BestRef 0 -nosplit 0 -outlier 1e-7 -endoutlier 1e-7 -f -maxmem 16 -maxthreads 8

# find molecules which have same aligned loci:
ADD_ME=$COMPLEX_START
# contig id, refstart, refend from WG alignment:
grep -v "#" "$SAMPLE"_fix_nested_tmp_results/$SAMPLE"_nested_WG.xmap" | cut -f2,6,7 > "$SAMPLE"_fix_nested_tmp_results/$SAMPLE"_nested_WG.xmap.filt"
# remove decimals (no floats):
sed 's/\.[0-9]//g' "$SAMPLE"_fix_nested_tmp_results/$SAMPLE"_nested_WG.xmap.filt" > "$SAMPLE"_fix_nested_tmp_results/$SAMPLE"_nested_WG.xmap.tmp"
# replace float version with non-float version:
rm -f "$SAMPLE"_fix_nested_tmp_results/$SAMPLE"_nested_WG.xmap.filt"
mv -f "$SAMPLE"_fix_nested_tmp_results/$SAMPLE"_nested_WG.xmap.tmp" "$SAMPLE"_fix_nested_tmp_results/$SAMPLE"_nested_WG.xmap.filt"

# contig id, refstart, refend from complex region specific alignment:
grep -v "#" results/"$SAMPLE"_separate_alignments/$XMAP_NESTED | cut -f2,6,7 > "$SAMPLE"_fix_nested_tmp_results/$XMAP_NESTED".filt"
# adding complex start padded to these coordinates to normalize to WG refstart and refend coordinates:
awk -v add="$ADD_ME" '{print $1"\t"$2+add-1"\t"$3+add-1}' "$SAMPLE"_fix_nested_tmp_results/$XMAP_NESTED".filt" > "$SAMPLE"_fix_nested_tmp_results/$XMAP_NESTED".tmp"
rm -f "$SAMPLE"_fix_nested_tmp_results/$XMAP_NESTED".filt"
mv -f "$SAMPLE"_fix_nested_tmp_results/$XMAP_NESTED".tmp" "$SAMPLE"_fix_nested_tmp_results/$XMAP_NESTED".filt"

while read line
do
    grep -w "$line" "$SAMPLE"_fix_nested_tmp_results/$SAMPLE"_nested_WG.xmap.filt"
done <"$SAMPLE"_fix_nested_tmp_results/$XMAP_NESTED".filt" | cut -f1 | sort | uniq | sort -k1n > "$SAMPLE"_fix_nested_tmp_results/$SAMPLE"_nested_ids_to_get"

grep "#" results/"$SAMPLE"_separate_alignments/$XMAP_NESTED > "$SAMPLE"_fix_nested_tmp_results/$XMAP_NESTED".fix"
while read line
do
    grep -v "#" results/"$SAMPLE"_separate_alignments/$XMAP_NESTED | awk -v id="$line" '$2 == id' >> "$SAMPLE"_fix_nested_tmp_results/$XMAP_NESTED".fix"
done <"$SAMPLE"_fix_nested_tmp_results/$SAMPLE"_nested_ids_to_get"

grep "#" results/"$SAMPLE"_separate_alignments/$QCMAP_NESTED > "$SAMPLE"_fix_nested_tmp_results/$QCMAP_NESTED".fix"
while read line
do
    grep -v "#" results/"$SAMPLE"_separate_alignments/$QCMAP_NESTED | awk -v id="$line" '$1 == id' >> "$SAMPLE"_fix_nested_tmp_results/$QCMAP_NESTED".fix"
done <"$SAMPLE"_fix_nested_tmp_results/$SAMPLE"_nested_ids_to_get"


# clean-up
gzip $WG_REFERENCE
mkdir -p results/"$SAMPLE"_fix_nested
mv -f "$SAMPLE"_fix_nested_tmp_results/$XMAP_NESTED".fix" results/"$SAMPLE"_fix_nested/$XMAP_NESTED
mv -f "$SAMPLE"_fix_nested_tmp_results/$QCMAP_NESTED".fix" results/"$SAMPLE"_fix_nested/$QCMAP_NESTED
rm -rf "$SAMPLE"_fix_nested_tmp_results

