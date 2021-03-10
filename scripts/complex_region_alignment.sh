#!/bin/sh

# requirements: 
# samtools
# RefAligner
# multi-threaded environment


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

# download reference sequence:
reference_link="https://hgdownload.cse.ucsc.edu/goldenpath/"$COMPLEX_GENOME"/chromosomes/chr"$COMPLEX_CHR".fa.gz"
wget $reference_link
gzip -d "chr"$COMPLEX_CHR".fa.gz"

# subset the sequence:
samtools faidx "chr"$COMPLEX_CHR".fa" "chr"$COMPLEX_CHR":"$COMPLEX_START"-"$COMPLEX_END > "chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END".fa"

# in silico nick the fasta file:
perl scripts/fa2cmap_multi_color.pl -i "chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END".fa" -o . -e DLE1 1
mv -f "chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END"_DLE1_0kb_0labels.cmap" "chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END".cmap"

# align molecules to complex region, allowing unlimited multi- and split-maps:
reference="chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END".cmap"

./scripts/RefAligner \
-i $BNX \
-ref $reference \
-o $SAMPLE"_chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END \
-M 3 3 -FP 0.918057 -FN 0.099062 -sf 0.233588 -sd 0.090609 -S 0 \
-minlen 160 -minsites 10 -T 1e-11 -res 3.5 -resSD 0.7 -Mfast 0 -biaswt 0 \
-A 5 -BestRef 0 -nosplit 0 -outlier 1e-7 -endoutlier 1e-7 -f \
-maxmem 32 \
-TotalThreads 16 \
-maxthreads 16

mkdir -p results/$SAMPLE"_complex_alignments"
mv -f $SAMPLE"_chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END".xmap" ./results/$SAMPLE"_complex_alignments"
mv -f $SAMPLE"_chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END"_q.cmap" ./results/$SAMPLE"_complex_alignments"
mv -f $SAMPLE"_chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END"_r.cmap" ./results/$SAMPLE"_complex_alignments"
#mv -f $SAMPLE"_complex_alignments" results/

## clean-up:
rm -f "chr"$COMPLEX_CHR".fa"
rm -f "chr"$COMPLEX_CHR".fa.fai"
rm -f "chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END".fa"
rm -f "chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END"_DLE1_0kb_0labels_key.txt"
rm -f "chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END"_DLE1_0kb_0labels_summary.txt"
rm -f status.txt
rm -f "chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END".cmap"
rm -f $SAMPLE"_chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END".errbin"
rm -f $SAMPLE"_chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END".errbin"
rm -f $SAMPLE"_chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END".err"
rm -f $SAMPLE"_chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END".maprate"
rm -f $SAMPLE"_chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END"_intervals.txt"
rm -f $SAMPLE"_chr"$COMPLEX_CHR"_"$COMPLEX_START"_"$COMPLEX_END".map"




