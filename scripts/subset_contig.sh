#!/usr/bin/env bash

CONFIG=config
SAMPLE=$(grep "^SAMPLE" $CONFIG | sed 's/SAMPLE=//g')
CHR=$(grep "^COMPLEX_CHR" $CONFIG | sed 's/COMPLEX_CHR=//g')

if [ "$CHR" == "X" ]
then
    CHR=23
fi
if [ "$CHR" == "Y" ]
then
    CHR=24
fi

COMPLEX_START=$(grep "^COMPLEX_START" $CONFIG | sed 's/COMPLEX_START=//g')
COMPLEX_END=$(grep "^COMPLEX_END" $CONFIG | sed 's/COMPLEX_END=//g')
FIVE_PADDING=$(grep "^FIVE_PADDING" $CONFIG | sed 's/FIVE_PADDING=//g')
THREE_PADDING=$(grep "^THREE_PADDING" $CONFIG | sed 's/THREE_PADDING=//g')

# typical padding:
PADDED_START=$( echo "$COMPLEX_START" | awk -v p="$FIVE_PADDING" '{print $0-p}' )
PADDED_END=$( echo "$COMPLEX_END" | awk -v p="$THREE_PADDING" '{print $0+p}' )

## no padding (only SDs):
#PADDED_START=$( echo "$COMPLEX_START" )
#PADDED_END=$( echo "$COMPLEX_END" )

## 100kb padding:
#PADDED_START=$( echo "$COMPLEX_START" | awk '{print $0-100000}' )
#PADDED_END=$( echo "$COMPLEX_END" | awk '{print $0+100000}' )


for x in results/"$SAMPLE"_initial_genome_check/"$SAMPLE"_fullContig*_molecules.xmap
do
    
    contig_id=$(echo $x | cut -d'/' -f3 | cut -d'_' -f2 | sed 's/fullContig//g')

    REF_LABEL_START=$(grep -v "#" data/"$SAMPLE"_output/contigs/exp_refineFinal1/alignref_final/EXP_REFINEFINAL1_r.cmap | awk -v c="$CHR" -v ps="$PADDED_START" '$1 == c && $6 >= ps' | head -1 | cut -f4)
    REF_LABEL_END=$(grep -v "#" data/"$SAMPLE"_output/contigs/exp_refineFinal1/alignref_final/EXP_REFINEFINAL1_r.cmap | awk -v c="$CHR" -v pe="$PADDED_END" '$1 == c && $6 <= pe' | tail -1 | cut -f4)

    # first contig label <= $REF_LABEL_START:
    contig_id_start_tmp=$(awk -v id="$contig_id" '$2==id' results/"$SAMPLE"_initial_genome_check/"$SAMPLE"_fullContigs.xmap | cut -f14 | sed 's/(//g' | tr ')' '\n' | tr ',' '\t' | awk 'NF' | awk -v rs="$REF_LABEL_START" '$1 <= rs' | tail -1 | cut -f2)
    # first contig label >= $REF_LABEL_END:
    contig_id_end_tmp=$(awk -v id="$contig_id" '$2==id' results/"$SAMPLE"_initial_genome_check/"$SAMPLE"_fullContigs.xmap | cut -f14 | sed 's/(//g' | tr ')' '\n' | tr ',' '\t' | awk 'NF' | awk -v re="$REF_LABEL_END" '$1 >= re' | head -1 | cut -f2)

    # flipping start and end due to - orientation, else if + keep same:
    if [ "$contig_id_start_tmp" -gt "$contig_id_end_tmp" ]; then contig_id_start="$contig_id_end_tmp"; contig_id_end="$contig_id_start_tmp"; fi
    if [ "$contig_id_start_tmp" -lt "$contig_id_end_tmp" ]; then contig_id_start="$contig_id_start_tmp"; contig_id_end="$contig_id_end_tmp"; fi

    # preliminary CMAP: will have to fix the label coordinates and total length after the fact:
    grep "#" results/"$SAMPLE"_initial_genome_check/exp_refineFinal1_contig"$contig_id"_r.cmap > results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".cmap
    grep -v "#" results/"$SAMPLE"_initial_genome_check/exp_refineFinal1_contig"$contig_id"_r.cmap | awk -v s="$contig_id_start" '$4 >= s' | awk -v e="$contig_id_end" '$4 <= e' >> results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".cmap

    # begin subsetting QCMAP to make new CMAP for OMTools grouping

    # get all cols which will not change:
    grep -v "#" results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".cmap | cut -f1 > results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col1
    grep -v "#" results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".cmap | cut -f7- > results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".cols7_end

    # get all cols which will change:
    grep -v "#" results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".cmap | cut -f2 > results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col2
    grep -v "#" results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".cmap | cut -f3 > results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col3
    grep -v "#" results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".cmap | cut -f4 > results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col4
    grep -v "#" results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".cmap | cut -f5 > results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col5
    grep -v "#" results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".cmap | cut -f6 > results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col6

    #First row col6 becomes 20.0:
    head -1 results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col6 | awk -v FS='\t' -v OFS='\t' 'BEGIN {OFS = FS} {$1 = "20.0"; print}' > results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col6.tmp

    #Col6 becomes current row col6-previous row col6 (skips first row printing, which is good):
    awk 'NR==1{p=$1;next} {print $1-p; p=$1} END{print $1-p}' results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col6 | sed \$d >> results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col6.tmp
    
    #Col6 cumulative addition of all rows (current + previous until end):
    awk -v OFMT='%.1f' 'BEGIN {sum=0} {sum= sum+$1; OFMT="%f"; print sum}' results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col6.tmp | awk '{$1=sprintf("%.1f",$1)}7' | tr ' ' '\t' >> results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col6.tmp2
    rm -f results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col6.tmp
    mv -f results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col6.tmp2 results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col6.tmp
    rm -f results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col6
    mv -f results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col6.tmp results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col6

    #Col2 becomes last row col6:
    total_length=$(tail -1 results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col6)
    awk -v s="$total_length" '{print $0"\t"s}' results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col2 | cut -f2 > results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col2.tmp
    rm -f results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col2
    mv -f results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col2.tmp results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col2

    #Col3 becomes total number of non-header rows but subtract 1:
    total_labels=$(wc -l results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col3 | cut -d' ' -f1 | awk '{print $0-1}')
    awk -v s="$total_labels" '{print $0"\t"s}' results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col3 | cut -f2 > results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col3.tmp
    rm -f results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col3
    mv -f results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col3.tmp results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col3

    #Col4 becomes 1..total number of non-header rows iteration:
    add_one_total_labels=$(echo "$total_labels" | awk '{print $0+1}')
    for (( c=1; c<=$add_one_total_labels; c++ )); do echo $c; done > results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col4.tmp
    rm -f results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col4
    mv -f results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col4.tmp results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col4

    #Last row col5 becomes 0:
    sed \$d results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col5 > results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col5.tmp
    echo "0" >> results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col5.tmp
    rm -f results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col5
    mv -f results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col5.tmp results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col5

    # paste all together:
    paste results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col1 results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col2 results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col3 results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col4 results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col5 results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".col6 results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".cols7_end > results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".pasted

    # header:
    grep "#" results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".cmap > results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END"_final.cmap

    # cat into header:
    cat results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END".pasted >> results/"$SAMPLE"_initial_genome_check/"$contig_id"_subset_"$PADDED_START"_"$PADDED_END"_final.cmap

done	

