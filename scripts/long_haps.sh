#!/bin/sh

CONFIG="config"
SAMPLE=$(grep -w "^SAMPLE" $CONFIG | sed 's/SAMPLE=//g')
MAPPED=results/"$SAMPLE"_proportional_comparisons/"$SAMPLE"_removed_one_sd_below.tsv
# format (5 cols): 54292	three	244192	nested	8.366


# find all 5'-N or N-5' connections:
while read line
do

    echo $line | tr ' ' '\t' | \
    awk '$2 == "five" && $4 == "nested" {print $1"\t"$3}'
    
    echo $line | tr ' ' '\t' | \
    awk '$2 == "nested" && $4 == "five" {print $3"\t"$1}'
    
done <$MAPPED > "$SAMPLE".five

# find all N-N connections:
while read line
do

    echo $line | tr ' ' '\t' | \
    awk '$2 == "nested" && $4 == "nested" {print $1"\t"$3}'
    
done <$MAPPED > "$SAMPLE".nested

# find all 3'-N or N-3' connections:
while read line
do

    echo $line | tr ' ' '\t' | \
    awk '$2 == "three" && $4 == "nested" {print $3"\t"$1}'
    
    echo $line | tr ' ' '\t' | \
    awk '$2 == "nested" && $4 == "three" {print $1"\t"$3}'
    
done <$MAPPED > "$SAMPLE".three

# find all 5'-N-3' or 3'-N-5' by finding common N from above 2 groups:
while read line1 line2
do

    awk -v five="$line1" -v three="$line2" '$1 == three {print five"\t"three"\t"$2}' "$SAMPLE".three

done <"$SAMPLE".five > "$SAMPLE"_five_nested_three

# get all unique 5'-N:
cut -f1,2 "$SAMPLE"_five_nested_three | sort -k1,1n -k2,2n | uniq | sort -k1,1n -k2,2n > "$SAMPLE"_five_nested

# get all unique 5' ids from this file:
cut -f1 "$SAMPLE"_five_nested | sort | uniq | sort -k1n > "$SAMPLE".five.unique

# match 5'-N's N to N-N (either column):
while read line
do

    awk -v five="$line" '$1 == five' "$SAMPLE"_five_nested | cut -f2 | sort | uniq | sort -k1n > "$SAMPLE".five.nesteds

    while read n
    do
        awk -v nested="$n" '{if ($1 == nested) {print $2} else if ($2 == nested) {print $1}}' "$SAMPLE".nested > "$SAMPLE".nested.nesteds
        echo $n >> "$SAMPLE".nested.nesteds
        # is this: "$SAMPLE".five.nesteds the same as this: "$SAMPLE".nested.nesteds
        # if same, then this 5'-N connection is good and kept
        diffs=$(grep -Fxvf "$SAMPLE".nested.nesteds "$SAMPLE".five.nesteds | wc -l)
        #diffs=$(grep -Fxvf "$SAMPLE".five.nesteds "$SAMPLE".nested.nesteds | wc -l)
        if [[ $diffs -eq 0 ]]
        then
            echo $line $n | tr ' ' '\t'
        fi
    done <"$SAMPLE".five.nesteds
 
done <"$SAMPLE".five.unique > "$SAMPLE"_pruned_five_nesteds

# go back and get the 5'-N-3' from these:
while read line1 line2
do
    awk -v five="$line1" -v nested="$line2" '$1 == five && $2 == nested' "$SAMPLE"_five_nested_three
done <"$SAMPLE"_pruned_five_nesteds > "$SAMPLE"_five_nested_three_final


mkdir -p results/"$SAMPLE"_long_haps
mv -f "$SAMPLE"_five_nested_three_final results/"$SAMPLE"_long_haps/


# clean-up
rm -f "$SAMPLE".five
rm -f "$SAMPLE".three
rm -f "$SAMPLE".five.unique
rm -f "$SAMPLE".nested
rm -f "$SAMPLE".five.nesteds
rm -f "$SAMPLE".nested.nesteds
rm -f "$SAMPLE"_five_nested
rm -f "$SAMPLE"_five_nested_three
rm -f "$SAMPLE"_pruned_five_nesteds


