#!/bin/sh


# requirements: 
# multi-threaded environment


#################
CONFIG="config"

SAMPLE=$(grep -w "^SAMPLE" $CONFIG | sed 's/SAMPLE=//g')
FIVE_N_THREE=results/"$SAMPLE"_long_haps/"$SAMPLE"_five_nested_three_final
COMPLEX_SCORES=results/"$SAMPLE"_proportional_comparisons/"$SAMPLE"_removed_one_sd_below.tsv
#################

# FIVE_N_THREE:
#154866	20078	1792176

# COMPLEX_SCORES:
#20078	nested	2375746	nested	9.605


mkdir -p results/"$SAMPLE"_five_n_three_all_vs_all_scores

# all combinations of FIVE_N_THREE:
awk '{ a[$0] } END { for (i in a){ for (j in a){print (i "\t" j) } } }' $FIVE_N_THREE > results/"$SAMPLE"_five_n_three_all_vs_all_scores/"$SAMPLE"_tmp
#1086827	479903	197544	1086827	479903	54292

while read line1 line2 line3 line4 line5 line6
do

    #echo "$line1 $line2 $line3 $line4 $line5 $line6" | tr '\n' '\t'
    echo $line1"_"$line2"_"$line3" "$line4"_"$line5"_"$line6 | tr ' ' '\t' | tr '\n' '\t'

    # if 1=4 and 2=5 and 3=6:
    if [ "$line1" -eq "$line4" ] && [ "$line2" -eq "$line5" ] && [ "$line3" -eq "$line6" ]
    then
        echo "100.0 100.0 100.0" | tr ' ' '\t' | awk '{print $1+$2+$3}'
    # if 1=4 and 2=5:
    elif [ "$line1" -eq "$line4" ] && [ "$line2" -eq "$line5" ] && [ "$line3" -ne "$line6" ]
    then
        three=$(grep -w "$line3" $COMPLEX_SCORES | grep -w "$line6" | cut -f5)
        if [ -z "${three}" ]; then three=0.0; fi
        echo "100.0 100.0 $three" | tr ' ' '\t' | awk '{print $1+$2+$3}'
    # if 1=4 and 3=6:
    elif [ "$line1" -eq "$line4" ] && [ "$line2" -ne "$line5" ] && [ "$line3" -eq "$line6" ]
    then
        n=$(grep -w "$line2" $COMPLEX_SCORES | grep -w "$line5" | cut -f5)
        if [ -z "${n}" ]; then n=0.0; fi
        echo "100.0 $n 100.0" | tr ' ' '\t' | awk '{print $1+$2+$3}'
    # if 2=5 and 3=6:
    elif [ "$line1" -ne "$line4" ] && [ "$line2" -eq "$line5" ] && [ "$line3" -eq "$line6" ]
    then
        five=$(grep -w "$line1" $COMPLEX_SCORES | grep -w "$line4" | cut -f5)
        if [ -z "${five}" ]; then five=0.0; fi
        echo "$five 100.0 100.0" | tr ' ' '\t' | awk '{print $1+$2+$3}'
    # if 1=4:
    elif [ "$line1" -eq "$line4" ] && [ "$line2" -ne "$line5" ] && [ "$line3" -ne "$line6" ]
    then
        n=$(grep -w "$line2" $COMPLEX_SCORES | grep -w "$line5" | cut -f5)
        three=$(grep -w "$line3" $COMPLEX_SCORES | grep -w "$line6" | cut -f5)
        if [ -z "${n}" ]; then n=0.0; fi
        if [ -z "${three}" ]; then three=0.0; fi
        echo "100.0 $n $three" | tr ' ' '\t' | awk '{print $1+$2+$3}'
    # if 2=5:
    elif [ "$line1" -ne "$line4" ] && [ "$line2" -eq "$line5" ] && [ "$line3" -ne "$line6" ]
    then
        five=$(grep -w "$line1" $COMPLEX_SCORES | grep -w "$line4" | cut -f5)
        three=$(grep -w "$line3" $COMPLEX_SCORES | grep -w "$line6" | cut -f5)
        if [ -z "${five}" ]; then five=0.0; fi
        if [ -z "${three}" ]; then three=0.0; fi
        echo "$five 100.0 $three" | tr ' ' '\t' | awk '{print $1+$2+$3}'
    # if 3=6:
    elif [ "$line1" -ne "$line4" ] && [ "$line2" -ne "$line5" ] && [ "$line3" -eq "$line6" ]
    then
        five=$(grep -w "$line1" $COMPLEX_SCORES | grep -w "$line4" | cut -f5)
        n=$(grep -w "$line2" $COMPLEX_SCORES | grep -w "$line5" | cut -f5)
        if [ -z "${five}" ]; then five=0.0; fi
        if [ -z "${n}" ]; then n=0.0; fi
        echo "$five $n 100.0" | tr ' ' '\t' | awk '{print $1+$2+$3}'
    elif [ "$line1" -ne "$line4" ] && [ "$line2" -ne "$line5" ] && [ "$line3" -ne "$line6" ]
    then
        five=$(grep -w "$line1" $COMPLEX_SCORES | grep -w "$line4" | cut -f5)
        n=$(grep -w "$line2" $COMPLEX_SCORES | grep -w "$line5" | cut -f5)
        three=$(grep -w "$line3" $COMPLEX_SCORES | grep -w "$line6" | cut -f5)
        if [ -z "${five}" ]; then five=0.0; fi
        if [ -z "${n}" ]; then n=0.0; fi
        if [ -z "${three}" ]; then three=0.0; fi
        echo "$five $n $three" | tr ' ' '\t' | awk '{print $1+$2+$3}'
    fi

done <results/"$SAMPLE"_five_n_three_all_vs_all_scores/"$SAMPLE"_tmp > results/"$SAMPLE"_five_n_three_all_vs_all_scores/"$SAMPLE"_tmp2

# groupby based on first column:
cut -f1 results/"$SAMPLE"_five_n_three_all_vs_all_scores/"$SAMPLE"_tmp2 | sort | uniq | sort -k1n | while read line; do echo $line | tr '\n' '\t'; awk -v id="$line" '$1 == id {print $3}' results/"$SAMPLE"_five_n_three_all_vs_all_scores/"$SAMPLE"_tmp2 | tr '\n' '\t'; echo ""; done <"${1:-/dev/stdin}" | sed 's/,$//g' > results/"$SAMPLE"_five_n_three_all_vs_all_scores/"$SAMPLE"_tmp3
# groups by one column: column 1, and prints all column 3’s with it as TSV in a 2-column TSV (first column is the unique column 1 value and second column is the tab-separated matches from column 3)


# clean-up:
#rm -f results/"$SAMPLE"_five_n_three_all_vs_all_scores/"$SAMPLE"_tmp
#rm -f results/"$SAMPLE"_five_n_three_all_vs_all_scores/"$SAMPLE"_tmp2

