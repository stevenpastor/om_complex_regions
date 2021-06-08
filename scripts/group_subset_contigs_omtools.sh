# Put all genomes' subset contigs into one map file (combined.cmap below)
# this naively uses results dir for all cmaps
# e.g., results/11029A_initial_genome_check/430_subset_103488556_104544470_final.cmap
# also requires RefAligner in root dir
find ./results/ -name "*_final.cmap" > input

rm -f combined.cmap; rm -f combined.idmap
./RefAligner -merge -if input -o combined
rm -f input

# NOTE: ref is same as query file
rm -f query*oma

#java -jar OMTools/OMTools.jar OMBlastMapper -refmapin combined.cmap -optmapin combined.cmap -optresout query1.oma -filtermode 1 -alignmentjoinmode 0 -thread 1 -exactmatch false -writeinfo false -writeunmap false -maxalignitem -1 -minconf 0 -fpp 15 -fnp 15 -ear 0.05 -meas 500 -minjoinscore 50 --allowequalrefquery false

#java -jar OMTools/OMTools.jar OMBlastMapper -refmapin combined.cmap -optmapin combined.cmap -optresout query2.oma -filtermode 1 -alignmentjoinmode 0 -thread 1 -exactmatch false -writeinfo false -writeunmap false -maxalignitem -1 -minconf 0 -fpp 10 -fnp 10 -ear 0.05 -meas 500 -minjoinscore 50 --allowequalrefquery false

#java -jar OMTools/OMTools.jar OMBlastMapper -refmapin combined.cmap -optmapin combined.cmap -optresout query3.oma -filtermode 1 -alignmentjoinmode 0 -thread 1 -exactmatch false -writeinfo false -writeunmap false -maxalignitem -1 -minconf 0 -fpp 5 -fnp 5 -ear 0.05 -meas 500 -minjoinscore 50 --allowequalrefquery false

#java -jar OMTools/OMTools.jar OMBlastMapper -refmapin combined.cmap -optmapin combined.cmap -optresout query4.oma -filtermode 1 -alignmentjoinmode 0 -thread 1 -exactmatch false -writeinfo false -writeunmap false -maxalignitem -1 -minconf 0 -fpp 50 -fnp 50 -ear 0.05 -meas 500 -minjoinscore 30 --allowequalrefquery false

#java -jar OMTools/OMTools.jar OMBlastMapper -refmapin combined.cmap -optmapin combined.cmap -optresout query5.oma -filtermode 1 -alignmentjoinmode 0 -thread 1 -exactmatch false -writeinfo false -writeunmap false -maxalignitem -1 -minconf 0 -fpp 50 -fnp 50 -ear 0.05 -meas 500 -minjoinscore 20 --allowequalrefquery false

java -jar OMTools/OMTools.jar OMBlastMapper -refmapin combined.cmap -optmapin combined.cmap -optresout query6.oma -filtermode 1 -alignmentjoinmode 0 -thread 1 -exactmatch false -writeinfo false -writeunmap false -maxalignitem -1 -minconf 0 -fpp 2 -fnp 2 -ear 0.05 -meas 500 -minjoinscore 50 --allowequalrefquery false

#echo """
#link query1.oma
#link query2.oma
#link query3.oma
#build allowrearrangement=true
#mergeproximity meas=500;ear=0.05
#merge query1.oma
#merge query2.oma
#merge query3.oma
#merge query4.oma
#merge query5.oma
#merge query6.oma
#mergeproximity meas=500;ear=0.05
#""" > maconfig

rm -f maconfig

echo """
link query6.oma
build allowrearrangement=true
mergeproximity meas=500;ear=0.05
merge query6.oma
mergeproximity meas=500;ear=0.05
""" > maconfig

# run the link, merge steps for multiple alignment:
rm -f matrixout_test; rm -f treeout_test; rm -f cb*_test
java -jar OMTools/OMTools.jar MultipleAlignment --optmapin combined.cmap --maconfig maconfig --maref combined.cmap --cblout cblout_test --cboout cboout_test

# Phylogenetic tree using the above:
java -jar OMTools/OMTools.jar UPGMATreeConstruction --matrixout matrixout_test --treeout treeout_test --optmapin combined.cmap --cblin cblout_test

## can now visualize the treeout_test in anything which can render a tree from Newick format (e.g., https://phylo.io/ or pick something in python or R)
## since this is only phylo tree, when want to see the blocks in maps, can then use OMView multiple alignment view to visualize the CNV/SD blocks




