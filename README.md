# om_complex_regions

## Requirements

* samtools >=1.9
* python >=3.7
* numpy >=1.18.1
* pandas >=1.2.2
* scikit-learn >=0.23.2
* scipy >=1.5.3
* Bionano's RefAligner (included in repo's scripts) but see here: https://bionanogenomics.com/support/software-downloads/
* Linux environment: multi-threaded

## Operations

* edit config file
* run in order and from the root directory (not scripts but one level up):

```
sh scripts/complex_region_alignment.sh

sh scripts/separate_alignment.sh

sh scripts/fix_nested.sh

sh scripts/extract_query_labels_anchor.sh

python scripts/find_anchor_and_complex_alignment_matches.py

sh scripts/filter_maps_from_complex.sh

sh scripts/merge_into_one_cmap_self_align.sh

python scripts/filter_molecules_self_alignment.py

python scripts/proportional_comparisons.py

sh scripts/long_haps.sh

sh scripts/all_vs_all_scores_matrix.sh

sh scripts/final_molecules.sh
```
