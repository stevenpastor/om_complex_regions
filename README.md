# om_complex_regions

## Requirements

* Samtools
* Python >=3.7
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
```
