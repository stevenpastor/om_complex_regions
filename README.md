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
scripts/complex_region_alignment.sh

scripts/separate_alignment.sh

scripts/fix_nested.sh

scripts/extract_query_labels_anchor.sh
```
