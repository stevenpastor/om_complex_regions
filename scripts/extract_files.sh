#!/bin/bash

SAMPLE=11029C
ASSEMBLY=$SAMPLE'_Assembly_pipeline_results.tar.gz'

## For tar.gz assemblies
tar -xvf $ASSEMBLY output/contigs/exp_refineFinal1_sv/merged_smaps/exp_refineFinal1_merged.xmap \
output/contigs/exp_refineFinal1_sv/merged_smaps/exp_refineFinal1_merged_q.cmap \
output/contigs/exp_refineFinal1_sv/merged_smaps/exp_refineFinal1_merged_r.cmap \
output/all.bnx.gz \
output/molecules.tar.gz &&
# tar -xvf output/all.bnx.gz &&
tar -xvf output/molecules.tar.gz &&
mv output data/$SAMPLE'_output'