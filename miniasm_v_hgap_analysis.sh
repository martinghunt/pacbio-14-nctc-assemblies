#!/usr/bin/env bash
set -e

# Make ACT-style figures
for x in `cat nctc_ids.txt`
do
#  ./act_cartoon.py --min_id 80 $x/miniasm.fa $x/hgap.fa $x/ref.fa $x/act_cartoon.miniasm_ref_hgap
  convert -resize 1200x360 -density 2000 $x/act_cartoon.miniasm_ref_hgap.svg $x/act_cartoon.miniasm_ref_hgap.png
done


# Gather assembly stats and make a few plots
# assembly-stats available here: https://github.com/sanger-pathogens/assembly-stats
assembly-stats -t */*.fa | awk 'BEGIN{OFS="\t"; getline; $1="sample\tassembler"; print $0} {split($1,a,"/"); $1=a[1]"\t"a[2]; print $0}' | sed 's/\.fa//' > miniasm_v_hgap.assembly_stats.tsv

echo -e "Sample\tAssembler\tLength_ratio" > miniasm_v_hgap.hit_length_compare.tsv
for x in `cat nctc_ids.txt`; do awk -vsample=$x '$1~/^[0-9]+/ {print sample"\tminiasm\t"100*$6/$5}' $x/act_cartoon.miniasm_ref_hgap.nucmer_before_v_ref.coords; done >> miniasm_v_hgap.hit_length_compare.tsv
for x in `cat nctc_ids.txt`; do awk -vsample=$x '$1~/^[0-9]+/ {print sample"\tHGAP\t"100*$6/$5}' $x/act_cartoon.miniasm_ref_hgap.nucmer_after_v_ref.coords; done >> miniasm_v_hgap.hit_length_compare.tsv

./miniasm_v_hgap_analysis_plots.R
rm -f Rplots.pdf
