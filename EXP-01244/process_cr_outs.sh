#!/usr/bin/bash

exp="EXP-01244"
proj="mm_bm_chip_organoids"

home="/home/workspace"
cr_outs="${home}/${proj}/${exp}/cr_outs/outs/per_sample_outs"
metrics="${home}/${proj}/${exp}/metrics_summaries"
webs="${home}/${proj}/${exp}/web_summary"

rm -rf "$metrics" "$webs"
mkdir -p "$metrics" "$webs"

for outs in $cr_outs/*; do
    id=$(echo $outs | cut -d'/' -f9)
    new_csv="${id}_metrics_summary.csv"
    new_html="${id}_web_summary.html"
    mv $outs/metrics* "$outs/${new_csv}"
    mv $outs/web* "$outs/${new_html}"
    cp "$outs/${new_csv}" "$metrics/${new_csv}"
    cp "$outs/${new_html}" "$webs/${new_html}"
done