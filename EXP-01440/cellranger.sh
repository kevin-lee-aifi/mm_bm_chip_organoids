#!/usr/bin/bash

home="/home/workspace"
config="${home}/mm_analysis/EXP-01440/EXP-01440_config.csv"
cr="${home}/cellranger-8.0.1/cellranger"
cr_outs="${home}/mm_analysis/EXP-01440/cr_outs"

${cr} multi --id='EXP-01440' \
    --csv=${config} \
    --output-dir=${cr_outs} \
    --localcores=60