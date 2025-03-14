#!/usr/bin/bash

hdir="/home/workspace/"
dir1="/home/workspace/mm_bm_chip_organoids/EXP-01244"
dir2="/home/workspace/mm_bm_chip_organoids/EXP-01244/cr_outs"

output="${dir1}/EXP-01244_cr_outs.tar.gz"

mkdir -p ${dir2}

gcloud auth login

gcloud storage cp gs://experimental-pipeline-data/EXP-01244/EXP-01244_cr_outs.tar.gz ${output}

tar -xzvf ${output} -C ${dir2}

rm -r ${output}

outs1="${hdir}/mm_bm_chip_organoids/EXP-01244/cr_outs/home"
outs2="${hdir}/mm_bm_chip_organoids/EXP-01244/cr_outs/home/jupyter/TissDiss/EXP-01244/EXP-01244_cr_outs"

mv ${outs2}/* ${dir2}

rm -fr ${outs1}