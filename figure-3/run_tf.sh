#!/usr/bin/env bash

conda env create -f environment_tf.yml

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

conda activate immunet_fig3_tf

Rscript setup_tf.R

make
