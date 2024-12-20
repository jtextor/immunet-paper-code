#!/usr/bin/env bash

conda create -n immunet_fig6 python=3.8 r-base=4.3.1 r-essentials

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

conda activate immunet_fig6
conda install --file requirements.txt

make
