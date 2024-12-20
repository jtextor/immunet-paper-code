#!/usr/bin/env bash

conda env create -f environment.yml

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh

conda activate immunet_supp_f4

make
