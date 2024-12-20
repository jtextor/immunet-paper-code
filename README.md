# Code for ImmuNet paper 

This repository will be used to publish code that underlies analyses and figures shown in the forthcoming paper "ImmuNet: A Segmentation-Free Machine Learning Pipeline for Immune Landscape Phenotyping in Tumors by Muliplex Imaging" by Sultan et al (see preprint here: https://www.biorxiv.org/content/10.1101/2021.10.22.464548v2).

# Steps to reproduce the figures

The data needed to reproduce the analysis perfromed in the paper and to build the figures are located on [Zenodo achrive](https://zenodo.org/records/8084976); the data hierarchy matches the hierarchy of the code in this repository. Every figure folder should contain all the assets necessary to build it: code, data, external images, a Makefile and a latex file that defines the layout of a figure. The typical file hierarchy inside a figure folder is:

```
- code
- data
- images
- plots
Makefile
figure.tex
requirements.txt
run.sh
```

The `code` folder contains the source code that reproduces the analysis and generates auxiliary data files and plots. The `data` folder contains the data used in the analysis. External assets that are needed to complete the figure should be located in `images` folder. The plots that visualise the results of data analysis are saved in the `plots` folder, which should exist before the code is executed. `Makefile` specifies the order in which different parts of the analysis must be run with the `figure.pdf` file as a final target, built with `figure.tex` file. For convenience, `run.sh` script that fully automates the workfolw is provided. First, a conda environment is created based on the dependencies (R and Python) listed in the `requirements.txt` and then `make` is executed.

Executing `run.sh` is a recommended way to run the code. This requires `conda` (minimum version 24.9.2) and `latex` (including XeLaTeX engine) to be installed in your environment.
First, put the data related to a figure in the `data` folder and uncompressed the related `images.tar.gz` into the `images` folder. All `.tar.gz` archives should be uncompressed, whereas files compressed as `.gz`, should be left as they are. Then, run `bash run.sh` and the figure will be built automatically. Note, that creation of a conda environment and execution of the code might take a while for some figures. 

Instructions specific to some figures are given below.

## Figure 3 and 4

By default the cached prediction (`prediction.tar.gz`) of ImmuNet and the baseline algorithm InForm will be used for the analysis. However, an option to run inference with ImmuNet on tissue images is provided. For this you also need to uncompressed a corresponding image archive into the `data` folder:

- `tilecache.tar.gz` - for figure 3
- `rois_tilecache.tar.gz` - for figure 4

Then, you need to use a separate scripts, named `run_tf.sh` that create conda environments in which TensorFlow 1 is installed. Because of the confilict bethween the Python version that supports TensorFlow 1 and necessary R version, R packages will be installed in your default R.

## Figure 5

This figure uses Bioconductor software that is problematic to install with conda, therefore needed R packages have to be installed manually. First, install:

```
irr=0.84.1
BiocManager=1.30.25
polyclip=1.10-7
remotes=2.5.0
```

Then, `flowCore` should be installed with `BiocManager::install("flowCore")` and `tiltools` with `remotes::install_github("jtextor/tiltools")`. After that `run.sh` can be executed.




