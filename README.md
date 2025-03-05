# Code for ImmuNet paper 

This repository will be used to publish code that underlies analyses and figures shown in the forthcoming paper "ImmuNet: A Segmentation-Free Machine Learning Pipeline for Immune Landscape Phenotyping in Tumors by Muliplex Imaging" by Sultan et al (see preprint here: https://www.biorxiv.org/content/10.1101/2021.10.22.464548v2).

# Requirements

For all figures except Figure 5, the following software is sufficient:

- `Conda`, minimum version 24.9.2
- `Latex`, including XeLaTeX engine

Then, a figure can be built by running `bash run.sh` from its root folder. This script automatically hangles creation of a virtual environment with all the depencencies installed, runs data analysis scripts in the correct order and eventually builds a figure with Latex. 

If you want to run inference on mIHC images to reproduce Figures 3 and 4 you need hardware that supports `TensorFlow 1.14.0` and `R` (minimum version 4.3.0) installed. By default the figures are built using the cached prediction. 

Figure 5:
- `R`, minimum version 4.3.0

# Steps to reproduce the figures

The data needed to reproduce the analysis performed in the paper and to build the figures are located on [Zenodo achrive](https://zenodo.org/records/8084976); the data hierarchy matches the hierarchy of the code in this repository. Every figure folder should contain all the assets necessary to build it: code, data, external images, a Makefile and a latex file that defines the layout of a figure. The typical file hierarchy inside a figure folder is:

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

The `code` folder contains the source code that reproduces the analysis and generates auxiliary data files and plots. The `data` folder contains the data required for the analysis. External assets that are needed to complete the figure should be located in `images` folder. The plots that visualise the results of data analysis are saved in the `plots` folder, which should exist before the code is executed. `Makefile` specifies the order in which different parts of the analysis must be run with the `figure.pdf` file as a final target, built with `figure.tex` file. `run.sh` script creates a conda environment based on the dependencies (R and Python) listed in the `requirements.txt` and then executes `make`.

Steps to reproduce a figure:
1. Place the data needed for a figure in the `data` folder. This include everything except `images.tar.gz`.
2. In the `data` folder, uncompressed all `.tar.gz` archives (but do not uncompressed `.gz` files).
3. Uncompress `images.tar.gz` into the `images` folder.
4. Run `bash run.sh`  

Note, that creation of a conda environment and execution of the code might take a while for some figures. 

Figure-specific instructions are given below.

## Figure 3 and 4

By default the cached prediction (`prediction.tar.gz`) of ImmuNet and the baseline algorithm InForm will be used for the analysis. However, an option to run inference with ImmuNet on mIHC images is provided. For this you also need to uncompressed a corresponding image archive into the `data` folder:

- `tilecache.tar.gz` - for figure 3
- `rois_tilecache.tar.gz` - for figure 4

Then, use a script, named `run_tf.sh` that creates a conda environment with TensorFlow 1.14. **Note**: because of the confilict between the Python version that supports TensorFlow 1 and necessary R version, R packages will be installed in your default R!

## Figure 5

This figure uses Bioconductor software that is problematic to install with conda, therefore needed R packages have to be installed manually. Steps to buld a figure:

1. Fill in the `data` and `images` folders as described above
2. Install R packages: `irr=0.84.1`, `BiocManager=1.30.25`, `polyclip=1.10-7`, `remotes=2.5.0`
3. Run `BiocManager::install("flowCore")`
4. Run `remotes::install_github("jtextor/tiltools")`
5. Run `make` in the root folder of the figure 

# Common issues

1. We use `Helvetica Neue` font family for text on figures. If it is not installed on your system, you can change fonts in the `settings.R` file located in the root folder of the repository and in the `\setmainfont` comman of each `figure.tex` file.
2. On some systems `conda` may fail to resolve the dependencies. Switching to [mamba](https://github.com/mamba-org/mamba) should fix the issue.  



