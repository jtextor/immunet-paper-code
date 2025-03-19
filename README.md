# Code for ImmuNet paper 

This repository contains the code that underlies analyses and figures shown in the paper "ImmuNet: A Segmentation-Free Machine Learning Pipeline for Immune Landscape Phenotyping in Tumors by Muliplex Imaging" by Sultan et al ([https://doi.org/10.1093/biomethods/bpae094](https://doi.org/10.1093/biomethods/bpae094)).

# Requirements

For all figures except Figure 5, the following software is sufficient:

- `Conda`, minimum version 24.9.2
- `Latex`, including XeLaTeX engine

Then, a figure can be built by running `bash run.sh` from its root folder. This script automatically creates a virtual environment with all dependencies installed, runs data analysis scripts in the correct order, and eventually builds a figure with Latex. 

By default, the code that evaluates the model performance shown in Figures 3 and 4 uses the cached prediction for time efficiency. We also provide an option to directly detect immune cells in in mIHC images with the final model used in the paper. This requires hardware supporting `TensorFlow 1.14.0` and `R` (minimum version 4.3.0).

Figure 5:
- `R`, minimum version 4.3.0

# Steps to reproduce the figures

The data needed to reproduce the analysis performed in the paper and to build the figures are located on [Zenodo archive](https://zenodo.org/records/15046015); the folder hierarchy in the data archive matches the folder hierarchy in this repository. Before building a figure, you need to make sure that its folder contains all the necessary assets: code, data, external images, a Makefile, and a Latex file that defines the figure layout. The typical file hierarchy inside a figure folder is:

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

The `code` folder contains the source code that reproduces the analysis and generates auxiliary data files and plots. The `data` folder contains the data required for the analysis. External assets that are needed to complete the figure should be located in the `images` folder. The plots that visualise the results of the data analysis are saved in the `plots` folder, which is created automatically. The `Makefile` specifies the order in which the analysis scripts must be run, with the `figure.pdf` file as the final target, built with the `figure.tex` file. The `run.sh` script creates a conda environment based on the dependencies (R and Python) listed in the `requirements.txt` and then executes `make`.

Steps to reproduce a figure:
1. Place the data needed for a figure in the `data` folder. This includes everything except `images.tar.gz`,
2. In the `data` folder, uncompress all the `.tar.gz` archives (but do not uncompress the `.gz` files),
3. Extract `images.tar.gz` into the `images` folder,
4. Run `bash run.sh`.  

Note, that creation of a conda environment and execution of the code might take a while for some figures. 

Figure specific instructions are given below.

## Figures 3 and 4

By default, the cached prediction (`prediction.tar.gz`) from ImmuNet and the baseline algorithm InForm will be used to reproduce the analysis. However, an option to run inference on mIHC images with ImmuNet is provided. To do this, you need to extract the corresponding image archive into the `data` folder:

- `tilecache.tar.gz` - for Figure 3
- `rois_tilecache.tar.gz` - for Figure 4

Then, use the `run_tf.sh` script which creates a conda environment with TensorFlow 1.14. **Note**: due to of the conflict between the Python version that supports TensorFlow 1 and the required R version, R packages will be installed in your default R!

## Figure 5

This figure uses Bioconductor software that is difficult to install with conda, therefore the required R packages have to be installed manually. Steps to build a figure:

1. Fill in the `data` and `images` folders as described above,
2. Install R packages: `irr=0.84.1`, `BiocManager=1.30.25`, `polyclip=1.10-7`, `remotes=2.5.0`,
3. Run `BiocManager::install("flowCore")`,
4. Run `remotes::install_github("jtextor/tiltools")`,
5. Run `make` in the root folder of the figure. 

# Common issues

1. We use the `Helvetica Neue` font family for text on figures. If it is not installed on your system, you can change the fonts in the `settings.R` file located in the root folder of the repository and in the `\setmainfont` command of each `figure.tex` file.
2. On some systems `conda` may fail to resolve the dependencies. Switching to [mamba](https://github.com/mamba-org/mamba) should fix the issue.  



