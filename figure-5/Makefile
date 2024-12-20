.SECONDARY:
.DELETE_ON_ERROR:

all : figure.pdf supporting-figure.pdf

figure.pdf : figure.tex plots/populations.pdf plots/biplots-facs.pdf plots/biplots.pdf
	latexmk -use-make -xelatex $<

supporting-figure.pdf : supporting-figure.tex plots/populations.pdf plots/foxp3-gating.pdf
	latexmk -use-make -xelatex $<

plots/foxp3-gating.pdf : code/plot-foxp3-gating.R data/tonsil01.rds
	Rscript $<

plots/populations.pdf : code/plot-populations.R data/populations.txt
	Rscript $<

data/populations.txt : code/prepare-populations.R data/cytoagars-inform.rds
	Rscript $<

plots/biplots.pdf : code/plot-biplots.R data/biplots.rds
	Rscript $<

data/biplots.rds : code/prepare-biplots.R external-data/inform/cytoagar01.fcs
	Rscript $<

plots/biplots-facs.pdf : code/plot-biplots-facs.R data/tonsil01.rds
	Rscript $<

data/tonsil01.rds : code/prepare-biplots-facs.R
	Rscript $<

data/cytoagars-inform.rds : code/prepare-inform-data.R external-data/inform/cytoagar01.fcs
	Rscript $<
