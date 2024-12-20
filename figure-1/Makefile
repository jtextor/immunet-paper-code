.SECONDARY:
.DELETE_ON_ERROR:

all : figure.pdf

figure.pdf : figure.tex plots/expressions.pdf plots/simulated-expression-medium-density-full.pdf plots/simulated-expression-low-density.pdf 
	latexmk -use-make -xelatex $<

plots/simulated-expression-low-density.pdf : code/plot-expression-lo-vs-hi.R data/expression_data_low_density.csv.gz
	Rscript $<

plots/simulated-expression-medium-density-full.pdf : code/plot-expression-full-vs-membrane.R data/expression_data_full_vs_membrane.csv.gz
	Rscript $<

plots/expressions.pdf : code/plot-expressions.R data/expression_data.csv.gz
	Rscript $<

