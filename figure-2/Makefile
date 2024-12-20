.SECONDARY:
.DELETE_ON_ERROR:

all : figure.pdf

figure.pdf : figure.tex 
	latexmk -use-make -xelatex $<

