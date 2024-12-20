.SECONDARY:
.DELETE_ON_ERROR:

all : figure.pdf

figure.pdf : figure.tex images_witness.txt
	latexmk -use-make -xelatex $<


images_witness.txt : code/label_maps_vis.py data/annotations.json.gz images/composite.png 
	python3 $< && touch images_witness.txt
