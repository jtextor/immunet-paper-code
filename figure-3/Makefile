.SECONDARY:
.DELETE_ON_ERROR:

all : figure.pdf supporting-figure.pdf

supporting-figure.pdf : supporting-figure.tex images/accuracy.pdf 
	latexmk -use-make -xelatex $<

figure.pdf : figure.tex images/tuning.pdf images/accuracy-markers.pdf images/accuracy.pdf images/cartoon_hit_criterion.pdf
	latexmk -use-make -xelatex $<

images/accuracy.pdf : code/plot-accuracy.R performance_csvs_witness.txt
	Rscript $<

images/accuracy-markers.pdf : code/plot-accuracy-markers.R prediction_tsvs_witness.txt
	Rscript $<

images/tuning-pseudocolor-thresholds.pdf images/tuning.pdf : code/plot-tuning.R ../settings.R prediction_tsvs_witness.txt
	Rscript $<

performance_csvs_witness.txt : code/evaluation.py code/panels.py code/utils.py prediction_tsvs_witness.txt
	python3 $< && touch performance_csvs_witness.txt

prediction_tsvs_witness.txt : code/detection.py code/panels.py code/utils.py data/annotations_train.json.gz data/annotations_val.json.gz
	python3 $< && touch prediction_tsvs_witness.txt
