all : figure.pdf plots/f_score.table.pdf supporting-figure.pdf

figure.pdf : figure.tex plots/counts_density_histogram.pdf plots/f_score_r_value.pdf plots/density_f_score_r_value.pdf
	latexmk -use-make -xelatex $<

plots/f_score.table.pdf : code/table.r data/rois_info.csv data/detection_stat.csv data/cd45ro_stat.csv
	Rscript code/table.r

supporting-figure.pdf plots/density_f_score_r_value.pdf : code/plot_density_accuracy.r data/rois_info.csv data/detection_stat.csv data/cd45ro_stat.csv  ../settings.R
	Rscript code/plot_density_accuracy.r

plots/f_score_r_value.pdf : code/fscores_per_tissue.r data/rois_info.csv data/detection_stat.csv data/cd45ro_stat.csv ../settings.R
	Rscript code/fscores_per_tissue.r

data/rois_info.csv data/detection_stat.csv data/cd45ro_stat.csv : code/analyse.py
	python3 code/analyse.py

plots/counts_density_histogram.pdf : code/counts.r data/counts.csv ../settings.R
	Rscript code/counts.r

data/counts.csv : code/counts.py
	python3 code/counts.py
