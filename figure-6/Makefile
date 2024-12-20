all : figure.pdf supporting-figure.pdf

figure.pdf : figure.tex plots/nn.combined.pdf plots/stroma_tumor_counts.pdf infiltration_insets_witness.txt plots/low_density.png plots/high_density.png
	latexmk -use-make -xelatex $<

supporting-figure.pdf : supporting-figure.tex plots/dice_scores.pdf
	latexmk -use-make -xelatex $<

plots/dice_scores.pdf : code/dice_scores.r slide_dataframes_witness.txt data/inform_threshold_comp.csv data/inform_threshold_comp.stroma.csv
	Rscript code/dice_scores.r


infiltration_insets_witness.txt : code/plot_infiltration.r segmentation_images_witness.txt data/low_high.csv data/dims.csv data/offset.csv
	Rscript code/plot_infiltration.r
	touch infiltration_insets_witness.txt

plots/low_density.png plots/high_density.png : code/construct_overview_image.py data/low_high.csv segmentation_images_witness.txt
	python3 code/construct_overview_image.py

data/dims.csv : code/get_dims.py
	python3 code/get_dims.py

segmentation_images_witness.txt : code/make_seg_images.py data/low_high.csv
	python3 code/make_seg_images.py
	touch segmentation_images_witness.txt

data/low_high.csv : code/low_high.r data/infiltration.csv data/regions_per_distance.csv
	Rscript code/low_high.r

plots/stroma_tumor_counts.pdf : code/combine_infiltration.r data/infiltration.csv data/regions_per_distance.csv data/regions_per_distance.csv
	Rscript code/combine_infiltration.r


plots/nn.combined.pdf : code/combined_nn_plots.r data/neighbor_workspace.RData data/neighbor_per_cell_workspace.RData
	Rscript code/combined_nn_plots.r

data/neighbor_workspace.RData : code/neighbor.preprocess.r neighbor_dataframes_witness.txt data/regions_per_distance.csv
	Rscript code/neighbor.preprocess.r

data/neighbor_per_cell_workspace.RData : code/neighbor_per_cell.preprocess.r neighbor_dataframes_witness.txt data/infiltration.csv
	Rscript code/neighbor_per_cell.preprocess.r

data/regions_per_distance.csv : code/get_areas.py
	python3 code/get_areas.py

data/infiltration.csv data/distance_spread.csv : code/preprocess.r slide_dataframes_witness.txt
	Rscript code/preprocess.r

# idea of using witness files (sort of) taken from 
# https://www.gnu.org/software/automake/manual/html_node/Multiple-Outputs.html
slide_dataframes_witness.txt : code/make_slide_df.py
	python3 code/make_slide_df.py
	touch slide_dataframes_witness.txt

neighbor_dataframes_witness.txt : slide_dataframes_witness.txt code/process_df.py
	python3 code/process_df.py
	touch neighbor_dataframes_witness.txt


