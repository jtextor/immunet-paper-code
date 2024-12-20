rm(list=ls())

library(flowCore)

color.cols <- c("Entire Cell CD3 (Opal 520) Mean (Normalized Counts, Total Weighting)",
	"Entire Cell CD20 (Opal 570) Mean (Normalized Counts, Total Weighting)")

tonsils <- list(
	exprs(read.FCS(paste0("external-data/inform/cytoagar01.fcs"), truncate_max_range=FALSE))[,color.cols],
	exprs(read.FCS(paste0("external-data/inform/tonsil01.fcs"), truncate_max_range=FALSE))[,color.cols],
	exprs(read.FCS("external-data/immunet/cytoagar01.fcs"))[,c("CD3", "CD20")],
	exprs(read.FCS("external-data/immunet/tonsil01.fcs"))[,c("CD3", "CD20")]
)

saveRDS(tonsils, file="data/biplots.rds")


