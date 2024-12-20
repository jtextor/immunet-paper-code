rm(list=ls())

library(flowCore)
suppressMessages(library(polyclip))

source("code/tools.R")

gates.facs <- readRDS("data/gates.rds")
data.facs <- readRDS("data/tonsils-facs.rds")

gates.inform.cytoagar <- readRDS("data/gates-inform-cytoagar.rds")
data.inform.cytoagar <- readRDS( "data/cytoagars-inform.rds" )

gates.inform.tissue <- readRDS("data/gates-inform-tissue.rds")
data.inform.tissue <- readRDS( "data/tonsils-inform.rds" )

inform.mapping.tissue <- inform.mapping

manual.gating <- list(list(n="CD20", a="CD3", b="CD20"), 
                      list(n="CD3", a="CD3", b="CD20"),
                      list(n="CD8", a="CD3", b="CD8"),
                      list(n="CD45RO", a="CD3", b="CD45RO"),
                      list(n="FOXP3", a="CD3", b="FOXP3"))

sink("data/populations.txt")
options(digits=3)
for (m in markers) {
	cat(m, "_F ", sep="")
}

for (suff in c("I","M")) {
	for(wh in c("T","C")) {
		for(m in markers) {
			cat(m, "_", suff, "_", wh, " ", sep="")
		}
	}
}

cat("\n")
for(i in 1:6) {
	cat( 
		mean(apply_facs_gates(data.facs[[i]], gates.facs, "CD20","-CD3","-CD8","-FOXP3")), " ",
		mean(apply_facs_gates(data.facs[[i]], gates.facs, "CD3","-CD8","-FOXP3","-CD20")), " ",
		mean(apply_facs_gates(data.facs[[i]], gates.facs, "CD3","CD8","-FOXP3","-CD20")), " ",
		mean(apply_facs_gates(data.facs[[i]], gates.facs, "CD3","CD45RO","-CD8","-CD20")), " ",
		mean(apply_facs_gates(data.facs[[i]], gates.facs, "CD3","FOXP3","-CD8","-CD20")), " "
	)

	for(wh in c("tissue","cytoagar")) {
		xx3.ph <- get(paste0("data.inform.",wh))[[i]][,1]
		xx3.ph <- xx3.ph[is.finite(xx3.ph)]

		cat( mean(xx3.ph < 1.5)," ") # B cells
		cat( mean(xx3.ph > 1.5 & xx3.ph < 3.5 )," ") # Th  (CD45RO+ / -)
		cat( mean(xx3.ph > 3.5 & xx3.ph < 5.5 )," ") # CTL  (CD45RO+ / -)
		cat( mean(xx3.ph > 2.5 & xx3.ph < 3.5 )," ") # Th  (CD45RO+)
		cat( mean(xx3.ph > 5.5 & xx3.ph < 7.5 )," ") # Treg  (CD45RO+ / -)
	}

	marker.thresholds <- c(CD3=0.4, CD20=0.4, CD8=0.4, FOXP3=0.4, CD45RO=0.25)

	for(wh in c("tonsil", "cytoagar")) {
		xx <- exprs(read.FCS(paste0("data/immunet/", wh, "0", i, ".fcs")))

		cat( 
			mean(apply_immunet_thresholds(xx, marker.thresholds, "CD20", "-CD3", "-CD8", "-FOXP3")), " ",
			mean(apply_immunet_thresholds(xx, marker.thresholds, "CD3", "-CD8", "-FOXP3", "-CD20")), " ",
			mean(apply_immunet_thresholds(xx, marker.thresholds, "CD3", "CD8", "-FOXP3", "-CD20"))," ",
			mean(apply_immunet_thresholds(xx, marker.thresholds, "CD3", "CD45RO", "-CD8", "-CD20"))," ",
			mean(apply_immunet_thresholds(xx, marker.thresholds, "CD3", "FOXP3", "-CD8", "-CD20"))," "
		)
	}

	cat("\n")
}
sink()

