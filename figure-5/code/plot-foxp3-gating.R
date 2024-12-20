rm(list=ls())

library(polyclip)
library(tiltools) 

source("../settings.R")

pdf_out("plots/foxp3-gating.pdf", width=6, height=1 )
par(mar=c(1.2,1.2,0,0), mfrow=c(1,6),
    mgp=c(0.2,0,0), font.main=1, cex.main=1)

par(cex=1)


if (!exists("tonsils")) {
	tonsils <- readRDS("data/tonsils-facs.rds")
}

color_names <- c(
  CD3="BV421-A",
	CD20="PE-A",
	CD8="PerCP-Cy5.5-A",
	FOXP3="FITC-A",
	CD45RO="APC-A"
)

gate_dims <- list(
	CD20=c("CD3", "CD20"),
	CD3=c("CD3", "CD20"),
	CD8=c("CD3", "CD8"),
	CD45RO=c("CD3", "CD45RO"),
	FOXP3=c("CD3", "FOXP3")
)

for (i in 1:6) {

	ds <- tonsils[[i]]
	ds <- ds[!is.na(ds[,1]),]

	## plot all cells
	gates.facs <- readRDS("data/gates.rds")
	gg <- gates.facs[["FOXP3"]][[i]]

	densityplot2d(ds[,color_names["CD3"]], ds[,color_names["FOXP3"]], 
		xlim=c(-1,5.5), ylim=c(-.5,5),
		xlab="CD3", ylab="FOXP3", bty="l",
		xaxt="n", yaxt="n", nbin=100, cex=.1, pch=16 )
	
	polygon(gg, border=2)
}

dev.off()

