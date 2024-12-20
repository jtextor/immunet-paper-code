rm(list=ls())

library(flowCore)

compensation.matrix <- "external-data/compensation-matrix.csv"
comp.mat <- read.csv(compensation.matrix, header=TRUE, check.names=FALSE, row.names=1)
rownames(comp.mat) <- gsub(" ::.*", "", rownames( comp.mat ))
colnames(comp.mat) <- gsub(" ::.*", "", colnames( comp.mat ))

color.cols <- 7:14

read.tonsil <- function(i=1) {
	print(i)
	x <- read.FCS(paste0("external-data/facs/tonsil0", i, ".fcs"))
	x <- compensate(x, comp.mat)
	xe <- exprs(x)
	tf <- transformList(from=colnames(x)[color.cols], tfun=logicleTransform())
	xt <- tf %on% x
	xte <- exprs(xt)[,color.cols]
	xte <- signif(xte,4)
}

tonsils <- list(read.tonsil(1))

saveRDS( tonsils, file="data/tonsil01.rds" )
