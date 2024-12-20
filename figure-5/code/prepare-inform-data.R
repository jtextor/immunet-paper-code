rm(list=ls())

library(flowCore)

compensation.matrix <- "data/compensation-matrix.csv"
comp.mat <- read.csv(compensation.matrix, header=TRUE, check.names=FALSE, row.names=1)
rownames(comp.mat) <- gsub(" ::.*", "", rownames( comp.mat ))
colnames(comp.mat) <- gsub(" ::.*", "", colnames( comp.mat ))

read.facs <- function(i=1, color.cols=7:14) {
	print(i)
	x <- read.FCS(paste0("data/facs/tonsil0", i, ".fcs"))
	x <- compensate(x, comp.mat)
	xe <- exprs(x)
	tf <- transformList(from=colnames(x)[color.cols], tfun=logicleTransform())
	xt <- tf %on% x
	xte <- exprs(xt)[,color.cols]
	xte <- signif(xte,4)
}

read.cytoagar <- function(i=1, color.cols=4:10) {
	print(i)
	x <- read.FCS(paste0("data/inform/cytoagar0",i,".fcs"), truncate_max_range=FALSE)
	xte <- exprs(x)[,color.cols]
	xte
}

read.tonsil <- function(i=1, color.cols=4:10) {
	print(i)
	x <- read.FCS(paste0("data/inform/tonsil0",i,".fcs"), truncate_max_range=FALSE)
	xte <- exprs(x)[,color.cols]
	xte
}

saveRDS(lapply(1:6, read.cytoagar), file="data/cytoagars-inform.rds")
saveRDS(lapply(1:6, read.tonsil), file="data/tonsils-inform.rds")
