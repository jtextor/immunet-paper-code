rm(list=ls())

library(tiltools)

source("../settings.R")

if (!exists("tonsils")) {
	tonsils <- readRDS("data/tonsil01.rds")
}

cd3 <- "BV421-A"
cd20 <- "PE-A"
cd8 <- "PerCP-Cy5.5-A"
foxp3 <- "FITC-A"
cd45ro <- "APC-A"


pdf_out("plots/biplots-facs.pdf", width=1, height=1)
par(mai=c(1-.7874,1-.7874,0,0),
    mgp=c(0,0,0), font.main=1, cex.main=1)

plt <- function(i, xpop=cd3, ypop=cd20, xlab="", ylab="", xth=3, yth=2.3) {
	ds <- tonsils[[i]]
	ds <- ds[!is.na(ds[,1]),]

	x <- ds[,xpop]
	y <- ds[,ypop]

	is.t <- x > xth & y <= yth
	is.b <- y > yth & x <= xth
	is.tt <- x > xth & y > yth

	cat(sum(is.t), "\t", sum(is.b), "\n", sep="")
	cat(max(x))

	densityplot2d(x, y,
	              xlim=c(-1,5.5), ylim=c(-2,5),
	              xlab=xlab, ylab=ylab, bty="l",
	              xaxt="n", yaxt="n", nbin=100, cex=.1, pch=16)
	abline(h=yth, col=2)
	abline(v=xth, col=2)
	legend("topleft", paste0(round(100 * mean(is.b)) / 1, "%"), bty="n", inset=c(-.12, -.09))
	legend("bottomright", paste0(round(100 * mean(is.t)) / 1, "%"), bty="n", inset=c(-.05, -.05))
}

plt(1, xlab="CD3", ylab="CD20")
dev.off()

