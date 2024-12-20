rm(list=ls())

library( tiltools )

source("../settings.R")

pdf_out("plots/biplots.pdf", width=2.12, height=2)
par(mar=c(1.2,1.3,0,0), mfrow=c(2,2), omi=c(0, 0, 0, 0.12),
	mgp=c(0.2, 0,0), font.main=1, cex.main=1, xpd=NA)

par(cex=1)

plt <- function(i, xpop=1, ypop=2, xlim=c(0, 100), ylim=c(0, 100),
		xlab="", ylab="", xth=3, yth=2.3, logtransform=TRUE) {
  
	ds <- readRDS("data/biplots.rds")[[i]]
	print(colnames(ds))

	ds <- ds[!is.na(ds[,1]),]

	x <- ds[,xpop];
	y <- ds[,ypop];

	if (logtransform) {
		x <- log(1 + x)
		y <- log(1 + y)
	}

	x[x > xlim[2]] <- xlim[2] * 0.99
 	y[y > ylim[2]] <- ylim[2] * 0.99

	is.t <- x > xth & y <= yth
	is.b <- y > yth & x <= xth
	is.tt <- x > xth & y > yth

	cat(sum(is.t), "\t", sum(is.b), "\n", sep="")

	cat(max(x))

	densityplot2d(x, y,
		xlim=xlim, ylim=ylim,
		xlab=xlab, ylab=ylab, bty="l",
		tf.fun=log,
		xaxt="n", yaxt="n", nbin=100, cex=0, pch=16 )
}

plt(1, xlab="CD3", ylab="CD20", xlim=c(0,4), ylim=c(0,6.5))
mtext("AgarCyto", 3, line=-.9)
plt(2, xlab="CD3", ylab="CD20", xlim=c(0,4), ylim=c(0,6.5))
mtext("tissue", 3, line=-.9)
text("inForm", x=4.5, y=3.25, srt=270)
plt(3, xlab=expression( psi*"CD3" ), ylab=expression( psi*"CD20" ),
	xlim=c(-.2,1.2), ylim=c(-.25,1.2), logtransform=FALSE)

plt(4, xlab=expression( psi*"CD3" ), ylab=expression( psi*"CD20" ),
	xlim=c(-.2,1.2), ylim=c(-.25,1.2), logtransform=FALSE)
text("ImmuNet", x=1.375, y=0.7, srt=270)

dev.off()


