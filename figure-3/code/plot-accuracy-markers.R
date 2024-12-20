rm(list=ls())

source("../settings.R")

library(polycor)

`%notin%` <- Negate(`%in%`)


d <- read.csv( "data/prediction-immunet-train.tsv", sep="\t")
d <- d[d$tissue != "cytoagar",]
d$tissue <- factor(d$tissue)

pdf_out( "plots/accuracy-markers.pdf", width=3, height=2 )

par( mar=c(3,1.6,1,.2), mgp=c(1.8,.6,0), cex.sub=1, cex.lab=1 )
layout( matrix( c(1,2,3), nrow=1 ), widths=c(1,1,1.5) )
par( cex.sub=1, cex.lab=1, cex=1 )
dt <- d[d$ann_type == "T cell",]

print( head(dt) )

fp <- function( mrk="CD8", col=5 ){

 	mru <- toupper( mrk )
	dtp <- dt[dt[,mru]>3,]
	dtn <- dt[dt[,mru]<3,]

	pseudo_mrk <- paste0("p", mru)

	xdp <- density( dtp[,pseudo_mrk], bw=0.02 )
	mxp <- max(xdp$y)
	xdp$y <- xdp$y / mxp

	xdn <- density( dtn[,pseudo_mrk], bw=0.02 )
	mxn <- max(xdn$y)
	xdn$y <- xdn$y / mxn

	plot( NA, type='l', xlim=c(0,2.25), xaxt="n", yaxt="n", ylab="", xlab="", bty="n", ylim=c(0,1) )

	axis( 2, at=c(0,1) )
	mtext( bquote( psi * .(mru) ), 2, line=0.2 )

	polygon( .5+xdn$y/2, xdn$x, col=clrs[col], border=NA, xpd=TRUE )
	polygon( .5-xdn$y/2, xdn$x, col=clrs[col], border=NA, xpd=TRUE )

	polygon( 1.75+xdp$y/2, xdp$x, col=clrs[col], border=NA, xpd=TRUE )
	polygon( 1.75-xdp$y/2, xdp$x, col=clrs[col], border=NA, xpd=TRUE )


	abline( h=0.4 )
	text( 0, 0.3, paste0( round( 100*mean( dtn[,pseudo_mrk] < 0.4 ), 1 ), "%" ), adj=c(0,0) )
	text( 2.25, 0.5, paste0( round( 100*mean( dtp[,pseudo_mrk] > 0.4 ), 1 ), "%" ), adj=c(.9,1), xpd=TRUE )

	mtext( mru, 1, line=1 )

	text( 1.25*c(0.5,1.5), -.1, labels=c("-","+"), xpd=TRUE )
}

fp("CD8",col=4)
fp("FOXP3")

mru <- "CD45RO"

plot( NA, type='l', xlim=c(0,6), xaxt="n", yaxt="n", ylab="", xlab="", bty="n", ylim=c(0,1) )

for( i in 1:5 ){
	dtp <- dt[dt[,mru]==i,]
	pseudo_mrk <- paste0("p", mru)
	xdp <- density( dtp[,pseudo_mrk], bw=0.02 )
	mxp <- max(xdp$y)
	xdp$y <- xdp$y / mxp
	polygon( xdp$y/2+i*1.25-.5, xdp$x, col=clrs[6], border=NA, xpd=TRUE )
	polygon( xdp$y/-2+i*1.25-.5, xdp$x, col=clrs[6], border=NA, xpd=TRUE )

}

axis( 2, at=c(0,1) )

text( 0+1.25*(0:4)+1.25/2, -.1, labels=c("--","-","o","+","++"), xpd=TRUE )

mtext( bquote( psi * .(mru) ), 2, line=0.2 )
mtext( mru, 1, line=1 )

pscor <- signif(polyserial(dt$pCD45RO, dt$CD45RO), 2)
legend("topleft", legend=bquote(rho == .(pscor)), bty="n", inset=c(-.1, -.02))

dev.off()
