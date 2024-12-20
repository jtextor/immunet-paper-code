rm(list=ls())

source("../settings.R")

clrs <- clrs[c(7,2,3)]

`%notin%` <- Negate(`%in%`)

d <- read.csv("data/prediction-immunet-train.tsv", sep="\t" )

d$ann_type[d$ann_type %notin% c("T cell", "B cell")] <- "Other cell"
d$dist.mu <- d$dist / 2
d$dist.mu[d$dist.mu>=6] <- 6

h1 <- hist(d$dist.mu[d$ann_type=="T cell"], breaks=seq(-.5, 6.5), plot=FALSE)
h2 <- hist(d$dist.mu[d$ann_type=="B cell"], breaks=seq(-.5, 6.5), plot=FALSE)
h3 <- hist(d$dist.mu[d$ann_type=="Other cell"], breaks=seq(-.5, 6.5), plot=FALSE)

h <- 100 * rbind(h1$density, h2$density, h3$density)

cc <- c("0-1", "1-2", "2-3", "3-4", "4-5", "5-6", ">6")

colnames(h) <- cc

pdf_out("plots/tuning.pdf", width=2.7, height=2)

par(mar=c(3,3,1,0), mgp=c(1.8,.6,0))

xps <- barplot(h, beside=T, ylim=c(0,100), border=NA, col=clrs, xaxt="n",
               xlab=, ylab="% of annotated cells")
title(xlab=bquote("distance to nearest found cell (" * mu * "m)"), line=1.8)

mtext(side=1, cc, at=xps[2,], line=.6)
legend("topright", legend=c("T cells","B cells","other cells"),
	bty="n", pch=15, col=clrs, xpd=TRUE, inset=c(.1,0))

dev.off()

pdf_out("plots/tuning-pseudocolor-thresholds.pdf", width=2, height=2)

par(mar=c(3,3,0,0), mgp=c(1.8,.6,0))

ds <- d[d$dist.mu < 3.5 & d$ann_type %in% c("T cell", "B cell"),]

plot(ds$pCD3, ds$pCD20, xlim=c(-.2, 1.2), ylim=c(-.3, 1.4), pch=21,
     col=NA, bg=c(rgb(.9,0,.9,.1), rgb(.4,.1,.1,.1))[1+(ds$ann_type=="T cell")],
     cex=.5, xlab=bquote(psi * "CD3"),
     ylab=bquote(psi * "CD20"), bty="l")

legend("topright", pch=19, col=clrs[c(1,2)], legend=c("T cell", "B cell"), bty="n")

abline( v = 0.4 )
abline( h = 0.4 )

sdig <- function(x) paste0(round(100 * x, 1), "%")

legend("bottomright", legend=
	sdig(mean(ds$pCD3[ds$ann_type == "T cell"] > 0.4 & ds$pCD20[ds$ann_type == "T cell"] < 0.4)),
	text.col=clrs[1], bty="n", inset=c(0,-.03))

legend("topleft", legend=
	sdig(mean( ds$pCD3[ds$ann_type == "B cell"] < 0.4 & ds$pCD20[ds$ann_type == "B cell"] > 0.4 )),
	text.col=clrs[2],
	bty="n", inset=c(-.12,-.06))

dev.off()
