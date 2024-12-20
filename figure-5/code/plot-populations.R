rm(list=ls())

source("code/tools.R")
source("../settings.R")

d <- read.table("data/populations.txt", header=T)

logit <- function(x) log(x / (1 - x))
identity <- function(x) x

plt <- function(a="F", b="M_C", tf=log, xlab="% FCM", ylab="% mIHC ML",
                pops=c("CD3", "FOXP3", "CD8", "CD20"), conf.int=FALSE, xmin=-5,
                show.numbers=FALSE) {
	na <- list()
	nb <- list()
	for (m in markers) {
		na[[m]] <- paste0(m, "_", a)
		nb[[m]] <- paste0(m, "_", b)
	}

	den.x <- (d[[na[["CD3"]]]] + d[[na[["CD20"]]]] + d[[na[["CD8"]]]] + d[[na[["FOXP3"]]]])
	den.y <- (d[[nb[["CD3"]]]] + d[[nb[["CD20"]]]] + d[[na[["CD8"]]]] + d[[na[["FOXP3"]]]])

	clrmp <- c(CD3="Th", CD8="CTL", FOXP3="Treg", CD20="Bcell", CD45RO="Cd45ro")

	x <- tf(d[[na[[pops[1]]]]] / den.x)
	y <- tf(d[[nb[[pops[1]]]]] / den.y)

	xlim <- NULL
	ylim <- NULL

	if (length(pops) > 1) {
		xlim <- c(xmin, 0); ylim <- c(xmin, 0)
	} else {
		xlim <- c(min(c(x, y)) -0.2, max(c(x, y)) + 0.2)
		ylim <- xlim
	}

	plot(x, y, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, col=clrs[clrmp[pops[1]]],
	     xaxt="n", yaxt="n", asp=1 )

	for (axs in c(1, 2)) {
		tcks <- c(seq(.01, .09, by=.01), seq(.1, .9, by=.1))
		lbls <- rep(NA, length(tcks))
		lbls[1] <- 1
		lbls[length(lbls) / 2 + 1] <- 10
		lbls[length(lbls)] <- 100
		axis(axs, at=tf(tcks),
		labels=lbls,
		gap.axis=-100)
	}

	if (show.numbers) {
		text(x, y, as.character(seq_len(nrow(d))), pos=4)
	}

	i <- 2

	for (m in tail(pops, -1)) {
		xx <- tf(d[[na[[m]]]] / den.x)
		yy <- tf(d[[nb[[m]]]] / den.y)
		points(xx, yy, col=clrs[clrmp[m]])
		x <- c(x,xx)
		y <- c(y,yy)
		i <- i + 1
	}

	abline(0, 1)
	aggr <- agreement(cbind(x, y), conf.int=TRUE)

	if (conf.int) {
		legend("topleft", paste0("icc=", aggr[1], "\n(", aggr[2], "-", aggr[3], ")"), 
		       bty="n",inset=c(-.15,-.12), xpd=TRUE)
	} else {
		legend("topleft", paste0("icc=", aggr[1]), bty="n",
		       inset=c(-.15,-.12), xpd=TRUE)
	}
}

pdf_out( "plots/populations.pdf", width=6, height=6.4/4 + .2 )

## plot 1, with TREGs

par(mfrow=c(1, 4), mai=c(.42, .42, .25, .05), mgp=c(1.6, .6, 0), bty="l", pch=19)
par(cex=1)
plt(b="I_C", ylab="% mIHC AgarCyto", conf.int=TRUE)
mtext("baseline", 3, line=.4)
legend("bottomright", rev(c("B cell", "Th")), pch=16, col=rev(clrs[c("Bcell", "Th")]), bty="n", x.intersp=.4, inset=c(-.1, 0), xpd=TRUE)

plt(b="M_C", ylab="", conf.int=TRUE)
mtext("ImmuNet", 3, line=.4)
legend("bottomright", rev(c("CTL", "Treg")), pch=16, col=rev(clrs[c("CTL", "Treg")]), bty="n", x.intersp=.4, inset=c(-.1,0), xpd=TRUE)

plt(b="I_T", ylab="% mIHC tissue", conf.int=TRUE)
mtext("baseline", 3, line=.4)

plt(b="M_T", ylab="", conf.int=TRUE)
mtext("ImmuNet", 3, line=.4)

## plot 2, without TREGs

plt(b="I_C", ylab="% mIHC AgarCyto", pops=c("CD3","CD8","CD20"), xmin=-3.5, conf.int=TRUE)
mtext("baseline", 3, line=.4)
legend("bottomright", c("B cell", "Th", "CTL"), pch=16, col=clrs[c("Bcell","Th","CTL")], 
       bty="n", x.intersp=.4, inset=c(-.1, 0), xpd=TRUE)

plt(b="M_C", ylab="", conf.int=TRUE, xmin=-3.5, pops=c("CD3", "CD8", "CD20"))
mtext("ImmuNet", 3, line=.4)

plt(b="I_T", ylab="% mIHC tissue", pops=c("CD3", "CD8", "CD20"), xmin=-3.5, conf.int=TRUE)
mtext("baseline", 3, line=.4)

plt(b="M_T", ylab="", pops=c("CD3", "CD8", "CD20"), xmin=-3.5, conf.int=TRUE)
mtext("ImmuNet", 3, line=.4)

plt(a="I_C", b="M_C", xlab="% baseline", ylab="% ImmuNet", conf.int=TRUE)
mtext("mIHC AgarCyto", 3, line=.4)

plt(a="I_T", b="M_T", xlab="% baseline", ylab="% ImmuNet", conf.int=TRUE)
mtext("mIHC tissue", 3, line=.4)

plt(a="I_C", b="I_T", xlab="% mIHC AgarCyto", ylab="% mIHC Tissue", conf.int=TRUE)
mtext("baseline", 3, line=.4)
legend("bottomright", rev(c("B cell", "Th")), pch=16, col=rev(clrs[c("Bcell", "Th")]), bty="n", x.intersp=.4, inset=c(-.1, 0), xpd=TRUE)

plt(a="M_C", b="M_T", xlab="% mIHC AgarCyto", ylab="% mIHC Tissue", conf.int=TRUE)
mtext("ImmuNet", 3, line=.4)
legend("bottomright", rev(c("CTL", "Treg")), pch=16, col=rev(clrs[c("CTL", "Treg")]), bty="n", x.intersp=.4, inset=c(-.1,0), xpd=TRUE)

dev.off()

