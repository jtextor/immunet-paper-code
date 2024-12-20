rm(list=ls())

source("../settings.R")

plt <- function(file="data/immunet_perf_val_subtypes.csv", ymax=25) {
  par(mar=c(3, 3, 1, 0), mgp=c(1.8, .6, 0))
  
  ci_lower <- function(x) {
    100 * binom.test(x[1], x[2])$conf.int[1]
  }
  
  ci_upper <- function(x) {
    100*binom.test(x[1], x[2])$conf.int[2]
  }

  rr0 <- read.csv(file, row.names=1)
  rr0 <- rr0[rr0$tissue != "cytoagars",]
  
  rr0$err_rate <- 100*rr0$err_rate
  rr0$ci_lower <- apply(rr0[,c("errors", "cases")], 1, ci_lower)
  rr0$ci_upper <- apply(rr0[,c("errors", "cases")], 1, ci_upper)
  
  
  rr.c <- reshape(rr0[,c("ann_subtype", "tissue", "err_rate")], idvar="ann_subtype", timevar="tissue", direction="wide")
  rownames(rr.c) <- rr.c[,1]
  rr.c <- rr.c[,-1]
  
  rr.c_hi <- reshape(rr0[,c("ann_subtype", "tissue", "ci_upper")], idvar="ann_subtype", timevar="tissue", direction="wide")
  rownames(rr.c_hi) <- rr.c_hi[,1]
  rr.c_hi <- rr.c_hi[,-1]
  
  rr.c_lo <- reshape( rr0[,c("ann_subtype","tissue","ci_lower")], idvar="ann_subtype", timevar="tissue", direction="wide")
  rownames(rr.c_lo) <- rr.c_lo[,1]
  rr.c_lo <- rr.c_lo[,-1]
  
  colnames(rr.c) <- gsub("err_rate.", "", colnames( rr.c ))
  
  print(rr.c)
  
  phh <- c("B cell", "Helper T cell", "Cytotoxic T cell", "Regulatory T cell", "Other cell")
  rr.c <- rr.c[phh,]
  rr.c_hi <- rr.c_hi[phh,]
  rr.c_lo <- rr.c_lo[phh,]
  
  ccol <- clrs[c(2,1,4,5,3)]
  
  xps <- barplot(as.matrix(rr.c),
                 beside=T, ylim=c(0,ymax), space=c(0.08,1),
                 ylab="errors (%)", border=ccol, xaxt="n",
                 col=ccol)
  
  axis(1, at=colMeans(xps), labels=colnames(rr.c), gap.axis=-1000, col=NA, col.ticks=NA)
  segments(xps, as.matrix(rr.c_lo), xps, as.matrix(rr.c_hi))
  
  legend("topleft", pch=15, col=ccol,
         legend = phh,
         inset=c(0,-.1),
         bty="n", xpd=TRUE)
  
  	print(table(rr.c < 10))
  	print(table(rr.c < 5))
}

pdf_out("plots/accuracy.pdf", width=3.5, height=2)

plt()
plt("data/immunet_perf_train_subtypes.csv", 30)
plt("data/inform_perf_val_subtypes.csv", 100)
plt("data/inform_perf_train_subtypes.csv", 100)
dev.off()
