markers <- c("CD20", "CD3", "CD8", "CD45RO", "FOXP3")

phenotypes <- c("B", "Th", "cT", "Treg")

facs.mapping <- c(
	CD20="PE-A",
	CD3="BV421-A",
	CD8="PerCP-Cy5.5-A",
	CD45RO="APC-A",
	FOXP3="FITC-A"
)

inform.mapping.cytoagar <- c(CD3="Entire Cell CD3 (Opal 520) Mean (Normalized Counts, Total Weighting)",
	CD8="Entire Cell CD8 (Opal 690) Mean (Normalized Counts, Total Weighting)",
	CD20="Entire Cell CD20 (Opal 570) Mean (Normalized Counts, Total Weighting)",
	FOXP3="Entire Cell FOXP3 (Opal 540) Mean (Normalized Counts, Total Weighting)",
	CD45RO="Entire Cell CD45ro (Opal 620) Mean (Normalized Counts, Total Weighting)")

inform.mapping <- c(CD3="Entire Cell CD3 (Opal 520) Mean (Normalized Counts, Total Weighting)",
	CD8="Entire Cell CD8 (Opal 690) Mean (Normalized Counts, Total Weighting)",
	CD20="Entire Cell CD20 (Opal 570) Mean (Normalized Counts, Total Weighting)",
	FOXP3="Entire Cell FOXP3 (Opal 540) Mean (Normalized Counts, Total Weighting)",
	CD45RO="Entire Cell CD45RO (Opal 620) Mean (Normalized Counts, Total Weighting)")

agreement <- function(d, digits = 2, conf.int=FALSE) {
	library(irr)
	r <- icc(d, model="twoway", type="agreement", unit="single")
	cat(signif(r$value,digits), " (95% CI: ",signif( r$lbound, digits ), "-", signif( r$ubound, digits ), ")\n", sep="")
	if(conf.int) { 
		c(signif(r$value, digits), signif(r$lbound, digits), signif(r$ubound, digits))
	} else {
		signif(r$value, digits) 
	}
}

gate_dims <- list(
	CD20=c("CD3","CD20"),
	CD3=c("CD3","CD20"),
	CD8=c("CD3","CD8"),
	CD45RO=c("CD3","CD45RO"),
	FOXP3=c("CD3","FOXP3") 
)

apply_facs_gates <- function(ds, gates, ...) {
		r <- rep(TRUE, nrow(ds))
		for (wh in list(...)) {
			negated <- FALSE
			if (substr(wh,1,1) == "-") {
				negated <- TRUE
				wh <- substr(wh, 2, nchar(wh))
			}
			ingate <- as.logical(pointinpolygon( 
				list(x=ds[,facs.mapping [gate_dims[[wh]][1]]],
				     y=ds[,facs.mapping [gate_dims[[wh]][2]]]),
				as.list(gates[[wh]][[i]]) )
			)
			if (negated) {
				r <- r & !ingate
			} else {
				r <- r & ingate
			}
		}
		r
}

apply_immunet_thresholds <- function(ds, thresholds, ...) {
		r <- rep(TRUE, nrow(ds))
		for( wh in list(...) ){
			negated <- FALSE
			if (substr(wh,1,1) == "-") {
				negated <- TRUE
				wh <- substr(wh,2,nchar(wh))
			}
			ingate <- ds[,wh]>thresholds[wh]
			if (negated) {
				r <- r & !ingate
			} else {
				r <- r & ingate
			}
		}
		r
}

