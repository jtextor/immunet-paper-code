
clrs <- c(Th=rgb(.6,0,0),Bcell=rgb(.9,0,.9),other=rgb(.3,.3,.3),CTL=rgb(0,.8,.8,.7),Treg=rgb(0,.9,0,.7),Cd45ro=rgb(.9,.9,0,.7), 
	Tcell=rgb(.4,.1,.1,1) )

quartzFonts(sans = quartzFont(c("Helvetica Neue Light", "Helvetica Neue", "Helvetica Light","Helvetica")))

if( capabilities()['aqua'] ) {
  default_font <- "Helvetica Neue Light"
} else {
  default_font <- "Helvetica"
}

default_pointsize <- 9

pdf_out <- function( file, pointsize=default_pointsize, 
	dpi=NA_real_, family=default_font, ... ){
	if( capabilities()['aqua'] ){
		# preferred device type with configurable fonts
		quartz( type="pdf", file=file, pointsize=pointsize, dpi=dpi, ...)
		par( family="sans" )

	} else {
		# fallback device type with default Helvetica font
		pdf( file, pointsize=pointsize, ...)
	}
}
