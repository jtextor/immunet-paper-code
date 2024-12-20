rm(list=ls())

library(dplyr)
library(readr)
library(ggplot2)
library(RColorBrewer)

source("../settings.R")

df <- read_csv("data/expression_data_low_density.csv.gz")

myColor <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
myColor <- c(NA,NA,NA,NA, "white","grey","darkblue","lightblue", "green", "yellow","orange", "red")
myColor <- c("grey","darkblue","lightblue", "green", "yellow","orange", "red")
myColor_scale_fill <- scale_fill_gradientn(colours = myColor, trans="sqrt")

p <- df %>%
ggplot(aes(cd3, cd20)) +
  stat_binhex(bins=75) +
  xlab("CD3") +
  ylab("CD20") +
  myColor_scale_fill +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),  
        strip.text.x = element_blank(), 
        strip.background = element_blank(),
        text = element_text(size=8, family="sans"),
        legend.position = "none",
        axis.text = element_blank(), 
        axis.ticks = element_blank()) +
  scale_x_log10(limits=c(0.002,NA)) + 
  scale_y_log10(limits=c(0.002, NA))

pdf_out("plots/simulated-expression-low-density.pdf", width=1.5, height=1.5)
print(p)
dev.off()

df <- read_csv("data/expression_data_high_density.csv.gz")
 
p <- df %>%
ggplot( aes(cd3, cd20)) +
  stat_binhex(bins=75) +
  xlab("CD3") +
  ylab("CD20") +
  myColor_scale_fill +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),  
        strip.text.x = element_blank(), 
        strip.background = element_blank(),
        text = element_text(size=8, family="sans"),
        legend.position = "none",
        axis.text = element_blank(), 
        axis.ticks = element_blank()) +
  scale_x_log10(limits=c(0.002,NA)) + 
  scale_y_log10(limits=c(0.002, NA))


pdf_out("plots/simulated-expression-high-density.pdf", width=1.5, height=1.5)
print(p)
dev.off()
 
