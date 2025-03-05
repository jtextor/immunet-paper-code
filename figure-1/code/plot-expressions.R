rm(list=ls())

library(dplyr)
library(readr)
library(ggplot2)
library(RColorBrewer)

source("../settings.R")

df <- read_csv("data/expression_data.csv.gz")

myColor <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
myColor <- c("#D3D3D3", myColor)

myColor <- c("grey","darkblue","lightblue", "green", "yellow","orange", "red")
myColor_scale_fill <- scale_fill_gradientn(colours = myColor, trans="sqrt")

nrows <- nrow(df)
rands <- 1 + 0.3 * rnorm(nrows)
df$cd3 <- (df$cd3+0.01) * rands

rands <- 1 + 0.3 * rnorm(nrows)
df$cd20 <- (df$cd20+0.01) * rands

p <- df %>%
ggplot(aes(cd3, cd20)) +
  stat_binhex(bins=55) +
  xlab("CD3") +
  ylab("CD20") +
  facet_wrap(~simulation, nrow=1,labeller = labeller(simulation = 
    c(
      "0" = "10%",
      "1" = "20%",
      "2" = "30%",
      "3" = "40%",
      "4" = "50%",
      "5" = "60%",
      "6" = "70%",
      "7" = "80%",
      "8" = "90%",
      "9" = "100%"
      )
  )) +
  myColor_scale_fill +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),  
        strip.text.x = element_text(size=8, family="sans"), 
        strip.background = element_blank(),
        text = element_text(size=8, family="sans"),
        legend.position = "none",
        axis.text = element_blank(), 
        axis.ticks = element_blank()) +
  scale_x_log10(limits=c(0.002,NA)) + 
  scale_y_log10(limits=c(0.002, NA))

dir.create("plots", showWarnings = FALSE)

pdf_out("plots/expressions.pdf", width=6.5, height=1.2)
print(p)
dev.off()

