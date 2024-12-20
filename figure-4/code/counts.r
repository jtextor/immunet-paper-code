rm(list=ls())

library(tidyr)
library(dplyr)
library(readr)
library(irr)
library(ggplot2)
library(cowplot)
library(boot)

source("../settings.R")

set.seed(1009)

df <- read_csv("data/counts.csv")

f <- 15

rmse_stat <- function(d, i) {
  d1 <- d[i,]
  ann <- d1$ann
  pred <- d1$pred

  rmse_inf <- sqrt(mean((ann - pred)**2))
  rmse_inf
}

mean_cl_boot_bca <- function(x, stat_fun, conf.int=0.95, B=2000) {

  b <- boot(x, statistic = stat_fun, R = B)
  bootci <- boot.ci(b, conf = conf.int, "bca")
  m <- b$t0
  low <- bootci$bca[4]
  high <- bootci$bca[5]
  c(m, low, high)
}

metrics_label <- function(rmse, rmse_low, rmse_high, icc, icc_low, icc_high) {
  paste0("rmse=", signif(rmse, 3), " (", signif(rmse_low, 3) ,"-", signif(rmse_high, 3) , ")",
         "\nicc=", signif(icc, 2), " (", signif(icc_low, 2) ,"-", signif(icc_high, 2) , ")")
}

compute_icc <- function(x, y) {
  icc_res <- icc(cbind(x,y), model="twoway", type="agreement", unit="single")
  c(icc_res$value, icc_res$lbound, icc_res$ubound)
}

compute_rmse <- function(x, y) {
  rmse_sign <- mean_cl_boot_bca(tibble(ann=x, pred=y), rmse_stat)
  c(rmse_sign[1], rmse_sign[2], rmse_sign[3])
}

df <- df %>%
    mutate(algorithm=replace(algorithm, algorithm=="inform", "baseline")) %>%
    mutate(algorithm=replace(algorithm, algorithm=="NN", "ImmuNet")) %>%
    mutate(algorithm=factor(algorithm, levels=c("baseline","ImmuNet")))

metrics_df <- df %>%
  group_by(algorithm) %>%
  summarize(icc = list(compute_icc(annotations, detections)),
            rmse = list(compute_rmse(annotations, detections))) %>%
  unnest_wider(c(icc, rmse), names_sep="_") %>%
  rename(icc=icc_1, icc_low=icc_2, icc_high=icc_3,
         rmse=rmse_1, rmse_low=rmse_2, rmse_high=rmse_3)

df_labels <- metrics_df %>% group_by(algorithm) %>%
    summarize(metrics = metrics_label(rmse, rmse_low, rmse_high, icc, icc_low, icc_high))

counts_p <- df %>%
  ggplot(aes(x=annotations, y=detections, fill=algorithm, color=algorithm)) +
  geom_abline(slope=1, linetype="dashed") +
  geom_violin(aes(x=round(annotations/f)*f, group=factor(round(annotations/f)*f   ))) +
  scale_fill_manual(values=c("darkgrey", "black")) +
  scale_color_manual(values=c("darkgrey", "black")) +
  facet_grid(cols=vars(algorithm)) +
  coord_fixed() +
  theme(strip.text.y.left = element_text(angle=0, size=default_pointsize, family=default_font),
        axis.title=element_text(size=default_pointsize, family=default_font),
        axis.text=element_text(size=default_pointsize, color="black", family=default_font),
        axis.line = element_line(color="black"),
        strip.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.spacing.x = unit(0.8, "lines"),
        strip.placement = "outside",
        legend.position="none",
        plot.margin = margin(5.5, 10, 5.5, 5.5, unit="pt"))

pdf_out("plots/counts.pdf", width=3.6, height=2.0, dpi=72)
print(counts_p)
dev.off()

myColor <- c("grey","darkred")
myColor_scale_fill <- scale_fill_gradientn(colours = myColor, trans="sqrt", name="ROIs", c(2, 5, 10))

counts_density_p <- df %>%
  ggplot(aes(x=annotations+1, y=detections+1)) +
  geom_hex(bins=25) +
  geom_abline(slope=1, linetype="dashed") +
  scale_x_log10(breaks=c(1, 2, 11, 101), labels=c("0", "1", "10", "100"), name="annotations") +
  scale_y_log10(breaks=c(1, 2, 11, 101), labels=c("0", "1", "10", "100"), name="detections") +
  geom_text(data=df_labels, aes(label=metrics), x=-0.05, y=2.5, size=default_pointsize*0.36, color="black", lineheight=0.85, hjust=0, family=default_font) +
  geom_text(data=df_labels, aes(label=algorithm), x=-0.05, y=1.9, size=default_pointsize*0.36, color="black", lineheight=0.85, hjust=0, family="Helvetica Neue Medium") +
  myColor_scale_fill +
  facet_wrap(vars(algorithm), strip.position = "top") +
  coord_fixed(clip = "off") +
  guides(fill=guide_colorbar(ticks.colour = NA)) +
  theme(text=element_text(size=default_pointsize, family=default_font),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.key.width=unit(0.20,"cm"),
        legend.key.height=unit(0.5,"cm"),
        legend.margin=margin(0, 3, 0, 0),
        axis.title=element_text(size=default_pointsize, family=default_font),
        axis.text=element_text(size=default_pointsize, color="black", family=default_font),
        axis.line = element_line(color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.spacing.x = unit(1.3, "lines"),
        legend.position="right",
        plot.margin = margin(25, 0, 2, 0, unit="pt"))


pdf_out("plots/counts_density.pdf", width=3.8, height=1.8, dpi=72)
print(counts_density_p)
dev.off()

counts_histogram_p <- df %>% filter(algorithm == "ImmuNet") %>%
  ggplot(aes(x=annotations+1)) +
  geom_histogram(bins=9, color="lightgrey", fill="lightgrey") +
  coord_flip(clip = "off") +
  scale_x_log10(breaks=c(1, 2, 11, 101), labels=c("0", "1", "10", "100"), name="annotations") +
  scale_y_continuous(breaks=c(0, 25, 50), labels = c("0", "", "50"), limit=c(0, 55), name="ROIs") +
  theme(axis.title=element_text(size=default_pointsize, family=default_font),
        axis.text=element_text(size=default_pointsize, color="black", family=default_font),
        axis.line = element_line(color="black"),
        strip.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position="none",
        plot.margin = margin(25, 5, 5.5, 2, unit="pt"))

pdf_out("plots/counts_histogram.pdf", width=2.5, height=2.0, dpi=72)
print(counts_histogram_p)
dev.off()

counts_density_histogram_p <- plot_grid(counts_histogram_p, counts_density_p, ncol=2, rel_widths=c(1.2, 5), align="h", axis="b")

pdf_out("plots/counts_density_histogram.pdf", width=4.3, height=1.8, dpi=72)
print(counts_density_histogram_p)
dev.off()
