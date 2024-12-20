rm(list=ls())

library(tidyr)
library(dplyr)
library(cowplot)
library(readr)
library(grid)
library(gridExtra)
library(scales)
library(boot)
library(RColorBrewer)
library(egg)
library(polycor)

source("../settings.R")

set.seed(1610)

fscore_color <- brewer.pal(9, "Greens")[8]
precision_color <- brewer.pal(9, "Oranges")[8]
recall_color <- brewer.pal(9, "Blues")[8]

f_score_stat <- function(d, i) {
  d1 <- d[i,]
  tps <- sum(d1$tp)
  fns <- sum(d1$fn)
  fps <- sum(d1$fp)

  if (tps + fps + fns == 0) {
    return(1)
  } else {
    f_score <- 2*tps/(2*tps + fps + fns)
    return(f_score)
  }
}

recall_stat <- function(d, i) {
  d1 <- d[i,]
  tps <- sum(d1$tp)
  fns <- sum(d1$fn)

  if (tps + fns == 0) {
    return(1)
  } else {
    sensitivity <- tps / (tps + fns)
    return(sensitivity)
  }
}

precision_stat <- function(d, i) {
  d1 <- d[i,]
  tps <- sum(d1$tp)
  fps <- sum(d1$fp)

  if (tps + fps == 0) {
    return(1)
  } else {
    sensitivity <- tps / (tps + fps)
    return(sensitivity)
  }
}

cor_stat <- function(d, i) {
  d1 <- d[i,]

  annotated <- unlist(d1$a)
  detected <- unlist(d1$d)
  detected_continuous <- length(unique(detected)) > 2

  cor_value <- if (detected_continuous) polyserial(detected, annotated) else polychor (detected, annotated)
  cor_value
}

mean_cl_boot_bca <- function(x, stat_fun, conf.int=0.95, B=2000 ) {
    if (nrow(x) < 6) {
        m <- stat_fun(x, 1:nrow(x))
        low <- 0
        high <- 1

    }else{

    b <- boot(x, statistic = stat_fun, R = B)
    bootci <- boot.ci(b, conf = conf.int, "bca")
    m <- b$t0
    low <- bootci$bca[4]
    high <- bootci$bca[5]
    }

  c(m, low, high)
}

det_df <- read_csv("data/detection_stat.csv")
roi_df <- read_csv("data/rois_info.csv")

det_df <- left_join(det_df, roi_df, by="roi_id")

det_df <- det_df %>% mutate(algorithm=replace(algorithm, algorithm=="inform", "baseline")) %>%
  mutate(algorithm=replace(algorithm, algorithm=="NN", "ImmuNet")) %>%
  mutate(algorithm=factor(algorithm, levels=c("baseline","ImmuNet"))) %>%
  mutate(density = annotations / area)

det_df <- det_df[det_df[,"celltype"] != "T cell",]

det_df <- det_df %>% filter(! (celltype=="B cell" & dataset == "prostate")) %>%
  mutate(celltype=replace(celltype, celltype=="CD8", "Cytotoxic T cell")) %>%
  mutate(celltype=replace(celltype, celltype=="CD4", "Helper T cell")) %>%
  mutate(celltype=replace(celltype, celltype=="FOXP3", "Regulatory T cell"))

density_hist_p <- det_df %>%
  filter(algorithm == "ImmuNet") %>%
  ggplot(aes(x=density, y=tp)) +
  scale_colour_manual(values=c("#666666", "black")) +
  facet_grid(cols=vars(celltype), switch="y")  +
  geom_hline(yintercept = c(1), color="grey", size=.25) +
  stat_summary_bin(fun="length", geom="line", bins=10, color="red") +
  scale_x_log10(breaks=c(10^2, 10^3, 10^4, 10^5), labels = trans_format("log10", math_format(10^.x))) +
  ylab("count") +
  xlab("F-score (95% CI)") +
  theme(text=element_text(size=11, family=default_font),
        strip.text.y.left = element_text(angle=0, size=10, family=default_font),
        axis.title.y=element_text(size=10, family=default_font),
        axis.title.x=element_blank(),
        axis.text=element_text(size=10, color="black", family=default_font),
        axis.line = element_line(color="black"),
        strip.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.spacing.x = unit(0.8, "lines"),
        strip.placement = "outside",
        legend.position="none",
        plot.margin = margin(5.5, 10, 5.5, 5.5, unit="pt")
  )

pdf_out("plots/density_histogram.pdf", width=6, height=2.0)
print(density_hist_p)
dev.off()

# Find breaks for density binning
q <- 6

min.density <- min(det_df$density[det_df$density > 0])
max.density <- max(det_df$density)

print(min.density)
print(max.density)

breaks <- exp(seq(log(min.density - 0.1), log(max.density), length.out=q+1))
breaks <- breaks[c(T,F,F,T,T,T)]
print(breaks)

det_df <- det_df %>%
  mutate(new_bin = cut(density, breaks=breaks))

fscore_density_df <- det_df %>%
  group_by(celltype, algorithm, new_bin) %>%
  summarize(bs=list(mean_cl_boot_bca(tibble(tp=tp, fn=fn, fp=fp), f_score_stat)), num=length(tp), density=mean(density)) %>%
  unnest_wider(bs, names_sep="_") %>%
  rename(fs=bs_1, fs_low=bs_2, fs_high=bs_3) %>%
  filter(num > 6)

# Shift baseline values slightly to avoid overlap of error bars on plots
fscore_density_df <- fscore_density_df %>%
  mutate(density = if_else(algorithm=="baseline", density * 1.1, density))

density_fscore_p <- fscore_density_df %>%
  ggplot(aes(x=density, y=fs, shape=algorithm, group=algorithm)) +
  geom_line(color=fscore_color) +
  geom_point(color=fscore_color) +
  scale_shape_manual(values = c(22,15)) +
  geom_errorbar(aes(ymin=fs_low, ymax=fs_high), width=0.0, color=fscore_color) +
  facet_grid(cols=vars(celltype), switch="y")  +
  xlab(expression(cells/mm^2)) +
  ylab("F-score (95% CI) \n ") +
  scale_y_continuous( breaks = c(0.0, 0.5, 1.0), limits=c(0,1),
	labels = scales::percent_format(accuracy = 1)) +
  scale_x_log10(breaks=c(10^2, 10^3, 10^4, 10^5), limits=c(800, NA),  labels = trans_format("log10", math_format(10^.x))) +
  theme(text=element_text(size=default_pointsize, family=default_font),
        strip.text.y.left = element_text(angle=0, size= default_pointsize,
                                         family=default_font),
        strip.text.x = element_text(size=default_pointsize, family=default_font),
        axis.title.y=element_text(size=default_pointsize, family=default_font),
        axis.title.x=element_blank(),
        axis.text=element_text(size=default_pointsize, color="black", family=default_font),
        axis.line = element_line(color="black"),
        panel.background = element_blank(),
        strip.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0.8, "lines"),
        strip.placement = "inside",
        strip.clip = "off",
        legend.position="none",
        plot.margin = margin(0, 10, 0, 5.5, unit="pt"))

pdf_out("plots/density_vs_accuracy.pdf", width=6, height=2.0)
print(density_fscore_p)
dev.off()

cd45ro_df <- read_csv("data/cd45ro_stat.csv")

cd45ro_df <- cd45ro_df %>% mutate(algorithm=replace(algorithm, algorithm=="inform", "baseline")) %>%
  mutate(algorithm=replace(algorithm, algorithm=="NN", "ImmuNet")) %>%
  mutate(algorithm=factor(algorithm, levels=c("baseline","ImmuNet"))) %>%
  group_by(roi_id, algorithm) %>%
  summarise(ann = list(ann_cd45ro), det = list(det_cd45ro)) %>%
  ungroup()

cd45ro_df <- left_join(cd45ro_df, roi_df, by="roi_id") %>%
  mutate(density = annotations / area)

cd45ro_df <- cd45ro_df %>%
  mutate(new_bin = cut(density, breaks=breaks))

rvalue_density_df <- cd45ro_df %>%
  group_by(algorithm, new_bin) %>%
  summarize(bs=list(mean_cl_boot_bca(tibble(a = ann, d = det), cor_stat)), num=length(density), density=mean(density)) %>%
  unnest_wider(bs, names_sep="_") %>%
  rename(r=bs_1, r_low=bs_2, r_high=bs_3) %>%
  filter(num > 6)

# Shift baseline values slightly to avoid overlap of error bars on plots
rvalue_density_df <- rvalue_density_df %>%
  mutate(density = if_else(algorithm=="baseline", density*1.1, density))

density_rscore_p <- rvalue_density_df %>%
  ggplot(aes(x=density, y=r, shape=algorithm)) +
  geom_line() +
  geom_point() +
  scale_shape_manual(values = c(22,15)) +
  labs(y = bquote(rho)) +
  geom_errorbar(aes(ymin=r_low, ymax=r_high), width=0.0) +
  scale_x_log10(breaks=c(10^2, 10^3, 10^4, 10^5),  limits=c(1000, NA), labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous( breaks = c(0.0, 0.5, 1.0), limits=c(0,1),
	labels=scales::percent_format(accuracy = 1)) +
  ggtitle("CD45RO") +
  theme(plot.title = element_text(hjust = 0.5, size=default_pointsize)) +
  theme(text=element_text(size=11, family=default_font),
        strip.text.y.left = element_blank(),
        strip.text.x = element_text(size = default_pointsize, family=default_font),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=default_pointsize, color="black", family=default_font),
        axis.text=element_text(size=default_pointsize, color="black", family=default_font),
        axis.text.y=element_blank(),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        strip.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.spacing.x = unit(0.5, "lines"),
        strip.placement = "inside",
        legend.position="none",
        plot.margin = margin(0, 10, 0, 5.5, unit="pt")
    )

pdf_out("plots/density_vs_cd45ro_r.pdf", width=6, height=2.5)
print(density_rscore_p)
dev.off()

x.grob <- textGrob(expression(lymphocytes/mm^2),
                   gp=gpar(fontface="plain", col="black", fontsize=default_pointsize,
                           fontfamily=default_font))

#add to plot
density_frscore_p <- plot_grid(density_fscore_p, density_rscore_p, ncol=2, rel_widths=c(5,1.1), align="h", axis="tb")

g <- arrangeGrob(density_frscore_p, bottom = x.grob, padding=unit(0.1, "line"))

pdf_out("plots/density_f_score_r_value.pdf", width=6, height=1.3)
grid.draw(g)
dev.off()

# Precision
precision_density_df <- det_df %>%
  group_by(celltype, algorithm, new_bin) %>%
  summarize(bs=list(mean_cl_boot_bca(tibble(tp=tp, fp=fp), precision_stat)), num=length(tp), density=mean(density)) %>%
  unnest_wider(bs, names_sep="_") %>%
  rename(prec=bs_1, prec_low=bs_2, prec_high=bs_3) %>%
  filter(num > 6)

# Shift baseline values slightly to avoid overlap of error bars on plots
precision_density_df <- precision_density_df %>%
  mutate(density = if_else(algorithm=="baseline", density*1.1, density))

density_precision_p <- precision_density_df %>%
  ggplot(aes(x=density, y=prec, shape=algorithm, group=algorithm)) +
  geom_line(color=precision_color) +
  geom_point(color=precision_color) +
  scale_shape_manual(values = c(22,15)) +
  geom_errorbar(aes(ymin=prec_low, ymax=prec_high), width=0.0, color=precision_color) +
  facet_grid(cols=vars(celltype), switch="y")  +
  xlab(expression(cells/mm^2)) +
  ylab("Precision") +
  scale_y_continuous(breaks = c(0.0, 0.5, 1.0), limits=c(0,1),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_x_log10(breaks=c(10^2, 10^3, 10^4, 10^5), limits=c(800, NA),  labels = trans_format("log10", math_format(10^.x))) +
  theme(text=element_text(size=default_pointsize, family=default_font),
        strip.text.y.left = element_text(angle=0, size= default_pointsize,
                                         family=default_font),
        strip.text.x = element_text(size=default_pointsize, family=default_font),
        axis.title.y=element_text(size=default_pointsize, family=default_font),
        axis.title.x=element_blank(),
        axis.text=element_text(size=default_pointsize, color="black", family=default_font),
        axis.line = element_line(color="black"),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        strip.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0.8, "lines"),
        strip.placement = "outside",
        strip.clip = "off",
        legend.position="none",
        plot.margin = margin(5, 10, 5, 5.5, unit="pt")
  )

# Recall
recall_density_df <- det_df %>%
  group_by(celltype, algorithm, new_bin) %>%
  summarize(bs=list(mean_cl_boot_bca(tibble(tp=tp, fn=fn), recall_stat)), num=length(tp), density=mean(density)) %>%
  unnest_wider(bs, names_sep="_") %>%
  rename(rec=bs_1, rec_low=bs_2, rec_high=bs_3) %>%
  filter(num > 6)

# Shift baseline values slightly to avoid overlap of error bars on plots
recall_density_df <- recall_density_df %>%
  mutate(density = if_else(algorithm=="baseline", density*1.1, density))

density_recall_p <- recall_density_df %>%
  ggplot(aes(x=density, y=rec, shape=algorithm, group=algorithm)) +
  geom_line(color=recall_color) +
  geom_point(color=recall_color) +
  scale_shape_manual(values = c(22,15)) +
  geom_errorbar(aes(ymin=rec_low, ymax=rec_high), width=0.0, color=recall_color) +
  facet_grid(cols=vars(celltype), switch="y")  +
  xlab(expression(lymphocytes/mm^2)) +
  ylab("Recall") +
  scale_y_continuous(breaks = c(0.0, 0.5, 1.0), limits=c(0,1),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_x_log10(breaks=c(10^2, 10^3, 10^4, 10^5), limits=c(800, NA),  labels = trans_format("log10", math_format(10^.x))) +
  theme(text=element_text(size=default_pointsize, family=default_font),
        strip.text = element_blank(),
        axis.title=element_text(size=default_pointsize, family=default_font),
        axis.text=element_text(size=default_pointsize, color="black", family=default_font),
        axis.line = element_line(color="black"),
        panel.background = element_blank(),
        strip.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0.8, "lines"),
        strip.clip = "off",
        legend.position="none",
        plot.margin = margin(5, 10, 5, 5.5, unit="pt")
  )

pdf_out("supporting-figure.pdf", width=6, height=2.5)
ggarrange(density_precision_p, density_recall_p, nrow=2, newpage=FALSE)
dev.off()

