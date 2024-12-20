rm(list=ls())

library(tidyr)
library(dplyr)
library(stringr)
library(readr)
library(ggplot2)
library(cowplot)
library(boot)
library(RColorBrewer)
library(polycor)

source("../settings.R")

set.seed(1610)

precision_color <- brewer.pal(9, "Oranges")[8]
recall_color <- brewer.pal(9, "Blues")[8]
fscore_color <- brewer.pal(9, "Greens")[8]

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
  tps <- sum(d1[,1])
  fns <- sum(d1[,2])

  if (tps + fns == 0) {
    return(1)
  } else {
    recall <- tps / (tps + fns)
    return(recall)
  }
}

precision_stat <- function(d, i) {
  d1 <- d[i,]
  tps <- sum(d1[,1])
  fps <- sum(d1[,3])

  if (tps + fps == 0) {
    return(1)
  } else {
    precision <- tps / (tps + fps)
    return(precision)
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

mean_cl_boot_bca <- function(x, stat_fun, conf.int=0.95, B=2000) {

    b <- boot(x, statistic = stat_fun, R = B)
    bootci <- boot.ci(b, conf = conf.int, "bca")
    m <- b$t0
    low <- bootci$bca[4]
    high <- bootci$bca[5]
    c(m, low, high)
}

det_df <- read_csv("data/detection_stat.csv", col_types = cols())
roi_df <- read_csv("data/rois_info.csv", col_types = cols())

det_df <- left_join(det_df, roi_df, by="roi_id")

det_df <- det_df %>% mutate(algorithm=replace(algorithm, algorithm=="inform", "baseline")) %>%
  mutate(algorithm=replace(algorithm, algorithm=="NN", "ImmuNet")) %>%
  mutate(algorithm=factor(algorithm, levels=c("baseline","ImmuNet")))

det_df <- det_df[det_df[,"celltype"] != "T cell",]

det_df <- det_df %>%
  filter(! (celltype=="B cell" & dataset == "prostate")) %>%
  mutate(celltype=replace(celltype, celltype=="CD8", "Cytotoxic T cell")) %>%
  mutate(celltype=replace(celltype, celltype=="CD4", "Helper T cell")) %>%
  mutate(celltype=replace(celltype, celltype=="FOXP3", "Regulatory T cell")) %>%
  mutate(dataset=replace(dataset, dataset=="lung-bcell", "lung"))

df_fscore <- det_df %>%
    group_by(dataset, celltype, algorithm) %>%
    summarize(bs=list(mean_cl_boot_bca(tibble(tp=tp, fn=fn, fp=fp), f_score_stat))) %>%
    unnest_wider(bs, names_sep="_") %>%
    rename(fscore=bs_1, fs_low=bs_2, fs_high=bs_3)

write_csv(df_fscore, "data/fscore.csv")

f_score_p <- df_fscore %>%
  ggplot(aes(x=algorithm, y=fscore, color=celltype, shape=algorithm)) +
  geom_point() +
  geom_errorbar(aes(ymin=fs_low, ymax=fs_high), data=df_fscore, width=0.0) +
  facet_grid(rows=vars(dataset), cols=vars(celltype), switch="y")  +
  scale_y_continuous(breaks = c(0.5, 1.0), labels = scales::percent_format(accuracy = 1)) +
  coord_flip(ylim=c(0.18, 1)) +
  scale_color_manual(values=unname(clrs)[c(2,4,1,5)]) +
  scale_shape_manual(values = c(22,15)) +
  ylab("F-score (95% CI)") +
  theme(text=element_text(size=10, family=default_font),
        strip.text.y.left = element_text(angle=0,
                                         size=default_pointsize, family=default_font),
        strip.text.x = element_text(size=default_pointsize, family=default_font),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=default_pointsize, family=default_font),
        axis.text=element_text(size=default_pointsize, color="black", family=default_font),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.background=element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line(color="black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.spacing.x = unit(1.0, "lines"),
        strip.placement = "outside",
        legend.position="none",
        strip.clip = "off",
        plot.margin = margin(5.5, 10, 5.5, 5.5, unit="pt"))

pdf_out("plots/f_score.pdf", width=6, height=2.5)
print(f_score_p)
dev.off()

cd45ro_df <- read_csv("data/cd45ro_stat.csv",col_types = cols() )
cd45ro_df <- left_join(cd45ro_df, roi_df, by="roi_id")

cd45ro_df <- cd45ro_df %>% mutate(algorithm=replace(algorithm, algorithm=="inform", "baseline")) %>%
  mutate(algorithm=replace(algorithm, algorithm=="NN", "ImmuNet")) %>%
  mutate(algorithm=factor(algorithm, levels=c("baseline","ImmuNet"))) %>%
  group_by(roi_id, algorithm, dataset) %>%
  summarise(ann = list(ann_cd45ro), det = list(det_cd45ro)) %>%
  ungroup()

df_cd45ro_r <- cd45ro_df %>%
  group_by(dataset, algorithm) %>%
  summarize(bs=list(mean_cl_boot_bca(tibble(a = ann, d = det), cor_stat))) %>%
  unnest_wider(bs, names_sep="_") %>%
  rename(r=bs_1, r_low=bs_2, r_high=bs_3)

r_value_p <- df_cd45ro_r %>%
  ggplot(aes(x=algorithm, y=r, shape=algorithm)) +
  geom_point() +
  geom_errorbar(aes(ymin=r_low, ymax=r_high), data=df_cd45ro_r, width=0.0) +
  facet_wrap(vars(dataset), nrow=nrow(df_cd45ro_r) / 2, strip.position="left")  +
  scale_y_continuous( breaks = c(0.0, 1.0), labels = scales::percent_format(accuracy = 1)) +
  coord_flip(ylim=c(0.0, 1)) +
  ggtitle("CD45RO") +
  labs(y = bquote(rho)) +
  scale_shape_manual(values = c(22,15)) +
  theme(text=element_text(size=11, family=default_font),
        plot.title = element_text(hjust = 0.5, size=default_pointsize, family=default_font),
        strip.text.y.left = element_blank(),
        strip.text.x = element_text(size = default_pointsize, family=default_font),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=default_pointsize, family=default_font),
        axis.text=element_text(size=default_pointsize, color="black", family=default_font),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.background=element_blank(),
        axis.line.x = element_line(color="black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.spacing.x = unit(1, "lines"),
        strip.placement = "outside",
        legend.position="none",
        plot.margin = margin(5.5, 10, 5.5, 5.5, unit="pt"))

pdf_out("plots/r_value.pdf", width=6, height=2.5)
print(r_value_p)
dev.off()

f_score_r_value_p <- plot_grid(f_score_p, r_value_p, ncol=2, rel_widths=c(5,0.96), align="h", axis="tb")

pdf_out("plots/f_score_r_value.pdf", width=6, height=2.2)
print(f_score_r_value_p)
dev.off()

# Precision
df_precision <- det_df %>%
  group_by(dataset, celltype, algorithm) %>%
  summarize(bs=list(mean_cl_boot_bca(data.frame(tp, fn, fp), precision_stat))) %>%
  unnest_wider(bs, names_sep="_") %>%
  rename(prec=bs_1, prec_low=bs_2, prec_high=bs_3)

write_csv(df_precision, "data/precision.csv")

precision_p <- df_precision %>%
  ggplot(aes(x=algorithm, y=prec, color=celltype, shape=algorithm)) +
  geom_point() +
  geom_errorbar(aes(ymin=prec_low, ymax=prec_high), data=df_precision, width=0.0) +
  facet_grid(rows=vars(dataset), cols=vars(celltype), switch="y")  +
  scale_y_continuous( breaks = c(0.5, 1.0), labels = scales::percent_format(accuracy = 1)) +
  coord_flip(ylim=c(0.18, 1)) +
  scale_color_manual(values=unname(clrs)[c(2,4,1,5)]) +
  scale_shape_manual(values = c(22,15)) +
  ylab("Precision (95% CI)") +
  theme(text=element_text(size=10, family=default_font),
        strip.text.y.left = element_text(angle=0,
                                         size=default_pointsize, family=default_font),
        strip.text.x = element_text(size=default_pointsize, family=default_font),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=default_pointsize, family=default_font),
        axis.text=element_text(size=default_pointsize, color="black", family=default_font),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.background=element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line(color="black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.spacing.x = unit(1.0, "lines"),
        strip.placement = "outside",
        legend.position="none",
        strip.clip = "off",
        plot.margin = margin(5.5, 10, 5.5, 5.5, unit="pt"))

# Recall
df_recall <- det_df %>%
  group_by(dataset, celltype, algorithm) %>%
  summarize(bs=list(mean_cl_boot_bca(data.frame(tp, fn, fp), recall_stat))) %>%
  unnest_wider(bs, names_sep="_") %>%
  rename(rec=bs_1, rec_low=bs_2, rec_high=bs_3)

write_csv(df_recall, "data/recall.csv")

recall_p <- df_recall %>%
  ggplot(aes(x=algorithm, y=rec, color=celltype, shape=algorithm)) +
  geom_point() +
  geom_errorbar(aes(ymin=rec_low, ymax=rec_high), data=df_recall, width=0.0) +
  facet_grid(rows=vars(dataset), cols=vars(celltype), switch="y")  +
  scale_y_continuous( breaks = c(0.5, 1.0), labels = scales::percent_format(accuracy = 1)) +
  coord_flip(ylim=c(0.18, 1)) +
  scale_color_manual(values=unname(clrs)[c(2,4,1,5)]) +
  scale_shape_manual(values = c(22,15)) +
  ylab("Recall (95% CI)") +
  theme(text=element_text(size=10, family=default_font),
        strip.text.y.left = element_text(angle=0,
                                         size=default_pointsize, family=default_font),
        strip.text.x = element_text(size=default_pointsize, family=default_font),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=default_pointsize, family=default_font),
        axis.text=element_text(size=default_pointsize, color="black", family=default_font),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.background=element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line(color="black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.spacing.x = unit(1.0, "lines"),
        strip.placement = "outside",
        legend.position="none",
        strip.clip = "off",
        plot.margin = margin(5.5, 10, 5.5, 5.5, unit="pt"))

df_fscore <- df_fscore %>%
  rename(mid=fscore, low=fs_low, high=fs_high) %>%
  mutate(type="fscore")

df_precision <- df_precision %>%
  rename(mid=prec, low=prec_low, high=prec_high) %>%
  mutate(type="precision")

df_recall <- df_recall %>%
  rename(mid=rec, low=rec_low, high=rec_high) %>%
  mutate(type="recall")


df_all <- rbind(df_fscore, df_precision, df_recall)
df_all <- df_all %>%
  mutate(shape=str_c(algorithm, type)) %>%
  mutate(line=str_c(dataset, type)) %>%
  mutate(shape=factor(shape, levels=c("ImmuNetrecall", "baselinerecall", "ImmuNetprecision", "baselineprecision",
                                      "ImmuNetfscore", "baselinefscore")))

row_lbl <- function(str) {
  should_show <- grepl("precision", str, fixed = TRUE)
  if_else(should_show, str_remove(str, "precision"), "")
}

all_p <- df_all %>%
  ggplot(aes(x=algorithm, y=mid, color=type, shape=algorithm)) +
  geom_point() +
  geom_errorbar(aes(ymin=low, ymax=high), data=df_all, width=0.0) +
  facet_grid(line ~ celltype, switch="y", labeller = labeller(line=row_lbl))  +
  scale_y_continuous( breaks = c(0.5, 1.0), labels = scales::percent_format(accuracy = 1)) +
  coord_flip(ylim=c(0.18, 1)) +
  scale_color_manual(values=c(fscore_color, precision_color, recall_color)) +
  scale_shape_manual(values = c(22, 15)) +
  ylab("Scores (95% CI)") +
  theme(text=element_text(size=10, family=default_font),
        strip.text.y.left = element_text(angle=0,
                                         size=default_pointsize, family=default_font),
        strip.text.x = element_text(size=default_pointsize, family=default_font),
        axis.title.y=element_blank(),
        axis.title.x=element_text(size=default_pointsize, family=default_font),
        axis.text=element_text(size=default_pointsize, color="black", family=default_font),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.background=element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line(color="black"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.spacing.x = unit(1.0, "lines"),
        legend.position="rigth",
        strip.clip = "off",
        plot.margin = margin(0, 0, 0, 0, unit="pt")
  )

all_r_value_p <- plot_grid(all_p, r_value_p, ncol=2, rel_widths=c(5, 0.96), align="h", axis="tb")

pdf_out("plots/roi_evaluation.pdf", width=6, height=3.5)
print(all_r_value_p)
dev.off()
