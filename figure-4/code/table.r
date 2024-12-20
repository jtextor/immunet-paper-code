rm(list=ls())

library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(boot)
library(readr)

source("../settings.R")

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

pearson_cor_stat <- function(d, i) {
  d1 <- d[i,]

  var1 <- unlist(d1[,1])
  var2 <- unlist(d1[,2])

  r_value <- cor(var1, var2, method="pearson")
  return(r_value)
}

mean_cl_boot_bca <- function(x, stat_fun, conf.int=0.95, B=5000) {

  b <- boot(x, statistic = stat_fun, R = B)
  bootci <- boot.ci(b, conf = conf.int, "bca")
  m <- b$t0
  low <- bootci$bca[4]
  high <- bootci$bca[5]
  c(m, low, high)
}

det_df <- read_csv("data/detection_stat.csv")
roi_df <- read_csv("data/rois_info.csv")

det_df <- left_join(det_df, roi_df, by="roi_id")

det_df <- det_df %>% mutate(algorithm=replace(algorithm, algorithm=="inform", "baseline")) %>%
  mutate(algorithm=replace(algorithm, algorithm=="NN", "ImmuNet")) %>%
  mutate(algorithm=factor(algorithm, levels=c("baseline","ImmuNet")))


det_df <- det_df %>%
  filter(! (celltype=="B cell" & dataset == "prostate")) %>%
  mutate(celltype=replace(celltype, celltype=="T cell", "CD3+")) %>%
  mutate(celltype=replace(celltype, celltype=="B cell", "CD20+")) %>%
  mutate(celltype=replace(celltype, celltype=="CD8", "CD8+")) %>%
  mutate(celltype=replace(celltype, celltype=="CD4", "CD3+CD8-FOXP3-")) %>%
  mutate(celltype=replace(celltype, celltype=="FOXP3", "FOXP3+")) %>%
  mutate(dataset=replace(dataset, dataset=="lung-bcell", "lung"))

df_means <- det_df %>%
  group_by(dataset, celltype, algorithm) %>%
  summarize(bs=list(mean_cl_boot_bca(tibble(tp=tp, fn=fn, fp=fp), f_score_stat))) %>%
  unnest_wider(bs, names_sep="_") %>%
  rename(fscore=bs_1, fscore.min=bs_2, fscore.max=bs_3) %>%
  mutate(fscore.label = paste0(
    format(round(fscore, digits=2), nsmall=1), " (",
    format(round(fscore.min, digits=2), nsmall=1), "-",
    format(round(fscore.max, digits=2), nsmall=1), ")"
  ))

cd45ro_df <- read_csv("data/cd45ro_stat.csv")
cd45ro_df <- left_join(cd45ro_df, roi_df, by="roi_id")

cd45ro_df <- cd45ro_df %>% mutate(algorithm=replace(algorithm, algorithm=="inform", "baseline")) %>%
  mutate(algorithm=replace(algorithm, algorithm=="NN", "ImmuNet")) %>%
  mutate(algorithm=factor(algorithm, levels=c("baseline","ImmuNet"))) %>%
  group_by(roi_id, algorithm, dataset) %>%
  summarise(ann = list(ann_cd45ro), det = list(det_cd45ro)) %>%
  ungroup()

df_cd45ro_r <- cd45ro_df %>%
  group_by(dataset, algorithm) %>%
  summarize(bs=list(mean_cl_boot_bca(tibble(a = ann, b = det), pearson_cor_stat))) %>%
  unnest_wider(bs, names_sep="_") %>%
  rename(fscore=bs_1, fscore.min=bs_2, fscore.max=bs_3)

df_cd45ro_r <- df_cd45ro_r %>%
  mutate(celltype="CD45RO") %>%
  mutate(fscore.label = paste0(
    format(round(fscore, digits=2), nsmall=1), " (",
    format(round(fscore.min, digits=2), nsmall=1), "-",
    format(round(fscore.max, digits=2), nsmall=1), ")"
  ))

df <- bind_rows(df_means, df_cd45ro_r)

df$celltype = factor(df$celltype, levels=c("CD20+", "CD3+", "CD3+CD8-FOXP3-", "CD8+", "FOXP3+", "CD45RO"))

write_csv(df, "data/f_scores.csv")

f_score_p <- df %>%
    ggplot(aes(y=algorithm, x=0, label=fscore.label, color=algorithm)) +
    geom_text(size=8*0.36) +
    facet_grid(rows=vars(dataset), cols=vars(celltype), switch="y") +
    scale_colour_manual(values=c("darkgrey", "black")) +
    theme(text=element_text(size=11, family=default_font),
          strip.text.y.left = element_text(angle=0, size=10, family=default_font),
          axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          strip.background=element_blank(),
          axis.line = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.spacing.x = unit(1.0, "lines"),
          strip.placement = "outside",
          legend.position="none",
          plot.margin = margin(5.5, 10, 5.5, 5.5, unit="pt")
    )


pdf_out("plots/f_score.table.pdf", width=8, height=2.7)
print(f_score_p)
dev.off()

# Precision
df_means <- det_df %>%
  group_by(dataset, celltype, algorithm) %>%
  summarize(bs=list(mean_cl_boot_bca(tibble(tp=tp, fp=fp), precision_stat))) %>%
  unnest_wider(bs, names_sep="_") %>%
  rename(prec=bs_1, prec.min=bs_2, prec.max=bs_3) %>%
  mutate(prec.label = paste0(
    format(round(prec, digits=2), nsmall=1), " (",
    format(round(prec.min, digits=2), nsmall=1), "-",
    format(round(prec.max, digits=2), nsmall=1), ")"
  ))

df_means$celltype = factor(df_means$celltype, levels=c("CD20+", "CD3+", "CD3+CD8-FOXP3-", "CD8+", "FOXP3+"))

precision_p <- df_means %>%
  ggplot(aes(y=algorithm, x=0, label=prec.label, color=algorithm)) +
  geom_text(size=8*0.36) +
  facet_grid(rows=vars(dataset), cols=vars(celltype), switch="y") +
  scale_colour_manual(values=c("darkgrey", "black")) +
  theme(text=element_text(size=11, family=default_font),
        strip.text.y.left = element_text(angle=0, size=10, family=default_font),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        strip.background=element_blank(),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.spacing.x = unit(1.0, "lines"),
        strip.placement = "outside",
        legend.position="none",
        plot.margin = margin(5.5, 10, 5.5, 5.5, unit="pt")
  )


pdf_out("plots/prec.table.pdf", width=8, height=2.7)
print(precision_p)
dev.off()


# Recall
df_means <- det_df %>%
  group_by(dataset, celltype, algorithm) %>%
  summarize(bs=list(mean_cl_boot_bca(tibble(tp=tp, fn=fn), recall_stat))) %>%
  unnest_wider(bs, names_sep="_") %>%
  rename(rec=bs_1, rec.min=bs_2, rec.max=bs_3) %>%
  mutate(rec.label = paste0(
    format(round(rec, digits=2), nsmall=1), " (",
    format(round(rec.min, digits=2), nsmall=1), "-",
    format(round(rec.max, digits=2), nsmall=1), ")"
  ))

df_means$celltype = factor(df_means$celltype, levels=c("CD20+", "CD3+", "CD3+CD8-FOXP3-", "CD8+", "FOXP3+"))

rec_p <- df_means %>%
  ggplot(aes(y=algorithm, x=0, label=rec.label, color=algorithm)) +
  geom_text(size=8*0.36) +
  facet_grid(rows=vars(dataset), cols=vars(celltype), switch="y") +
  scale_colour_manual(values=c("darkgrey", "black")) +
  theme(text=element_text(size=11, family=default_font),
        strip.text.y.left = element_text(angle=0, size=10, family=default_font),
        axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        strip.background=element_blank(),
        axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.spacing.x = unit(1.0, "lines"),
        strip.placement = "outside",
        legend.position="none",
        plot.margin = margin(5.5, 10, 5.5, 5.5, unit="pt")
  )


pdf_out("plots/rec.table.pdf", width=8, height=2.7)
print(rec_p)
dev.off()
