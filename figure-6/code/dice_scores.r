rm(list=ls())

library(dplyr)
library(readr)
library(ggplot2)
library(ggbeeswarm)

source("code/config.r")
source("../settings.R")

t <- theme(panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(),
           panel.background = element_blank(), 
           axis.line = element_line(colour = "black"),
           axis.text = element_text(size = 10),
           axis.title = element_text(size = 10),
           strip.background = element_blank(),
           legend.position="none")


load_df <- function(filename, tissue) {

  df <- read_csv(filename)
  
  
  df_thresh <- df %>% group_by(dataset, threshold) %>%
      summarize(dice_score = mean(dice_score)) %>%
      group_by(dataset) %>%
      filter(dice_score == max(dice_score)) %>%
      select(dataset, threshold)
  
  print(tissue)
  print(df_thresh)
  
  df <- inner_join(df, df_thresh) %>%
      mutate(tissue=tissue)
  df
}

tumor_df <- load_df(paste0(data_folder, "/inform_threshold_comp.csv"), "tumor")
stroma_df <- load_df(paste0(data_folder, "/inform_threshold_comp.stroma.csv"), "stroma")

df <- bind_rows(tumor_df, stroma_df) %>%
    mutate(dataset = recode(dataset, "2020-01-31-phenotyping-paper-bladder" = "bladder", "2020-01-31-phenotyping-paper-melanoma" = "melanoma", "2020-01-31-phenotyping-paper-prostate" = "prostate", "2020-02-12-phenotyping-paper-lung-bcell" = "lung"))


print(df)

dice_scores_p <- df %>% ggplot(aes(x=dataset, y=dice_score)) + 
  geom_quasirandom() +
  ylab("dice score") +
  facet_grid(cols=vars(tissue)) + 
  theme(text=element_text(size=11, family = default_font),
        axis.text.x=element_text(size=10)) +
  t

pdf_out(paste0(images_folder, "/dice_scores.pdf"), width=6.5, height=1.5, dpi=72)
print(dice_scores_p)
dev.off()
