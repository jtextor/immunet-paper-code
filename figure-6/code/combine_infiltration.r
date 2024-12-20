rm(list=ls())

library(scales)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

source("code/config.r")
source("../settings.R")

x.spacing <- 0.70

df_infiltration <- read_csv(paste0(data_folder, "/infiltration.csv"), show_col_types = FALSE)

df_areas <- read_csv(paste0(data_folder, "/regions_per_distance.csv"), show_col_types = FALSE) %>%
    filter(distance==100) %>%
    pivot_longer(c("tumor", "stroma"), names_to="region", values_to="area")


df <- inner_join(df_infiltration, df_areas) %>%
    mutate(density=count/area * 1000 * 1000) %>%
    pivot_wider(id_cols=c("dataset", "slide", "celltype"), names_from="region", values_from="density") %>%
    mutate(dataset = recode(dataset, 
                            "2020-01-31-phenotyping-paper-bladder" = "bladder", 
                            "2020-01-31-phenotyping-paper-melanoma" = "melanoma", 
                            "2020-01-31-phenotyping-paper-prostate" = "prostate", 
                            "2020-02-12-phenotyping-paper-lung-bcell" = "lung"))


df <- df %>% filter(!((celltype == "CD20") & (dataset=="prostate")))

stroma_tumor_counts_p <- ggplot(df, aes(x=stroma, y=tumor, color=celltype)) +
    scale_color_manual(values=c("#E8308A", "#800000")) +
    geom_abline(slope=1) +
    geom_line(aes(group=slide), color="grey") +
    geom_point(size=0.75) +
    xlab(expression(cells/mm^2~"in stroma"))+
    ylab(expression(cells/mm^2~"in tumor"))+
    facet_grid(cols=vars(dataset)) +
    xlim( 1,1e4 ) +
    scale_x_log10(breaks=c(10, 100, 1000, 10000), limits=c(1,1e5), labels = 
                  trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks=c(10, 100, 1000, 10000), limits=c(5,1e4), labels = 
                  trans_format("log10", math_format(10^.x))) +
    theme(legend.position="none", 
          strip.background = element_blank(),
          text=element_text(size=default_pointsize, family = default_font),
          axis.text=element_text(size=default_pointsize, color="black", family=default_font),
          panel.spacing.x = unit(x.spacing, "lines"),
          strip.text.x = element_text(size=default_pointsize, family = default_font), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")
    )

pdf_out(paste0(images_folder, "/stroma_tumor_counts.pdf"), width=7, height=1.6, dpi=72)
print(stroma_tumor_counts_p)
dev.off()
