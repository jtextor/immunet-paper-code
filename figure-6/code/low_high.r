rm(list=ls())
source("code/config.r")

library(readr)
library(dplyr)
library(tidyr)

slides_low <- c("bladder10", 
                "melanoma09", 
                "prostate08", 
                "lung07")

slides_high <- c("bladder01",
                  "melanoma22",
                  "prostate05",
                  "lung01")

slides_to_show <- c(slides_low, slides_high)

df_infiltration <- read_csv(paste0(data_folder, "/infiltration.csv"), show_col_types = FALSE)

df_areas <- read_csv(paste0(data_folder, "/regions_per_distance.csv"), show_col_types = FALSE) %>%
    filter(distance==500) %>%
    pivot_longer(c("tumor", "stroma"), names_to="region", values_to="area")

print(df_infiltration)
print(df_areas)

df <- inner_join(df_infiltration, df_areas) %>%
    mutate(density=count/area * 1000 * 1000) %>%
    pivot_wider(id_cols=c("dataset", "slide", "celltype"), names_from="region", values_from="density") %>%
    mutate(tissue = recode(dataset, "2020-01-31-phenotyping-paper-bladder" = "bladder", "2020-01-31-phenotyping-paper-melanoma" = "melanoma", "2020-01-31-phenotyping-paper-prostate" = "prostate", "2020-02-12-phenotyping-paper-lung-bcell" = "lung"))
    

df_min_max <- df %>% group_by(dataset, celltype) %>% 
    summarize(min_til = min(tumor, na.rm=T), max_til=max(tumor, na.rm=T))
print(df_min_max)

df <- inner_join(df, df_min_max) %>%
    filter(celltype == "CD3", slide %in% slides_to_show) %>%
    mutate(low=slide %in% slides_low)

print(df)

write_csv(df, paste0(data_folder, "/low_high.csv"))

filenames = paste0(df$dataset, ".", df$slide, ".csv")
print(filenames)
