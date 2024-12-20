rm(list=ls())
source("code/config.r")

library(stringr)
library(stringi)
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(forcats)

slides_folders <- list.dirs(path = slides_data_folder, full.names = TRUE, recursive = FALSE)

slides_csv_globs <- slides_folders %>%
  paste0(., "/*.csv")

files <- slides_csv_globs %>% 
  Sys.glob(.) %>% 
  str_subset(pattern="neighbors", negate=F)

df <- files %>%
  set_names() %>% 
  map_dfr(function(filename) read_csv(filename, col_types = cols()), .id="source") %>%
  rowwise() %>%
  mutate(dataset=stri_split_fixed(source, "/")[[1]][3], slide=stri_split_fixed(basename(source), ".")[[1]][1]) %>%
  ungroup() %>%
  select(-source) 

df <- df %>%
    mutate(neigh = CD8_neigh + CD4_neigh + FOXP3_neigh + CD20_neigh,
           CD8_neigh = CD8_neigh + 0.1,
           CD4_neigh = CD4_neigh + 0.1,
           FOXP3_neigh = FOXP3_neigh + 0.1,
           CD20_neigh = CD20_neigh + 0.1,
           CD45RO_neigh = CD45RO_neigh + 0.01,
           ) %>%
    filter(neigh>0, distance > -100, distance<100) 

avg_neigh_df <- df %>%
    group_by(dataset, slide) %>%
    summarize(
              avg_CD4_neigh = mean(CD4_neigh),
              avg_CD8_neigh = mean(CD8_neigh),
              avg_FOXP3_neigh = mean(FOXP3_neigh),
              avg_CD20_neigh = mean(CD20_neigh),
              avg_CD45RO_neigh = mean(CD45RO_neigh),
    ) 

df_cd45ro <- df %>% 
  filter(celltype %in% c("CD3", "CD4", "CD8", "FOXP3")) %>%
  select(-celltype) %>%
  rename(celltype = cd45ro_expression)

df <- bind_rows(df, df_cd45ro)

df <- inner_join(df, avg_neigh_df) %>%
    mutate(
           CD20_neigh = CD20_neigh / avg_CD20_neigh,
           CD4_neigh = CD4_neigh / avg_CD4_neigh,
           CD8_neigh = CD8_neigh / avg_CD8_neigh,
           FOXP3_neigh = FOXP3_neigh / avg_FOXP3_neigh,
           CD45RO_neigh = CD45RO_neigh / avg_CD45RO_neigh,
    ) %>% 
    select(-starts_with("avg")) %>%
    pivot_longer(cols=ends_with("_neigh"), names_to="neighbor_type", values_to="neighbors") %>%
    mutate(celltype = fct_relevel(celltype, 
                              c("hi", "mid", "low", "FOXP3", "CD8", 
                              "CD4", "CD20"))) %>%
    mutate(neighbor_type = recode(neighbor_type, "CD20_neigh" = "CD20", 
                                  "CD4_neigh" = "CD4", "CD8_neigh" = "CD8", 
                                  "FOXP3_neigh" = "FOXP3", 
                                  "CD45RO_neigh" = "CD45RO")) %>%
    mutate(neighbor_type = fct_relevel(neighbor_type, 
                              "CD20", "CD4", "CD8", 
                              "FOXP3", "CD45RO")) %>%
    select(dataset, slide, celltype, neighbor_type, neighbors, distance) 


df_infiltration <- read_csv(paste0(data_folder, "/infiltration.csv"), show_col_types = FALSE)

df_areas <- read_csv(paste0(data_folder, "/regions_per_distance.csv"), show_col_types = FALSE) %>%
    filter(distance==100) %>%
    pivot_longer(c("tumor", "stroma"), names_to="region", values_to="area")

df_inf <- inner_join(df_infiltration, df_areas) %>%
    mutate(density=count/area * 1000 * 1000) %>%
    pivot_wider(id_cols=c("dataset", "slide", "celltype"), names_from="region", values_from="density") %>%
    mutate(infiltration = tumor/stroma) %>% 
    filter(celltype=="CD3") %>%
    select(dataset, slide, infiltration, tumor, stroma)

df <- inner_join(df_inf, df) %>%
    mutate(dataset = recode(dataset, "2020-01-31-phenotyping-paper-bladder" = "bladder", "2020-01-31-phenotyping-paper-melanoma" = "melanoma", "2020-01-31-phenotyping-paper-prostate" = "prostate", "2020-02-12-phenotyping-paper-lung-bcell" = "lung"))

save.image(file = paste0(data_folder, "/neighbor_per_cell_workspace.RData"))
