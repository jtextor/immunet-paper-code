rm(list=ls())
source("code/config.r")

library(tidyr)
library(stringr)
library(stringi)
library(readr)
library(dplyr)
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

# data is already transformed into um, i.e. distance and x,y are in correct units
df <- df %>%
    mutate(neigh = CD8_neigh + CD4_neigh + FOXP3_neigh + CD20_neigh,
           CD8_neigh = CD8_neigh + 0.1,
           CD4_neigh = CD4_neigh + 0.1,
           FOXP3_neigh = FOXP3_neigh + 0.1,
           CD20_neigh = CD20_neigh + 0.1,
           CD45RO_neigh = CD45RO_neigh + 0.01,
           ) %>%
    filter(neigh > 0, distance > -100, distance < 100) 

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
    mutate(neighbor_type = recode(neighbor_type, "CD20_neigh" = "CD20", "CD4_neigh" = "CD4", "CD8_neigh" = "CD8", "FOXP3_neigh" = "FOXP3", "CD45RO_neigh" = "CD45RO")) %>%
    mutate(neighbor_type = fct_relevel(neighbor_type, 
                              "CD20", "CD4", "CD8", 
                              "FOXP3", "CD45RO")) %>%
    group_by(dataset, slide, celltype, neighbor_type) %>%
    summarize(neighbors = mean(neighbors)) %>%
    group_by(dataset, celltype, neighbor_type) %>%
    filter(neighbors > 0) %>%
    summarize(neighbor_spread = sd(log(neighbors)), neighbors = mean(neighbors))  %>%
    mutate(dataset = recode(dataset, "2020-01-31-phenotyping-paper-bladder" = "bladder", "2020-01-31-phenotyping-paper-melanoma" = "melanoma", "2020-01-31-phenotyping-paper-prostate" = "prostate", "2020-02-12-phenotyping-paper-lung-bcell" = "lung")) %>%
    mutate(celltype = recode(celltype, 
                             "CD20"="CD20^'+'",
                             "CD4"="CD3^'+'~CD8^'-'",
                             "CD8"="CD8^'+'",
                             "FOXP3"="FOXP3^'+'",
                             "low"="CD3^'+'~CD45RO^'low'",
                             "mid"="CD3^'+'~CD45RO^'mid'",
                             "hi"="CD3^'+'~CD45RO^'hi'",
                             )) %>%
    mutate(neighbor_type = recode(neighbor_type, 
                             "CD20"="CD20^'+'",
                             "CD4"="CD3^'+'~CD8^'-'",
                             "CD8"="CD8^'+'",
                             "FOXP3"="FOXP3^'+'",
                             ))

    
save.image(file = paste0(data_folder, "/neighbor_workspace.RData"))
