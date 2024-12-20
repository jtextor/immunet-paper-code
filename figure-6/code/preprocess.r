rm(list=ls())
source("code/config.r")

library(stringr)
library(stringi)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)

slides_folders <- list.dirs(path = slides_data_folder, full.names = TRUE, recursive = FALSE)

slides_csv_globs <- slides_folders %>%
  paste0(., "/*.csv")

files <- slides_csv_globs %>% 
  Sys.glob(.) %>% 
  str_subset(pattern="neighbors", negate=T)

# data is already transformed into um, i.e. distance and x,y are in correct units

df <- files %>%
  set_names() %>% 
  map_dfr(function(filename) read_csv(filename, col_types = cols()), .id="source") %>%
  rowwise() %>%
  mutate(dataset=stri_split_fixed(source, "/")[[1]][3], slide=stri_split_fixed(basename(source), ".")[[1]][1]) %>%
  ungroup() %>%
  select(-source) %>%
  mutate(celltype = factor(case_when(celltype %in% c("CD3", "CD4", "CD8", "FOXP3") ~ "CD3",
                                     celltype == "CD20" ~ "CD20",
                                     TRUE ~ "other"))) %>%
  filter(celltype != "other") 

infiltration_df <- df %>%
  mutate(region = ifelse(distance < 0, "tumor", "stroma"), close = distance <100 & distance > -100) %>% 
  filter (close == T) %>%
  group_by(dataset, slide, celltype, region) %>% 
  summarize(count = n())

print(infiltration_df)
write_csv(infiltration_df, paste0(data_folder, "/infiltration.csv"))
  

distance_df <- df %>%
  mutate(close = distance <100 & distance > -100) %>% 
  filter (close == T) %>%
  group_by(dataset, slide, celltype) %>% 
  summarize(
            spread=sd(distance), random_spread=sd(random_distance),
            distance = mean(distance), random_distance=mean(random_distance), 
            count=n())

print(distance_df)
write_csv(distance_df, paste0(data_folder, "/distance_spread.csv"))
  

