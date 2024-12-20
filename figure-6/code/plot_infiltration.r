rm(list=ls())
source("code/config.r")

library(readr)
library(dplyr)
library(ggplot2)
library(png)
library(grid)
library(gridExtra)


create_infiltration_plot <- function(dataset, slide, x_min, x_max, y_min, y_max, x_offset, y_offset) {
    images_name <- paste0(dataset, ".", slide)
    filename <- paste0(images_folder, "/", images_name, ".segmentation.png")
    print(filename)
    m <- readPNG(filename)
    img <- matrix(rgb(m[,,1],m[,,2],m[,,3]), nrow=dim(m)[1]) #0.2 is alpha
    g <- rasterGrob(img, interpolate=TRUE)

    filename <- paste0(images_folder, "/", images_name, ".border.png")
    print(filename)
    m <- readPNG(filename)
    h <- rasterGrob(m, interpolate=TRUE)

    filename <- paste0(slides_data_folder, "/", dataset, "/", slide, ".csv")
    print(filename)
    df <- read_csv(filename, show_col_types=F) %>%
        mutate(celltype = factor(
          case_when(celltype %in% c("CD3", "CD4", "CD8", "FOXP3") ~ "CD3",
                    celltype == "CD20" ~ "CD20",
                    TRUE ~ "other"))
          ) %>%
        filter(celltype != "other")

    x_range <- max(df$x) - min(df$x)
    y_range <- max(df$y) - min(df$y)
    
    bincount <- 300

    df %>% ggplot(aes(x=x, y =-y, group=celltype, color=celltype)) +
        annotation_custom(g, xmin=x_min, xmax=x_max, ymin=-y_max, ymax=y_min) +
        coord_fixed(xlim=c(x_offset,x_offset+2500), ylim=c(-y_offset-2500, -y_offset), expand=F) +
        annotation_custom(h, xmin=x_min, xmax=x_max, ymin=-y_max, ymax=y_min) +
        geom_point(size=0.4, stroke=0.0, shape=16, alpha=1) +
        scale_color_manual(values=c("#E8308A", "#801A1A")) +
        theme_void() + theme(legend.position="none")
        
    ggsave(paste0(images_folder, "/", images_name, ".infiltration.pdf"), width=1.2, height=1.2, useDingbats=F)
    ggsave(paste0(images_folder, "/", images_name, ".infiltration.png"), width=1.2, height=1.2, dpi=600)
}

min_max_df <- read_csv(paste0(data_folder, "/low_high.csv"))
dims_df <- read_csv(paste0(data_folder, "/dims.csv"))
offset_df <- read_csv(paste0(data_folder, "/offset.csv"))

files_df <- inner_join(min_max_df, dims_df)
files_df <- inner_join(files_df, offset_df)
print(files_df)

for (row in 1:nrow(files_df)) {
  dataset <- files_df[row, "dataset"]
  slide <- files_df[row, "slide"]
  x_min <- files_df[row, "x_min"][[1]]
  x_max <- files_df[row, "x_max"][[1]]
  y_min <- files_df[row, "y_min"][[1]]
  y_max <- files_df[row, "y_max"][[1]]
  x_offset <- files_df[row, "x_offset"][[1]]
  y_offset <- files_df[row, "y_offset"][[1]]
  create_infiltration_plot(dataset, slide, x_min, x_max, y_min, y_max, x_offset, y_offset)
}
