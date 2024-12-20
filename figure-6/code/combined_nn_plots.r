rm(list=ls())

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
library(ggpmisc)
library(scales)
library(forcats)

source("code/config.r")
source("../settings.R")

x.spacing <- 0.70
my.theme <- theme(
          strip.background =element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text = element_blank(),
	  strip.clip = "off"
)

load(paste0(data_folder, "/neighbor_workspace.RData"))

df <- df %>% filter(!(celltype %in% c("CD3^'+'~CD45RO^'low'","CD3^'+'~CD45RO^'mid'","CD3^'+'~CD45RO^'hi'")))
df <- df %>% filter(neighbor_type != "CD45RO")

for( cn in c("celltype","neighbor_type") ){
	df[,cn] <- as.character(df[,cn,drop=TRUE])
	df[df[,cn,drop=TRUE]=="CD20^'+'",cn] <- "B"
	df[df[,cn,drop=TRUE]=="FOXP3^'+'",cn] <- "Treg"
	df[df[,cn,drop=TRUE]=="CD8^'+'",cn] <- "CTL"
	df[df[,cn,drop=TRUE]=="CD3^'+'~CD8^'-'",cn] <- "Th"
	df[,cn] <- factor(df[,cn,drop=TRUE])
}

df_blanked <- df %>% filter((dataset == "prostate") & (celltype == "B" | neighbor_type == "B")) 
df <- df %>% filter((dataset != "prostate") | (celltype != "B" & neighbor_type != "B")) 

df$celltype <- ordered(df$celltype, levels=c("Treg", "Th", "CTL", "B"))
df$neighbor_type <- ordered(df$neighbor_type, levels=c("B",  "CTL", "Th", "Treg"))

nearest_neighbor_interactions_p <- df %>% 
  ggplot(aes(y=celltype, x=neighbor_type, fill = log(neighbors))) +
  facet_wrap(~dataset, nrow=1) +
  geom_tile() +
  scale_x_discrete(name="neighbors", guide = guide_axis(angle = -45), 
                   labels = levels(df$neighbor_type), position="top") +
  scale_y_discrete(name="cell type", labels = levels(df$celltype)) +
  geom_text(aes(label=format(round(neighbors, digits=1), nsmall=1)), family = default_font, size=default_pointsize*0.36) +
  geom_text(data=df_blanked, aes(label="-"), color="lightgrey", family = default_font, size=default_pointsize*0.36) +
  scale_fill_distiller(palette="RdBu", limits = c(-1,1)*max(abs(log(df$neighbors))))+
  my.theme +
  theme(legend.position="none",
        strip.placement = "outside",
        strip.text.x = element_text(size=default_pointsize, family = default_font),
        text=element_text(size=default_pointsize, family = default_font),
        axis.title.x.top=element_text(size=default_pointsize, color="black", margin=margin(b=-1.0), family = default_font),
        axis.title.y=element_text(size=default_pointsize, color="black", family = default_font),
        axis.text=element_text(size=default_pointsize, color="#6f6f6f", family = default_font),
        panel.spacing.x = unit(x.spacing, "lines")
  )

pdf_out(paste0(images_folder, "/nearest_neighbor_interactions.pdf"), width=7.5, height=2.3, dpi=72)
print(nearest_neighbor_interactions_p)
dev.off()

load(paste0(data_folder, "/neighbor_per_cell_workspace.RData"))
x.spacing <- 0.85

neighbor_vs_infiltration.all_cells_p <- df  %>% 
    filter(celltype == "FOXP3", neighbor_type=="CD8"  | neighbor_type=="CD4") %>% 
    mutate(location= ifelse(distance < 0, "tumor", "stroma")) %>%
    mutate(location = fct_relevel(location, "tumor", "stroma")) %>%
    ggplot(aes(x=location, y=neighbors, group=neighbor_type, color=neighbor_type)) +
    scale_colour_manual(values = c("#800000", "#1092B8")) +
    ylab("Treg neighbor\npreference ") +
    facet_grid(~dataset) +
    geom_hline(yintercept=1, size=0.5, color="lightgrey", linetype="longdash") +
    stat_summary(fun="mean", geom="line", color="darkgrey") +
  stat_fit_glance(method = "t.test", label.x="left", label.y= c(0.05, 0.30),
                  aes(label = sprintf('italic(p)~"="~\'%1.3f\'', stat(..p.value..))),
                  parse = TRUE, family = default_font, size=default_pointsize * 0.36) +
    stat_summary(fun="mean", geom="point") +
    stat_summary(fun.data = mean_se, geom = "linerange", position=position_dodge(width=0.07)) +
    my.theme +
    theme(legend.position="none",
          panel.spacing.x = unit(x.spacing, "lines"),
          text=element_text(size=default_pointsize, family = default_font),
          axis.title.x=element_blank(),
          axis.text=element_text(size=default_pointsize, family = default_font),
          axis.line = element_line(colour = "black")
          
    )

pdf_out(paste0(images_folder, "/neighbor_vs_infiltration.all_cells.pdf"),
	width=7, height=1.3, dpi=72)
print(neighbor_vs_infiltration.all_cells_p)
dev.off()

neighbor_vs_location.all_cells_p <- df  %>% 
    filter(celltype == "FOXP3", neighbor_type=="CD8") %>% 
    ggplot(aes(x=distance)) +
    facet_grid(~dataset) +
    xlab("distance to tumor") +
    ylab("Tregs") +
    geom_rect(xmin=-Inf, xmax = 0, ymin=-Inf, ymax=Inf, fill="lightgrey") +
    geom_density(fill="#909090", color="#909090") +
    geom_vline(xintercept=0, size=0.5) +
    my.theme +
    theme(legend.position="none",
          text=element_text(size=default_pointsize, family = default_font),
          axis.text.x=element_text(size=default_pointsize),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.spacing.x = unit(x.spacing, "lines"),
          axis.line = element_line(colour = "black")
          
    )

pdf_out(paste0(images_folder, "/neighbor_vs_location.all_cells.pdf"), width=7, height=1.5, dpi=72)
print(neighbor_vs_location.all_cells_p)
dev.off()

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
    geom_point() +
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

# For some reason it does not work without setting a null device 
# https://wilkelab.org/cowplot/reference/set_null_device.html
set_null_device("png")
nn.combined_p <- plot_grid(
	stroma_tumor_counts_p,
	NULL, 
	nearest_neighbor_interactions_p, 
	NULL, 
	neighbor_vs_infiltration.all_cells_p, 
	neighbor_vs_location.all_cells_p, 
	nrow=6, rel_heights=c(3.2,1,2.5, 0.25, 2, 2), align="v", axis="lr")

pdf_out(paste0(images_folder, "/nn.combined.pdf"), width=7.5, height=5.5, dpi=72)
print(nn.combined_p)
dev.off()

