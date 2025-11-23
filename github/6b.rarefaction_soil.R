##read-in packages
library(tidyverse)
library(vegan)
library(ggpubr)
library(RColorBrewer)

spe_soil <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\spe_soil_24.rds")
spe_needle <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\spe_needle_24.rds")
env_s <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\env1_soil_24.rds")
soil_rare <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\soilrare_24.rds")
needle_rare <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\needlerare_24.rds")
env_s_sub <- subset(env_s, select = c(4,6,10))
env_s_sub <- rownames_to_column(env_s_sub, var = "Group")


##create a function to rarefy otu table and produce a table of each site, otu, and abundance after one rarefy iteration
rarefy <- function(x){
  rare_table <- rrarefy(x, raremax)
  rare_tibble <- (as_tibble(rare_table, rownames = NA) %>%
                    rownames_to_column(var = "Group") %>%
                    pivot_longer(-Group))
  return(rare_tibble)
}

##test function and observe output
(raremax <- min(rowSums(spe_needle)))
rarefy(spe_soil)
raremax
##iterate function to go through 1000 passes of rarefying the otu table
rarefy_iterations <- map_dfr(1:1000, ~rarefy(spe_soil), .id = "iteration")

##create table with mean values for each site, otu, abundance combination
##detach(package:Hmisc)
summary_iterations <- (rarefy_iterations %>%
                         group_by(Group, name) %>%
                         summarize(mean_abund = mean(value),
                                   .groups = "drop"))
##reformat the table so it becomes a species x site matrix again
new_otu <- (summary_iterations %>%
              pivot_wider(names_from = name, values_from = mean_abund))
new_otu <- (new_otu %>%
              column_to_rownames(var = "Group"))

##create rarefaction curves to visualize the species accumulation.
rarecurvedata <- rarecurve(spe_soil, step = 100)
map_dfr(rarecurvedata, bind_rows) %>%
  bind_cols(Group = rownames(spe_soil),.) %>%
  pivot_longer(-Group) %>%
  drop_na() %>%
  mutate(n_seqs = as.numeric(str_replace(name, "N", ""))) %>%
  select(-name) %>%
  ggplot(aes(x=n_seqs, y=value, group = Group)) + geom_vline(xintercept = raremax, color = "gray") + geom_line()

soil_table <- (map_dfr(rarecurvedata, bind_rows) %>%
  bind_cols(Group = rownames(spe_soil),.) %>%
  pivot_longer(-Group) %>%
  drop_na() %>%
  mutate(n_seqs = as.numeric(str_replace(name, "N", ""))) %>%
  select(-name))

table_join <- left_join(env_s_sub, soil_table, by = "Group")
my_palette_soil <- brewer.pal(name="YlOrBr",n=9)[4:9]
soil_rare <- (ggplot(table_join, aes(x=n_seqs, y=value, group = Group, color = elevation)) + geom_vline(xintercept = raremax, color = "gray", linewidth = 1.2) + geom_line(linewidth = 1.2) + labs(y = "# OTUs per sample", x= "# of sequences") + theme_classic(base_size=20) + scale_color_gradient(name = my_palette_soil, low = "#FEC44F", high = "#662506", breaks=c(500, 750, 1000), labels = c("500", "750","1000")) + guides(color=guide_colorbar(title="elevation (m)")) + theme(legend.position = c(0.8, 0.25)) + scale_y_continuous(limits = c(0, 700), breaks = seq(0, 700, by = 100)) +scale_x_continuous(limits = c(0, 150000), breaks = seq(0, 150000, by = 50000)))
soil_rare
saveRDS(soil_rare, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\soilrare_24.rds")
saveRDS(needle_rare, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\needlerare_24.rds")

rare_fig <- ggarrange(soil_rare, needle_rare, ncol=2, align = "hv")
ggsave("sup_rare.png", rare_fig, units = "px", width = 8000, height = 4200, dpi = 550)



##examine one iteration of rarefied otu table and compare it to function output. This checks that the function worked properly. The values should be similar.
rare <- rrarefy(spe_soil, raremax)
rare <- as.data.frame(rare)
new_otu<-new_otu[names(rare)]

rare <- rare[match(row.names(new_otu), row.names(rare)),]

##make distance matrix
dist <- vegdist(new_otu, dmethod= "bray")
saveRDS(dist, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\dist_soil_24.rds")

##save new_otu table to use for analysis
saveRDS(new_otu, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\new_otu_soil_24.rds")
saveRDS(table_join, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\table_join_soil_24.RDS")
table_join <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\table_join_soil_24.RDS")

##THE END##
