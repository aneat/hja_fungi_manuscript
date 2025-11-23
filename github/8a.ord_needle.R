##read in packages
library(tidyverse)
library(vegan)
library(ecole)
library(phyloseq)
library(RColorBrewer)
library(ggpubr)
library(gridExtra)

##read in new_otu table, env1, and trait matrix
new_otu <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\new_otu_needle_24.rds")
env1 <-  readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\env1_needle_24.rds")
tra <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\tra_needle_24.rds")
dist <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\dist_needle_24.rds")
taxa <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\taxa_needle_24.rds")
per.table <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\per.table.needle_24.rds")
per.table$species <- gsub("psme", "p. menziesii", per.table$species)
per.table$species <- gsub("tshe", "t. heterophylla", per.table$species)
per.table$species <- gsub("tabr", "t. brevifolia", per.table$species)

##ordinate
nmds <- metaMDS(dist)
nmds
m0 <- MDSrotate(nmds, as.numeric(per.table$climate))
ecole:::fitstats_nms(m0)
m0
plot(m0, "sites")
env_sub <- subset(per.table, select =  c(8))
e <- envfit(m0, env_sub, na.rm = T)
e
plot(e)

##Create a df in ggplot with NMS scores and environmental variables
m_ggplot <- as.data.frame(vegan::scores(m0, display = c("sites")))

##GGplot matrix with environmental scores
env.scores <- as.data.frame(scores(e, display = "vectors"))
head(env.scores)
env.scores <- cbind(env.scores, env.variables = rownames(env.scores))
env.scores <- cbind(env.scores, pval = e$vectors$pvals)
env.scores

##isolate only the significant environmental scores
env.scores <- subset(env.scores, pval <= .1)
head(env.scores)

##GGplot base plot
nms_ggplot <- ggplot(m_ggplot, aes(x=NMDS1, y=NMDS2)) + geom_point(aes(NMDS1, NMDS2, fill = per.table$species, size = per.table$climate, stroke =2), color = "black", pch = 21) + scale_size_continuous(range = c(2, 15)) + scale_alpha(guide="none") + stat_ellipse(geom = "polygon", aes(fill = per.table$species), alpha = 0.25) + coord_equal()
nms_ggplot
new_plot <- (nms_ggplot + scale_fill_brewer(palette = "Dark2") + theme_classic(base_size=15) + labs( fill = "tree host species", size = "climate PCA axis 1"))
new_plot
new_plot_1 <- (new_plot + geom_segment(data = env.scores, aes(x = 0, xend=NMDS1*1.5, y=0, yend=NMDS2*1.5), arrow = arrow(length = unit(0.6, "cm")), colour = "grey10", lwd=1.5) + ggrepel::geom_label_repel(data = env.scores, aes(x=NMDS1, y=NMDS2, label = env.variables), size = 10, cex = 4, segment.size = 0.5, nudge_x = 1))
new_plot_1
common_breaks <- seq(-1, 2, by = 0.5)
needle_plot <- (new_plot_1 + guides(fill = guide_legend(override.aes = list(size = 5)))) + coord_fixed(ratio = 1) + theme_light(base_size = 30) + scale_x_continuous(limits = c(-1,1.75), breaks = common_breaks)
needle_plot
ggsave("needle_plot_climate.png", needle_plot, units = "px", width = 5200, height = 3000, dpi = 350)
saveRDS(needle_plot,"\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\needle_plot_24.rds")
