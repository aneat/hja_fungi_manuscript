  ##read in packages
  library(tidyverse)
  library(vegan)
  library(ecole)
  library(RColorBrewer)
  library(ComplexHeatmap)

  ##read in new_otu table, env1, and trait matrix
  new_otu <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\new_otu_soil_24.rds")
  env1 <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\env1_soil_24.rds")
  tra <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\tra_soil_24.rds")
  dist <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\dist_soil_24.rds")
  per.table.soil <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\per.table.soil_24.rds")
  per.table.soil$species <- gsub("psme", "p. menziesii", per.table.soil$species)
  per.table.soil$species <- gsub("tshe", "t. heterophylla", per.table.soil$species)
  per.table.soil$species <- gsub("tabr", "t. brevifolia", per.table.soil$species)
  ##ordinate
  nmds <- metaMDS(dist)
  m0_soil <- MDSrotate(nmds, as.numeric(per.table.soil$climate))
  plot(m0_soil, "sites")
  m0_soil
  plot(m0_soil, "sites")

  ##cumulative variance explained is metric fit, R2m##
  # metric fit = squared correlation between original dissimilarities and ordination distances. Foreign to NMS.  Null: no linear relationship between dissimilarities and ordination distances.
  ecole:::fitstats_nms(m0_soil)

  env_sub <- subset(per.table.soil, select =  c(2,8))
  es <- envfit(m0_soil, env_sub, na.rm = T)
  plot(es)
  es
  ##Create a df in ggplot with NMS scores and environmental variables
  m_ggplot_soil <- as.data.frame(vegan::scores(m0_soil, display = "sites"))

  ##GGplot matrix with environmental scores
  env.scores.soil <- as.data.frame(vegan::scores(es, display = "vectors"))
  head(env.scores.soil)
  env.scores.soil <- cbind(env.scores.soil, env.variables = rownames(env.scores.soil))
  env.scores.soil <- cbind(env.scores.soil, pval = es$vectors$pvals)
  env.scores.soil

  ##isolate only the significant environmental scores
  env.scores.soil <- subset(env.scores.soil, env.variables == "climate")
  env.scores.soil

  ##GGplot base plot
  nms_ggplot_soil <- ggplot(m_ggplot_soil, aes(x=NMDS1, y=NMDS2)) + geom_point(aes(NMDS1, NMDS2, fill = per.table.soil$species, size = per.table.soil$climate, stroke =2), color = "black", pch = 21) +  scale_size_continuous(range = c(2, 15)) + scale_alpha(guide="none") + stat_ellipse(geom = "polygon", aes(fill = per.table.soil$species), alpha = 0.25) + coord_equal()

  nms_ggplot_soil
  new_plot_soil <- (nms_ggplot_soil + scale_fill_brewer(palette = "Dark2") + theme_classic(base_size=15) + labs( fill = "tree host species", size = "climate PCA axis 1"))
  new_plot_soil
  new_plot_1_soil <- (new_plot_soil + geom_segment(data = env.scores.soil, aes(x = 0, xend=NMDS1*1.5, y=0, yend=NMDS2*1.5), arrow = arrow(length = unit(0.6, "cm")), colour = "grey10", lwd=1.5) + ggrepel::geom_label_repel(data = env.scores.soil, aes(x=NMDS1*.4, y=NMDS2*.4, label = env.variables), size = 10, cex = 4, segment.size = 0.5))
  new_plot_1_soil
  soil_plot <- (new_plot_1_soil) + coord_fixed(ratio = 1) + theme_light(base_size = 30) + ylim(-1,1.0)
  soil_plot
  ggsave("soil_plot_climate.png", soil_plot, units = "px", width = 5000, height = 3000, dpi = 350)

