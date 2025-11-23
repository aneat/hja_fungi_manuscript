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
env1 <-  readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\env1_needle.rds")
tra <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\tra_needle.rds")
dist <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\dist_needle_24.rds")
taxa <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\taxa_needle.rds")
per.table <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\per.table.needle.rds")
isa_join <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\isa_join.rds")



##ABUNDANCE TABLE##
env1$elevation_band <- cut(env1$elevation, c(0,600,910, Inf), c("low", "medium", "high"))
env1 <- subset(env1, select = c("species", "elevation_band"))
otu_notho <- subset(new_otu, select = c("otu.71"))

env1 <- env1[match(row.names(otu_notho), row.names(env1)),]
final <- cbind(env1, otu_notho)

summary <- final %>%
  group_by(species, elevation_band) %>%
  summarise(mean = mean(otu.71))
summary$prop <- (summary$mean/sum(summary$mean))

summary <- subset(summary, select = -c(3))
summary$prop <- round(summary$prop, digits = 2)
summary_pivot <- pivot_wider(summary, names_from = elevation_band, values_from = prop)
summary_pivot <- column_to_rownames(summary_pivot, "species")

my_palette_needle <- brewer.pal(name="Greens",n=9)[2:8]
gg_heat <- ggplot(summary, aes(x = species, y = elevation_band, fill = prop)) + geom_tile( color = "black", size = 1.3) + scale_fill_gradient(name = my_palette_needle, low = "#C7E9C0", high = "#006D2C") + geom_text(aes(label= prop), size = 6) + xlab("tree host species") + ylab("site elevation") + theme_classic() + guides(fill = "none") + theme(text=element_text(size=20))
gg_heat
##notho_heat <- ComplexHeatmap::pheatmap(summary_pivot, color = my_palette_needle, treeheight_col = 0, treeheight_row = 0, cluster_rows = FALSE, cluster_cols = FALSE, legend = FALSE, display_numbers = TRUE, number_color = "#737373", fontsize = 15)
##notho_heat

##COMMON SPECIES TABLE##
common_sp <- as.data.frame(colSums(new_otu))
common_sp$read_count <- common_sp$`colSums(new_otu)`
common_sp$read_count <- round(common_sp$read_count, digits = 0)
common_sp <- arrange(common_sp, desc(read_count))
common_sp$n.gaeumanii <- "no"
common_sp$n.gaeumanii[24:28] <- "yes"
hist1 <- ggplot(common_sp, aes(x = read_count) ) + geom_histogram(bins = 30, color = "black", aes(fill = n.gaeumanii), size = 1.3) + scale_x_log10() + scale_fill_manual(values = c("grey","#006D2C")) + ylab("# of otus") + xlab("log(total rarefied read count)") + theme_classic() + guides(fill = "none") + theme(text=element_text(size=20)) + geom_segment(aes(x = 13000, y = 30, xend = 13000, yend = 9), arrow = arrow(), size = 1.5)
hist1

##SPECIFICITY INDEX TABLE##
isa_join$avg <- round(isa_join$avg, digits = 3)
isa_join <- arrange(isa_join, desc(avg))

isa_join$n.gaeumanii <- "no"
isa_join$n.gaeumanii[76:116] <- "yes"
hist2 <- ggplot(isa_join, aes(x = avg) ) + geom_histogram(bins = 30, color = "black", aes(fill = n.gaeumanii), size = 1.3) + scale_fill_manual(values = c("grey","#006D2C")) + ylab("# of otus") + xlab("taxon level specificity index") + theme_classic() + guides(fill = "none") + theme(text=element_text(size=20))
#+ scale_x_log10() + scale_fill_manual(values = c("#41AB5D","#00441B")) + annotate(geom = "text", x = 22000, y = 15, label = "N. gaeumanii", color = "#00441B", angle = 45) + ylab("# of otus") + xlab("log(total rarefied read count)")
hist2

notho <- ggarrange(hist1, hist2, gg_heat, nrow = 1)
ggsave("notho.png", notho, units = "px", dpi = 450, width = 8000, height = 3000)


##NOTHO Regression
merge <- merge(otu_notho, per.table, by = "row.names", all = TRUE)
merge

merge_sub <- subset(merge, species == "psme")

notho_lm <- lm(otu.71 ~ climate, data = merge_sub)

summary(notho_lm)
