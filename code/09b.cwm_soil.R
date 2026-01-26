##read in packages
library(tidyverse)
library(vegan)
library(ecole)
library(RColorBrewer)
library(ComplexHeatmap)
library(MuMIn)
library(AICcmodavg)

##read in new_otu table, env1, and trait matrix
new_otu <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\new_otu_soil_24.rds")
env1 <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\env1_soil_24.rds")
tra <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\tra_soil_24.rds")
per.table <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\per.table.soil_24.rds")

##perform cwm##
cwm <- makecwm(new_otu, tra, na.rm = TRUE)
cwm <- subset(cwm)
head(cwm)

##combine per.table and cwm
cwm.per <- merge(cwm, per.table, by = 0)
cwm.per <- column_to_rownames(cwm.per, "Row.names")
env1 <- subset(env1, select = -c(4,6,10,34))
cwm.per <- merge(cwm.per, env1, by = 0)
cwm.per <- column_to_rownames(cwm.per, "Row.names")
saveRDS(cwm.per, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\cwm.per.soil_24.rds")

##find counts associated with each guild
guild_counts_soil <- colSums(tra, na.rm = FALSE, dims = 1)

guild_counts_soil <- as.data.frame(guild_counts_soil)
write.csv(guild_counts_soil, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\guild_counts_soil_24.csv")
