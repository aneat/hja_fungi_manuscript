##install packages
library(tidyverse)
library(ecole)
library(RColorBrewer)
library(ComplexHeatmap)
library(pheatmap)

##read in new_otu table, env1, and trait matrix
new_otu <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\new_otu_needle_24.rds")
env1 <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\env1_needle_24.rds")
tra <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\tra_needle_24.rds")
taxa <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\taxa_needle_24.rds")

env_low <- subset(env1, elevation < 800)
env_hi <- subset(env1, elevation > 800)

otu_low <- subset(new_otu, rownames(new_otu) %in% rownames(env_low))
otu_low <- otu_low[,colSums(otu_low) > 0]

otu_hi <- subset(new_otu, rownames(new_otu) %in% rownames(env_hi))
otu_hi <- otu_hi[,colSums(otu_hi) > 0]

env_low <- env_low[match(row.names(otu_low), row.names(env_low)),]
env_hi <- env_hi[match(row.names(otu_hi), row.names(env_hi)),]


##save csv files##
write.csv(env_low, "\\Users\\abbey\\Desktop\\needles_env_low_24.csv")
write.csv(env_hi, "\\Users\\abbey\\Desktop\\needles_env_hi_24.csv")
write.csv(otu_low, "\\Users\\abbey\\Desktop\\needles_otu_low_24.csv")
write.csv(otu_hi, "\\Users\\abbey\\Desktop\\needles_otu_hi_24.csv")

##save RDS files##
saveRDS(env_low, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\needles_env_low_24.rds")
saveRDS(env_hi, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\needles_env_hi_24.rds")
saveRDS(otu_low, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\needles_otu_low_24.rds")
saveRDS(otu_hi, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\needles_otu_hi_24.rds")
