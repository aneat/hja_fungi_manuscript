##read in packages
library(tidyverse)
library(vegan)
library(ecole)
library(phyloseq)
library(RColorBrewer)
library(AICcmodavg)
library(MuMIn)

##read in new_otu table, env1, and trait matrix
new_otu <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\new_otu_needle_24.rds")
spe_needle <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\spe_needle_24.rds")
env1 <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\env1_needle_24.rds")
tra <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\tra_needle_24.rds")
taxa <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\taxa_needle_24.rds")
per.table <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\per.table.needle_24.rds")

##perform cwm##
cwm <- makecwm(new_otu, tra, na.rm = TRUE)
cwm <- subset(cwm, select = -c(8:9, 11:12,18:19))
head(cwm)

##combine per.table and cwm
cwm.per <- merge(cwm, per.table, by = 0)
cwm.per <- column_to_rownames(cwm.per, "Row.names")
env1 <- subset(env1, select = -c(4,6,10,34))
cwm.per <- merge(cwm.per, env1, by = 0)
cwm.per <- column_to_rownames(cwm.per, "Row.names")
saveRDS(cwm.per, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\cwm.per.needle_24.rds")

##find counts associated with each guild
guild_counts_needle <- colSums(tra, na.rm = FALSE, dims = 1)

guild_counts_needle <- as.data.frame(guild_counts_needle)
write.csv(guild_counts_needle, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\guild_counts_needle_24.csv")
