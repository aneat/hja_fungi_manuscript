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
cwm.per <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\cwm.per.needle_24.rds")

##subset guilds
path <- subset(tra, pathogen == 1)
sap <- subset(tra, saprobe == 1)

path_names <- row.names(path)
sap_names <- row.names(sap)

path_otu <- new_otu[path_names]
sap_otu <- new_otu[sap_names]

x_path <- specnumber(path_otu)
y_path <- diversity(path_otu)

x_sap <- specnumber(sap_otu)
y_sap <- diversity(sap_otu)

x <- specnumber(new_otu)
y <- diversity(new_otu)

##match order of matrices
cwm.per <- cwm.per[match(row.names(new_otu), row.names(cwm.per)),]

cwm.per$path.diversity <- y_path
cwm.per$path.richness <- x_path

cwm.per$sap.diversity <- y_sap
cwm.per$sap.richness <- x_sap

cwm.per$total.diversity <- y
cwm.per$total.richness <- x

saveRDS(cwm.per,"\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\needle.div_24.rds")
