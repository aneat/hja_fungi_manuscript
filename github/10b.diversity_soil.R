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
cwm.per <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\cwm.per.soil_24.rds")

##subset guilds
path <- subset(tra, pathogen == 1)
sap <- subset(tra, saprobe == 1)
ecm <- subset(tra, ecm == 1)
am <- subset(tra, am == 1)

path_names <- row.names(path)
sap_names <- row.names(sap)
ecm_names <- row.names(ecm)
am_names <- row.names(am)

path_otu <- new_otu[path_names]
sap_otu <- new_otu[sap_names]
ecm_otu <- new_otu[ecm_names]
am_otu <- new_otu[am_names]

x_path <- specnumber(path_otu)
y_path <- diversity(path_otu)

x_sap <- specnumber(sap_otu)
y_sap <- diversity(sap_otu)

x_ecm <- specnumber(ecm_otu)
y_ecm <- diversity(ecm_otu)

x_am <- specnumber(am_otu)
y_am <- diversity(am_otu)

x <- specnumber(new_otu)
y <- diversity(new_otu)

##match order of matrices
cwm.per <- cwm.per[match(row.names(new_otu), row.names(cwm.per)),]

cwm.per$path.diversity <- y_path
cwm.per$path.richness <- x_path

cwm.per$sap.diversity <- y_sap
cwm.per$sap.richness <- x_sap

cwm.per$ecm.diversity <- y_ecm
cwm.per$ecm.richness <- x_ecm

cwm.per$am.diversity <- y_am
cwm.per$am.richness <- x_am

cwm.per$total.diversity <- y
cwm.per$total.richness <- x

saveRDS(cwm.per, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\soil.diversity_24.rds")
