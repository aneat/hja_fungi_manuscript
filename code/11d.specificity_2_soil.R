##load packages
library(tidyverse)
library(indicspecies)
library(purrr)
library(ecole)
library(ggpubr)

new_otu <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\new_otu_soil_24.rds")
env1 <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\env1_soil_24.rds")
taxa <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\taxa_soil_24.rds")
tra <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\tra_soil_24.rds")

env_low <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\soils_env_low_24.rds")
env_hi <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\soils_env_hi_24.rds")
otu_low <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\soils_otu_low_24.rds")
otu_hi <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\soils_otu_hi_24.rds")

isa_low <- read.csv("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\isa_soils_low_24.csv")
isa_hi <- read.csv("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\isa_soils_hi_24.csv")

##LOW##
isa_low <- (isa_low %>% 
  mutate(across(
    everything(),
    ~ map_chr(.x, ~ gsub("\"", "", .x))
  )))
otu_low_og <- otu_low
otu_low[otu_low > 0] = 1
occur <- as.data.frame(colSums(otu_low))

isa_low <- column_to_rownames(isa_low, var = "out.id")
isa_low$p <- as.numeric(isa_low$p)
isa_low$value.IV <- as.numeric(isa_low$value.IV)
isa_low_merge <- merge(isa_low, occur, by = 0)
isa_low_merge$occurances <- isa_low_merge$`colSums(otu_low)`
isa_low_merge <- subset(isa_low_merge, select = -c(7))

p_low <- ggplot(isa_low_merge, aes(x = occurances, y = p)) + geom_point(aes(color = value.IV)) + geom_vline(xintercept = 4) + geom_hline(yintercept = 0.05) + scale_x_continuous(breaks = seq(0, 60, by=5))
p_low
isa_low_sub <- subset(isa_low_merge, occurances > 3)
isa_low_sub$index <- 1 - isa_low_sub$p
saveRDS(isa_low_sub, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\isa_low_sub_soils.rds" )

##perform cwm analysis
isa_low_list <- isa_low_sub$Row.names
isa_low_sub <- rownames_to_column(isa_low_sub, "numbers")
isa_low_sub <- column_to_rownames(isa_low_sub, "Row.names")
otu_isa_low <- otu_low_og[, (colnames(otu_low_og) %in% isa_low_list)]
isa_low_cwm <- makecwm(otu_isa_low, isa_low_sub$index)


##HIGH##
isa_hi <- (isa_hi %>% 
              mutate(across(
                everything(),
                ~ map_chr(.x, ~ gsub("\"", "", .x))
              )))

otu_hi[otu_hi > 0] = 1
occur <- as.data.frame(colSums(otu_hi))

isa_hi <- column_to_rownames(isa_hi, var = "out.id")
isa_hi$p <- as.numeric(isa_hi$p)
isa_hi$value.IV <- as.numeric(isa_hi$value.IV)
isa_hi_merge <- merge(isa_hi, occur, by = 0)
isa_hi_merge$occurances <- isa_hi_merge$`colSums(otu_hi)`
isa_hi_merge <- subset(isa_hi_merge, select = -c(7))

p_hi <- ggplot(isa_hi_merge, aes(x = occurances, y = p)) + geom_point(aes(color = value.IV)) + geom_vline(xintercept = 4) + geom_hline(yintercept = 0.05) + scale_x_continuous(breaks = seq(0, 60, by=5))
p_hi

isa_hi_sub <- subset(isa_hi_merge, occurances > 3)
isa_hi_sub$index <- 1 - isa_hi_sub$p

##perform cwm analysis
isa_hi_list <- isa_hi_sub$Row.names
isa_hi_sub <- rownames_to_column(isa_hi_sub, "numbers")
isa_hi_sub <- column_to_rownames(isa_hi_sub, "Row.names")
otu_isa_hi <- otu_hi[, (colnames(otu_hi) %in% isa_hi_list)]
isa_hi_cwm <- makecwm(otu_isa_hi, isa_hi_sub$index)

cwm_merge <- rbind(isa_low_cwm, isa_hi_cwm)

env1 <- env1[match(row.names(cwm_merge), row.names(env1)),]
env1$gen.host.spe <- cwm_merge$V1

gs <- ggplot(env1, aes(x= elevation, y = gen.host.spe)) + geom_point(size = 2) + ylim(.40,.90) + ylab("specificity cwm - all soil taxa") + geom_smooth(method = "lm", color = "#993404", size = 3) + ylim(.4,.9)
gs

##pathogen##
##LOW##
pathogen <- subset(tra, pathogen == "1")
pathogen_names <- rownames(pathogen)
isa_path_low <- subset(isa_low_sub, rownames(isa_low_sub) %in% pathogen_names)
otu_path_low <- otu_low_og[, (colnames(otu_low_og) %in% rownames(isa_path_low))]

path_cwm_low <- makecwm(otu_path_low, isa_path_low$index)

##HIGH##
pathogen <- subset(tra, pathogen == "1")
pathogen_names <- rownames(pathogen)
isa_path_hi <- subset(isa_hi_sub, rownames(isa_hi_sub) %in% pathogen_names)
otu_path_hi <- otu_hi[, (colnames(otu_hi) %in% rownames(isa_path_hi))]

path_cwm_hi <- makecwm(otu_path_hi, isa_path_hi$index)

cwm_merge <- rbind(path_cwm_low, path_cwm_hi)

env1 <- env1[match(row.names(cwm_merge), row.names(env1)),]
env1$path.host.spe <- cwm_merge$V1

ps <- ggplot(env1, aes(x= elevation, y = path.host.spe)) + geom_point(size = 2) + ylab("specificity cwm - soil pathogens") + geom_smooth(method = "lm", color = "#993404", size = 3)
ps


##ecm##
##LOW##
ecm <- subset(tra, ecm == "1")
ecm_names <- rownames(ecm)
isa_ecm_low <- subset(isa_low_sub, rownames(isa_low_sub) %in% ecm_names)
otu_ecm_low <- otu_low_og[, (colnames(otu_low_og) %in% rownames(isa_ecm_low))]

ecm_cwm_low <- makecwm(otu_ecm_low, isa_ecm_low$index)

##HIGH##
ecm <- subset(tra, ecm == "1")
ecm_names <- rownames(ecm)
isa_ecm_hi <- subset(isa_hi_sub, rownames(isa_hi_sub) %in% ecm_names)
otu_ecm_hi <- otu_hi[, (colnames(otu_hi) %in% rownames(isa_ecm_hi))]

ecm_cwm_hi <- makecwm(otu_ecm_hi, isa_ecm_hi$index)

cwm_merge <- rbind(ecm_cwm_low, ecm_cwm_hi)

env1 <- env1[match(row.names(cwm_merge), row.names(env1)),]
env1$ecm.host.spe <- cwm_merge$V1

es <- ggplot(env1, aes(x= elevation, y = ecm.host.spe)) + geom_point(size = 2) + ylab("specificity cwm - ecm") + geom_smooth(method = "lm", color = "#993404", size = 3) + ylim(.4,.9)
es


##am##
##LOW##
am <- subset(tra, am == "1")
am_names <- rownames(am)
isa_am_low <- subset(isa_low_sub, rownames(isa_low_sub) %in% am_names)
otu_am_low <- otu_low_og[, (colnames(otu_low_og) %in% rownames(isa_am_low))]

am_cwm_low <- makecwm(otu_am_low, isa_am_low$index)

##HIGH##
am <- subset(tra, am == "1")
am_names <- rownames(am)
isa_am_hi <- subset(isa_hi_sub, rownames(isa_hi_sub) %in% am_names)
otu_am_hi <- otu_hi[, (colnames(otu_hi) %in% rownames(isa_am_hi))]

am_cwm_hi <- makecwm(otu_am_hi, isa_am_hi$index)

cwm_merge <- rbind(am_cwm_low, am_cwm_hi)

env1 <- env1[match(row.names(cwm_merge), row.names(env1)),]
env1$am.host.spe <- cwm_merge$V1

as <- ggplot(env1, aes(x= elevation, y = am.host.spe)) + geom_point(size = 2) + ylab("specificity cwm - am") + geom_smooth(method = "lm", color = "#993404", size = 3) + ylim(.4,.9)
as

##saprobe##
##LOW##
sap <- subset(tra, saprobe == "1")
sap_names <- rownames(sap)
isa_sap_low <- subset(isa_low_sub, rownames(isa_low_sub) %in% sap_names)
otu_sap_low <- otu_low_og[, (colnames(otu_low_og) %in% rownames(isa_sap_low))]

sap_cwm_low <- makecwm(otu_sap_low, isa_sap_low$index)

##HIGH##
sap <- subset(tra, sap == "1")
sap_names <- rownames(sap)
isa_sap_hi <- subset(isa_hi_sub, rownames(isa_hi_sub) %in% sap_names)
otu_sap_hi <- otu_hi[, (colnames(otu_hi) %in% rownames(isa_sap_hi))]

sap_cwm_hi <- makecwm(otu_sap_hi, isa_sap_hi$index)

cwm_merge <- rbind(sap_cwm_low, sap_cwm_hi)

env1 <- env1[match(row.names(cwm_merge), row.names(env1)),]
env1$sap.host.spe <- cwm_merge$V1

ss <- ggplot(env1, aes(x= elevation, y = sap.host.spe)) + geom_point(size = 2) + ylab("specificity cwm - soil saprobes") + geom_smooth(method = "lm", color = "#993404", size = 3) + ylim(.4,.9)
ss

##combine all into one figure
env_guild_cwm <- subset(env1, select = c(54:58))
colnames(env_guild_cwm)[1] <- "all taxa"
colnames(env_guild_cwm)[2] <- "pathogens"
colnames(env_guild_cwm)[3] <- "ecm"
colnames(env_guild_cwm)[4] <- "am"
colnames(env_guild_cwm)[5] <- "saprobes"

egc_pivot <- pivot_longer(env_guild_cwm, col = c(1:4), names_to = "guild")
egc_pivot$guild <- factor(egc_pivot$guild,levels= c("all taxa", "pathogens", "ecm", "am"))
box_soil <- ggplot(egc_pivot, aes(x= guild, y = value)) + geom_boxplot() + geom_jitter(color = "#993404", size = 1.5) + ylim(.4, 1.00) + ylab("specificity index")
box_soil

ggarrange(box_soil,box_needle, ncol = 2)
ggarrange(gs, ps, es, as, g, p, ncol = 4, nrow = 2)

##save soil host specificifty RDS##
env_spec <- subset(env1, select = c(54:58))
saveRDS(env_spec,"\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\soil.spec_24.rds")

##THE END##

