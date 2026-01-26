##read in packages
library(ecole)
library(vegan)
library(tidyverse)

##read in otu matrix, taxonomy matrix, and metadata matrix for needle samples.
spe <- read.csv("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\new tests\\otu_24.csv", header = T, row.names=1)
spe_needle <- spe
env <- read.csv("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\meta_data_raw.csv", header = T, row.names=1)
taxa <- read.csv("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\new tests\\taxa_24.csv", header = T, row.names=1)

##look at funguilds and ordination. Obtain this trait matrix from the trait_table_raw.R code
tra <- read.csv("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\trait_table_raw.csv", header = T, row.names = 1)

##remove all soil samples. These are only needle samples. Also, remove outliers or dysfunctional samples.
env <- subset(env, substrate == "needle")
env <- subset(env, tag != "5268" & tag != "4472")

##match matrix row/column name order to each other
tra <- tra[match(row.names(taxa), row.names(tra)),]
spe_needle <- spe_needle[match(row.names(env), row.names(spe_needle)),]

##Make spe and env matrices match
rownames(spe_needle)
rownames(env)
spe_needle_rows <- rownames(spe_needle)
env_rows <- rownames(env)
common <- intersect(spe_needle_rows, env_rows)
spe_needle <- subset(spe_needle, rownames(spe_needle) %in% common)
env <- subset(env, rownames(env) %in% common)
identical(rownames(spe_needle),
          rownames(env))

##remove species that are not in the dataset
spe_needlecounts <- as.data.frame(colSums(spe_needle))
spe_needlecounts$counts <- spe_needlecounts$`colSums(spe_needle)`
spe_needle <- spe_needle[, spe_needlecounts$counts > 0]
spe_needlecounts <- subset(spe_needlecounts, counts >0)

c <- ggplot(spe_needlecounts, aes(x = counts)) + geom_bar() + xlim(0, 50)
c

spe_needle <- spe_needle[, spe_needlecounts$counts > 20]
mx_profile(spe_needle)

##Replace missing value with median values in env matrix
env <- replace(env, env=='', NA)
env <- replace(env, env == "na", NA)
env <- env %>% mutate_at(c(13:44), as.numeric)
env1 <- (env %>% 
           rownames_to_column("row.names") %>%
           group_by(env$species) %>% 
           mutate_at(c(10:45), 
                     ~replace_na(., 
                                 median(., na.rm = TRUE, row.names = TRUE))))
env1 <- data.frame(env1, row.names = 1)

##add climate data to env1
load("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\RefStand_Diversity_info_for_Abbey_wMicroclimateData_20240709.RData")
climate.data <- div.data6
env1$site <- toupper(env1$site)
env1$STANDID <- env1$site
env1$row.names <- rownames(env1)
env1 <- left_join(env1, climate.data, by = "STANDID")
env1 <- column_to_rownames(env1, "row.names")
env1 <- subset(env1, select = -c(45:48, 51:53))
env1$row.names <- rownames(env1)

##make spe and taxa matrices match
spe_needle_col <- colnames(spe_needle)
taxa_rows <- rownames(taxa)
head(taxa_rows)
head(spe_needle_col)
common <- intersect(spe_needle_col, taxa_rows)
taxa <- subset(taxa, rownames(taxa) %in% common)
identical(rownames(taxa),
          colnames(spe_needle))

##make spe and tra matrices match
spe_needle_col <- colnames(spe_needle)
tra_rows <- rownames(tra)
common <- intersect(spe_needle_col, tra_rows)
tra <- subset(tra, rownames(tra) %in% common)
identical(rownames(tra),
          colnames(spe_needle))

##Save filtered matrices
saveRDS(spe_needle, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\spe_needle_24.rds")
saveRDS(env1, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\env1_needle_24.rds")
saveRDS(taxa, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\taxa_needle_24.rds")
saveRDS(tra, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\tra_needle_24.rds")


##THE END##

