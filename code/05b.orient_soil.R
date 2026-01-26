##read in packages
library(ecole)
library(vegan)
library(tidyverse)

##read in otu matrix, taxonomy matrix, and metadata matrix for needle samples.
spe <- read.csv("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\new tests\\otu_24.csv", header = T, row.names=1)
spe_soil <- spe
env <- read.csv("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\meta_data_raw.csv", header = T, row.names=1)
taxa <- read.csv("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\new tests\\taxa_24.csv", header = T, row.names=1)

##look at funguilds and ordination. Obtain this trait matrix from the trait_table_raw.R code
tra <- read.csv("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\trait_table_raw.csv", header = T, row.names = 1)

##remove all foliar samples. These are only soil samples. Also, remove outliers or dysfunctional samples.
env <- subset(env, substrate == "soil")
env <- subset(env, tag != "5129" & tag != "5213")

##Remove samples below 5000 reads
read_count <- as.data.frame(sort(rowSums(spe_soil)))
##low_reads <- c("C11", "C10", "G8", "D2", "E2", "H4") 
##spe_soil <- subset(spe_soil, !(rownames(spe_soil) %in% low_reads))

##match matrix row/column name order to each other
tra <- tra[match(row.names(taxa), row.names(tra)),]
spe_soil <- spe_soil[match(row.names(env), row.names(spe_soil)),]

##Make spe and env matrices match
rownames(spe_soil)
rownames(env)
spe_soil_rows <- rownames(spe_soil)
env_rows <- rownames(env)
common <- intersect(spe_soil_rows, env_rows)
spe_soil <- subset(spe_soil, rownames(spe_soil) %in% common)
env <- subset(env, rownames(env) %in% common)
identical(rownames(spe_soil),
          rownames(env))

##remove species that are not in the dataset
spesoilcounts <- colSums(spe_soil)
spe_soil <- spe_soil[, spesoilcounts > 0]
mx_profile(spe_soil)

spe_soilcounts <- as.data.frame(colSums(spe_soil))
spe_soilcounts$counts <- spe_soilcounts$`colSums(spe_soil)`
spe_soil <- spe_soil[, spe_soilcounts$counts > 0]
spe_soilcounts <- subset(spe_soilcounts, counts >0)

c <- ggplot(spe_soilcounts, aes(x = counts)) + geom_bar() + xlim(0, 50)
c

spe_soil <- spe_soil[, spe_soilcounts$counts > 20]
mx_profile(spe_soil)

##Replace missing value with median values in env matrix
env <- replace(env, env=='', NA)
env <- replace(env, env == "na", NA)
env <- env %>% mutate_at(c(12:44), as.numeric)
env1 <- (env %>% 
           rownames_to_column("row.names") %>%
           group_by(env$species) %>% 
           mutate_at(c(10:45), 
                     ~replace_na(., 
                                 median(., na.rm = TRUE))))
env1 <- data.frame(env1, row.names = 1)

##upload climate data
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
spe_soil_col <- colnames(spe_soil)
taxa_rows <- rownames(taxa)
head(taxa_rows)
head(spe_soil_col)
common <- intersect(spe_soil_col, taxa_rows)
taxa <- subset(taxa, rownames(taxa) %in% common)
identical(rownames(taxa),
          colnames(spe_soil))

##make spe_soil and tra matrices match
spe_soil_col <- colnames(spe_soil)
tra_rows <- rownames(tra)
common <- intersect(spe_soil_col, tra_rows)
tra <- subset(tra, rownames(tra) %in% common)
identical(rownames(tra),
          colnames(spe_soil))

##make spe_soil and tra matrices match
spe_soil_col <- colnames(spe_soil)

##Save filtered matrices
saveRDS(spe_soil, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\spe_soil_24.rds")
saveRDS(env1, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\env1_soil_24.rds")
saveRDS(taxa, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\taxa_soil_24.rds")
saveRDS(tra, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\tra_soil_24.rds")


##THE END##

