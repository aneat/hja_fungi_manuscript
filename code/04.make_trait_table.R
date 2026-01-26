##read in packages
library("tidyverse")
library("stringr")

##read in fungaltraits datasheet
mega <- read.csv("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\FungalTraits 1.2_ver_16Dec_2020 - V.1.2.csv", header = T)

##This was made on R Studio Cloud with Funguild!!##
taxa_table <- read.csv("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\new tests\\guilds_24.csv", header = T)

#rename Genus column so that genus column in taxa matches Genus column in genus_group
mega <-  (mega %>% 
  rename(
    Genus = GENUS,
  ))

##join taxa and genus_group so that every otu has one row
join <- left_join(taxa_table, mega, 
           by = c("Genus"))
join<- (join %>%
  column_to_rownames(var="X"))

##make a new column in the data table with merged guild assignments. The first pass will come from the calls from the Tedersoo data table (mega). The second pass will come from Funguild (taxa).
join1 <- (join %>% mutate(guild_merge = coalesce(primary_lifestyle, guild)))

##edit nomenclature so that info from the two tables are consistent in this new table. This also replaces blank spaces witn NA
join1$guild_merge <- tolower(join1$guild_merge)
join1$guild_merge <- gsub("arbuscular_mycorrhizal", "arbuscular mycorrhizal", join1$guild_merge)
join1$guild_merge <- gsub("plant_pathogen", "plant pathogen", join1$guild_merge) 
join1[join1==""]<-NA

##create a new table with the traits we are interested in. If the trait is present for the otu, then it will be assigned a value of "1". If it is absent, it will be assigned a value of "0". If we are not sure, it remains as NA. 
subset <- join1[c(27,28,29,30,34)]
is.na(subset$ecm) <- is.na(subset$guild_merge)
subset$ecm <- ifelse(grepl("ectomycorrhizal", subset$guild_merge), "1", "0")
is.na(subset$ecm) <- is.na(subset$guild_merge)
subset$pathogen <- ifelse(grepl("plant pathogen", subset$guild_merge), "1", "0")
is.na(subset$pathogen) <- is.na(subset$guild_merge)
subset$saprobe <- ifelse(grepl("saprotroph", subset$guild_merge), "1", "0")
is.na(subset$saprobe) <- is.na(subset$guild_merge)
subset$wood.saprobe <- ifelse(grepl("wood_saprotroph", subset$guild_merge), "1", "0")
is.na(subset$wood.saprobe) <- is.na(subset$guild_merge)
subset$soil.saprobe <- ifelse(grepl("soil_saprotroph", subset$guild_merge), "1", "0")
is.na(subset$soil.saprobe) <- is.na(subset$guild_merge)
subset$litter.saprobe <- ifelse(grepl("litter_saprotroph", subset$guild_merge), "1", "0")
is.na(subset$litter.saprobe) <- is.na(subset$guild_merge)
subset$lichenized <- ifelse(grepl("lichenized", subset$guild_merge), "1", "0")
is.na(subset$lichenized) <- is.na(subset$guild_merge)
subset$am <- ifelse(grepl("arbuscular mycorrhizal", subset$guild_merge), "1", "0")
is.na(subset$am) <- is.na(subset$guild_merge)
subset$ecm.short.distance <- ifelse(grepl("short-distance", subset$Ectomycorrhiza_exploration_type_template), "1", "0")
is.na(subset$ecm.short.distance) <- is.na(subset$Ectomycorrhiza_exploration_typ_template)
subset$ecm.medium.distance <- ifelse(grepl("medium-distance", subset$Ectomycorrhiza_exploration_type_template), "1", "0")
is.na(subset$ecm.medium.distance) <- is.na(subset$Ectomycorrhiza_exploration_type_template)
subset$ecm.mat <- ifelse(grepl("mat", subset$Ectomycorrhiza_exploration_type_template), "1", "0")
is.na(subset$ecm.mat) <- is.na(subset$Ectomycorrhiza_exploration_type_template)
subset$ecm.contact <- ifelse(grepl("contact", subset$Ectomycorrhiza_exploration_type_template), "1", "0")
is.na(subset$ecm.contact) <- is.na(subset$Ectomycorrhiza_exploration_type_template)

subset$hymenium.type.smooth <- ifelse(grepl("smooth", subset$Hymenium_type_template), "1", "0")
is.na(subset$hymenium.type.smooth) <- is.na(subset$Hymenium_type_template)
subset$hymenium.type.gills <- ifelse(grepl("gills", subset$Hymenium_type_template), "1", "0")
is.na(subset$hymenium.type.gills) <- is.na(subset$Hymenium_type_template)
subset$hymenium.type.closed <- ifelse(grepl("closed", subset$Hymenium_type_template), "1", "0")
is.na(subset$hymenium.type.closed) <- is.na(subset$Hymenium_type_template)
subset$fruitbody.tremelloid <- ifelse(grepl("tremelloid", subset$Fruitbody_type_template), "1", "0")
is.na(subset$fruitbody.tremelloid) <- is.na(subset$Fruitbody_type_template)
subset$fruitbody.perithecium <- ifelse(grepl("perithecium", subset$Fruitbody_type_template), "1", "0")
is.na(subset$fruitbody.perithecium) <- is.na(subset$Fruitbody_type_template)
subset$fruitbody.gasteroid <- ifelse(grepl("gasteroid", subset$Fruitbody_type_template), "1", "0")
is.na(subset$fruitbody.gasteroid) <- is.na(subset$Fruitbody_type_template)
subset$fruitbody.corticioid <- ifelse(grepl("corticioid", subset$Fruitbody_type_template), "1", "0")
is.na(subset$fruitbody.corticoid) <- is.na(subset$Fruitbody_type_template)
subset$fruitbody.cleistothecium <- ifelse(grepl("cleistothecium", subset$Fruitbody_type_template), "1", "0")
is.na(subset$fruitbody.cleistothecium) <- is.na(subset$Fruitbody_type_template)
subset$fruitbody.clavarioid <- ifelse(grepl("clavarioid", subset$Fruitbody_type_template), "1", "0")
is.na(subset$fruitbody.clavarioid) <- is.na(subset$Fruitbody_type_template)
subset$fruitbody.apothecium <- ifelse(grepl("apothecium", subset$Fruitbody_type_template), "1", "0")
is.na(subset$fruitbody.apothecium) <- is.na(subset$Fruitbody_type_template)
subset$fruitbody.agaricoid <- ifelse(grepl("agaricoid", subset$Fruitbody_type_template), "1", "0")
is.na(subset$fruitbody.agaricoid) <- is.na(subset$Fruitbody_type_template)
subset$growthform.yeast <- ifelse(grepl("yeast", subset$Growth_form_template), "1", "0")
is.na(subset$growthform.yeast) <- is.na(subset$Growth_form_template)
subset$growthform.thallus.photosynthetic <- ifelse(grepl("thallus_photosynthetic", subset$Growth_form_template), "1", "0")
is.na(subset$growthform.thallus.photosynthetic) <- is.na(subset$Growth_form_template)
subset$growthform.filamentous <- ifelse(grepl("filamentous_mycelium", subset$Growth_form_template), "1", "0")
is.na(subset$growthform.filamentous) <- is.na(subset$Growth_form_template)
subset$growthform.dimorphic.yeast <- ifelse(grepl("dimorphic_yeast", subset$Growth_form_template), "1", "0")
is.na(subset$growthform.dimorphic.yeast) <- is.na(subset$Growth_form_template)

##add a column of everything identified as an antagonistic organism. IE any type of pathogen and parasite.
subset$gen.pathogen <- ifelse(grepl("pathogen", subset$guild_merge), "1", "0")
subset$gen.parasite <- ifelse(grepl("parasite", subset$guild_merge), "1", "0")
subset$gen.pathogen[subset$gen.pathogen == 0] <- NA
subset$gen.parasite[subset$gen.parasite == 0] <- NA

subset <- (subset %>% mutate(general.pathogen = coalesce(gen.pathogen, gen.parasite)))

##a final data table with only the columns of interest.
subsetfinal <- subset[ -c(1:5) ]
subsetfinal <- subsetfinal %>% replace(is.na(.), 0)
subsetfinal <- mutate_all(subsetfinal, function(x) as.numeric(as.character(x)))
subset_col <- colSums(subsetfinal)
write.csv(subsetfinal, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\trait_table_raw.csv")






