##read in packages
library(tidyverse)
library(factoextra)
library(AICcPermanova)
library(vegan)
library(ecole)

###FOR NEEDLES###

##read in matrices
new_otu <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\new_otu_needle_24.rds")
env1 <-  readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\env1_needle_24.rds")
tra <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\tra_needle_24.rds")
dist <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\dist_needle_24.rds")

##subset data frames for soil, foliar, and climate PCA. DO NOT INCLUDE SULFUR. IT LOOKS LIKE A DATA MISSENTRY.
soil.characteristics <- select(env1, c(27:44))
rownames(soil.characteristics) <- env1$row.names
foliar.nutrients <- select(env1, c(13:26))
rownames(foliar.nutrients) <- env1$row.names
climate <- select(env1, c(47:52))
rownames(climate) <- env1$row.names

##soil characteristics pca
data.pca <- princomp(soil.characteristics)
summary(data.pca)
soil.pca.scores <- data.pca$scores
soil.pca.scores <- as.data.frame(soil.pca.scores)
data.pca$loadings[, 1:2]
soil1 <- (fviz_pca_var(data.pca, col.var = "black", labelsize = 8, title = "", arrowsize = 1.5, repel = TRUE, ggtheme = theme_classic(base_size = 30)))
soil2 <- (fviz_cos2(data.pca, choice = "var", axes = 1, ggtheme = theme_classic(base_size=30), labelsize = 8, fill = "#993404", color = "black", title = FALSE))
soil2
##ggarrange(soil1, soil2, labels = "AUTO", font.label = list(size = 30))

##foliar nutrients pca
foliar.pca <- princomp(foliar.nutrients)
summary(foliar.pca)
foliar.pca.scores <- foliar.pca$scores
foliar.pca.scores <- as.data.frame(foliar.pca.scores)
foliar.pca$loadings[, 1:2]
ned1 <- (fviz_pca_var(foliar.pca, col.var = "black", ggtheme = theme_classic(base_size = 30), arrowsize = 1.5, labelsize = 8, title = "", repel = TRUE))
ned2 <- (fviz_cos2(foliar.pca, choice = "var", axes = 1, ggtheme = theme_classic(base_size = 30), labelsize = 8, fill = "#006D2C", color = "black", title = FALSE))
##ned1
##ned2
##ggarrange(ned1, ned2, labels = "AUTO", font.label = list(size = 30))

##climate pca
climate.pca <- princomp(climate)
summary(climate.pca)
climate.pca$loadings <- climate.pca$loadings*-1
climate.pca$scores <- climate.pca$scores * -1
climate.pca.scores <- climate.pca$scores
climate.pca.scores <- as.data.frame(climate.pca.scores)
climate.pca$loadings[, 1:2]
c1 <- (fviz_pca_var(climate.pca, col.var = "black", ggtheme = theme_classic(base_size = 30), labelsize = 8, title = "", arrowsize = 1.5, repel = TRUE))
c2 <- (fviz_cos2(climate.pca, choice = "var", axes = 1, ggtheme = theme_classic(base_size = 30), title = FALSE, labelsize = 8))
c1
#ggarrange(c1, c2, labels = "AUTO", font.label = list(size = 30))

##create permanova table
soil.pca.scores <- select(soil.pca.scores, c(1))
climate.pca.scores <- select(climate.pca.scores, c(1))
foliar.pca.scores <- select(foliar.pca.scores, c(1))
per.table <- select(env1, c(5,6,34,10,4))
per.table <- merge(per.table, soil.pca.scores,
                   by = 'row.names', all = TRUE)
per.table <- column_to_rownames(per.table, "Row.names")
per.table <- merge(per.table, foliar.pca.scores,
                          by = 'row.names', all = TRUE)
per.table <- column_to_rownames(per.table, "Row.names")
per.table <- merge(per.table, climate.pca.scores,
                   by = 'row.names', all = TRUE)
per.table <- column_to_rownames(per.table, "Row.names")
per.table$soil <- per.table$Comp.1.x
per.table$foliar <- per.table$Comp.1.y
per.table$climate <- per.table$Comp.1

per.table <- select(per.table, -c(6:8))
per.table <- per.table %>% filter(row.names(per.table) != "Needles-G5")
saveRDS(per.table, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\per.table.needle_24.rds")

##permanvoa
a <- adonis2(dist ~ species + climate + pH + foliar, data=per.table, by='margin')
b <- adonis2(dist ~ species + climate + pH, data=per.table, by='margin')
c <- adonis2(dist ~ species + climate, data=per.table, by='margin')
d <- adonis2(dist ~ species, data = per.table, by= "margin")
saveRDS(b, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\permanova_needle_24.rds")
a
b
c
d
##determine best fit model. The lower the AIC value, the better the model fit.
AICc_permanova2(a)
AICc_permanova2(b)
AICc_permanova2(c)
AICc_permanova2(d)
a
c
needle_perm_tab <- c
write.csv(needle_perm_tab, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\needle.perm_24.csv")

##THE END##

