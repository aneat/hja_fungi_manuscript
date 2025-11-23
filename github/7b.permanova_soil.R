##read in packages
library(tidyverse)
library(factoextra)
library(AICcPermanova)
library(vegan)
library(ecole)
library(Hmisc)
library(corrplot)
library(ggpubr)

###FOR SOILS###
##read in matrices
new_otu <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\new_otu_soil_24.rds")
env1 <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\env1_soil_24.rds")
tra <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\tra_soil_24.rds")
dist <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\dist_soil_24.rds")
#per.table <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\per.table.soil_24.rds")

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
soil1 <- (fviz_pca_var(data.pca, col.var = "#993404", labelsize = 5, title = "", arrowsize = 1.5, repel = TRUE, ggtheme = theme_classic(base_size = 20)))
soil_final <- (soil1 + coord_fixed(ratio = 1) + theme_light(base_size = 20))
soil2 <- (fviz_cos2(data.pca, choice = "var", axes = 1, ggtheme = theme_classic(base_size=20), labelsize = 8, fill = "#993404", color = "black", title = FALSE))
soil_final <- (soil1 + coord_fixed(ratio = 1) + theme_light(base_size = 20))

##foliar nutrients pca
foliar.pca <- princomp(foliar.nutrients)
summary(foliar.pca)
foliar.pca.scores <- foliar.pca$scores
foliar.pca.scores <- as.data.frame(foliar.pca.scores)
foliar.pca$loadings[, 1:2]
ned1 <- (fviz_pca_var(foliar.pca, col.var = "#006D2C", ggtheme = theme_classic(base_size = 20), arrowsize = 1.5, labelsize = 5, title = "", repel = TRUE))
ned2 <- (fviz_cos2(foliar.pca, choice = "var", axes = 1, ggtheme = theme_classic(base_size = 20), labelsize = 8, fill = "#006D2C", color = "black", title = FALSE))
ned_final <- (ned1 + xlim(-200, 50) + ylim(-200, 50) + theme_light(base_size = 20))

##ned2
##ggarrange(ned1, ned2, labels = "AUTO", font.label = list(size = 20))

##climate pca
climate.pca <- princomp(climate)
summary(climate.pca)
climate.pca$loadings <- climate.pca$loadings*-1
climate.pca$scores <- climate.pca$scores * -1
climate.pca.scores <- climate.pca$scores
climate.pca.scores <- as.data.frame(climate.pca.scores)
climate.pca$loadings[, 1:2]
c1 <- (fviz_pca_var(climate.pca, col.var = "royalblue", ggtheme = theme_classic(base_size = 20), labelsize = 5, title = "", arrowsize = 1.5, repel = TRUE))
c2 <- (fviz_cos2(climate.pca, choice = "var", axes = 1, ggtheme = theme_classic(base_size = 20), title = FALSE, labelsize = 8))
climate_final <- (c1 + coord_fixed(ratio = 1) + theme_light(base_size = 20) + xlim(-2, 1) + ylim(-2,1))
climate_final
climate.pca$scores

##create permanova table
soil.pca.scores <- select(soil.pca.scores, c(1))
climate.pca.scores <- select(climate.pca.scores, c(1))
foliar.pca.scores <- select(foliar.pca.scores, c(1))
per.table <- select(env1, c(5,6,10,34,4))
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
saveRDS(per.table, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\per.table.soil_24.rds")

corr_data <- subset(per.table, select = c(3,4,6:8))
p_values <- rcorr(as.matrix(corr_data))
print(p_values)
c <- cor(corr_data)
cor.plot <- corrplot(c, method = "number", type = "lower", tl.col = "black")
cor.plot
cl.ele <- ggplot(per.table, aes(x = climate, y = elevation)) + stat_smooth(method = "lm", size = 3, color = "royalblue") + theme_classic(base_size = 20) + xlab("climate PCA axis 1") + ylab("elevation (meters)") + annotate("text", x=-2, y=500, label= "r = 0.96", size = 5)
cl.ele


##MAKE PCAS SUPPLEMENTAL FIGURE##
pcas <- ggarrange(climate_final, soil_final, ned_final, cl.ele, ncol = 2, nrow = 2, align = c("hv"))
pcas
ggsave("pcas.png", pcas, units = "px", width = 4800, height = 5000, dpi = 550)


##permanvoa
a <- adonis2(dist ~ species + climate + soil + pH, data=per.table, by='margin')
b <- adonis2(dist ~ species + climate + pH, data=per.table, by='margin')
c <- adonis2(dist ~ species + climate, data=per.table, by='margin')
d <- adonis2(dist ~ climate, data = per.table, by= "margin")
saveRDS(b, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\permanova_soil_24.rds")

##determine bestfit model. The lower the AIC value, the better the model fit.

AICc_permanova2(a)
AICc_permanova2(b)
AICc_permanova2(c)
AICc_permanova2(d)
soil_perm_tab <- d
write.csv(soil_perm_tab, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.perm_24.csv")
d
##THE END##

