##load packages
library('ggpubr')
library('tidyverse')
library('broom')

needle_d <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\needle.div_24.rds")
soil_d <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\soil.diversity_24.rds")

needle_d <- subset(needle_d, select = c(30,33,87,88))
soil_d <- subset(soil_d, select = c(36,39,97,98))

needle_d <- rownames_to_column(needle_d, "row.names")
soil_d <- rownames_to_column(soil_d, "row.names")
combo <- full_join(needle_d, soil_d)
combo <- column_to_rownames(combo, "row.names")

richness <- (ggplot(combo, aes(x = climate, y=total.richness)) + geom_point() + geom_smooth(method = "lm", size = 3, aes(color = combo$substrate)) + scale_color_manual(values = c("#006D2C", "#993404"))) + guides(color = guide_legend(title = "substrate")) + ylab("total fungal richness (# of OTUs/sample)") + theme_classic(base_size=25) + xlab("climate PCA axis 1")
diversity <- (ggplot(combo, aes(x = climate, y=total.diversity)) + geom_point() + geom_smooth(method = "lm", size = 3, linetype = "twodash", aes(color = combo$substrate)) + scale_color_manual(values = c("#006D2C","#993404"))) + ylab("total fungal diversity (H')") + guides(color = guide_legend(title = "substrate")) + theme_classic(base_size=25) + xlab("climate PCA axis 1")
richness
diversity

ab  <- ggarrange(richness, diversity, ncol=2, common.legend = TRUE, legend="right", labels = "AUTO", align = "hv", font.label = list(size = 30))
ab
##regression analyses
rich_n <- lm(total.richness ~ climate, data = needle_d)
summary(rich_n)
rich_n_table <- tidy(rich_n)

div_n <- lm(total.diversity ~ climate, data = needle_d)
summary(div_n)
div_n_table <- tidy(div_n)

rich_s <- lm(total.richness ~ climate, data = soil_d)
summary(rich_s)
rich_s_table <- tidy(rich_s)

div_s <- lm(total.diversity ~ climate, data = soil_d)
summary(div_s)
div_s_table <- tidy(div_s)

write.csv(rich_n_table, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar.richness.csv")
write.csv(div_n_table, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar.diversity.csv")

write.csv(rich_s_table, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.richness.csv")
write.csv(div_s_table, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.diversity.csv")

