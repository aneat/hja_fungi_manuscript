##load packages
library(ggpubr)
library(tidyverse)


##load data
resid_soil <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\resid_soil_24.rds")
resid_foliar <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\resid_foliar_24.rds")

##SOIL GRAPHS##
##SOIL DIVERSITY##
spd <- ggplot(resid_soil, aes(x = climate, y = p_resid_d)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = substrate)) + ylab("shannon's diversity (H')") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.03", size = 5, vjust = -0.5, hjust = -0.2) + ggtitle("soil pathogens") + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5))
spd

sed <- ggplot(resid_soil, aes(x = climate, y = e_resid_d)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = substrate), linetype = "twodash") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.92", size = 5, vjust = -0.5, hjust = -0.2) + ggtitle("emf") + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + ylab("") + scale_y_continuous(expand = expansion(add = 0.4))
sed

sad <- ggplot(resid_soil, aes(x = climate, y = a_resid_d)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = substrate), linetype = "twodash") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ") + ylab("") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.30", size = 5, vjust = -0.5, hjust = -0.2) + ggtitle("amf") + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + scale_y_continuous(expand = expansion(add = 1))
sad

ssd <- ggplot(resid_soil, aes(x = climate, y = s_resid_d)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = substrate), linetype = "twodash") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.20", size = 5, vjust = -0.5, hjust = -0.2) + ggtitle("saprobe") + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + ylab("")
ssd

##SOIL CWM##
spcwm <- ggplot(resid_soil, aes(x = climate, y = p_resid_cwm)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = substrate)) + ylab("community weighted mean") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.02", size = 5, vjust = -0.5, hjust = -0.2) + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + scale_y_continuous(expand = expansion(add = 0.02))
spcwm

secwm <- ggplot(resid_soil, aes(x = climate, y = e_resid_cwm)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = substrate), linetype = "twodash") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.61", size = 5, vjust = -0.5, hjust = -0.2) + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + ylab("") + scale_y_continuous(expand = expansion(add = 0.1))
secwm

sacwm <- ggplot(resid_soil, aes(x = climate, y = a_resid_cwm)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = substrate), linetype = "twodash") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ") + ylab("") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.13", size = 5, vjust = -0.5, hjust = -0.2) + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + scale_y_continuous(expand = expansion(add = 0.0005))
sacwm

sscwm <- ggplot(resid_soil, aes(x = climate, y = s_resid_cwm)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = substrate), linetype = "twodash") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.45", size = 5, vjust = -0.5, hjust = -0.2) + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + ylab("") + scale_y_continuous(expand = expansion(add = 0.1))
sscwm


##SOIL SPECIFICITY##
sps <- ggplot(resid_soil, aes(x = climate, y = p_resid_spe)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = substrate)) + ylab("host specificity index") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.02", size = 5, vjust = -0.5, hjust = -0.2) + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + scale_y_continuous(expand = expansion(add = 0.02))
sps

ses <- ggplot(resid_soil, aes(x = climate, y = e_resid_spe)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = substrate)) + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ") + ylab("") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p < 0.00", size = 5, vjust = -0.5, hjust = -0.2) + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5))
ses

sas <- ggplot(resid_soil, aes(x = climate, y = a_resid_spe)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = substrate), linetype = "twodash") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ") + ylab("") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.80", size = 5, vjust = -0.5, hjust = -0.2) + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5))
sas

sss <- ggplot(resid_soil, aes(x = climate, y = s_resid_spe)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = substrate)) + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p < 0.00", size = 5, vjust = -0.5, hjust = -0.2) + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + ylab("")
sss

#FOLIAR GRAPHS##
##foliar DIVERSITY##
fpd <- ggplot(resid_foliar, aes(x = climate, y = p_resid_d)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = substrate), linetype = "twodash") + theme_classic(base_size=15) + scale_color_manual(values=c("#006D2C")) + guides(color = "none") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.90", size = 5, vjust = -0.5, hjust = -0.2) + ggtitle("foliar pathogens") + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + ylab("") + scale_y_continuous(expand = expansion(add = 0.2))
fpd

fsd <- ggplot(resid_foliar, aes(x = climate, y = s_resid_d)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = substrate), linetype = "twodash") + theme_classic(base_size=15) + scale_color_manual(values=c("#006D2C")) + guides(color = "none") + xlab(" ") + theme(text = element_text(size=15), element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.84", size = 5, vjust = -0.5, hjust = -0.2) + ggtitle("saprobe") + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + ylab("")
fsd

##foliar CWM##
fpcwm <- ggplot(resid_foliar, aes(x = climate, y = p_resid_cwm)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = substrate), linetype = "twodash") + theme_classic(base_size=15) + scale_color_manual(values=c("#006D2C")) + guides(color = "none") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.35", size = 5, vjust = -0.5, hjust = -0.2) + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + ylab("") + scale_y_continuous(expand = expansion(add = 0.1))
fpcwm

fscwm <- ggplot(resid_foliar, aes(x = climate, y = s_resid_cwm)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = substrate), linetype = "twodash") + theme_classic(base_size=15) + scale_color_manual(values=c("#006D2C")) + guides(color = "none") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.56", size = 5, vjust = -0.5, hjust = -0.2) + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + ylab("") + scale_y_continuous(expand = expansion(add = 0.1))
fscwm

##foliar SPECIFICITY##
fps <- ggplot(resid_foliar, aes(x = climate, y = p_resid_spe)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = substrate)) + ylab("") + theme_classic(base_size=15) + scale_color_manual(values=c("#006D2C")) + guides(color = "none") + xlab(" ") + theme(text = element_text(size=15), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p < 0.00", size = 5, vjust = -0.5, hjust = -0.2) + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + scale_y_continuous(expand = expansion(add = 0.05))
fps

fss <- ggplot(resid_foliar, aes(x = climate, y = s_resid_spe)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = substrate)) + theme_classic(base_size=15) + scale_color_manual(values=c("#006D2C")) + guides(color = "none") + xlab(" ") + theme(text = element_text(size=15), panel.background = element_rect(fill = "white"),
                                                                                                                                                                                                                                                                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y = element_text(size= 15, face = "bold")) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.21", size = 5, vjust = -0.5, hjust = -0.2) + theme(plot.title = element_text(size= 15, face="bold", hjust = 0.5)) + ylab("") + scale_y_continuous(expand = expansion(add = 0.02))
fss

plots <- ggarrange(spd, sed, sad, fpd, spcwm, secwm, sacwm, fpcwm, sps, ses, sas, fps,
                   ncol = 4, nrow = 3, align = c("hv"))
plots

