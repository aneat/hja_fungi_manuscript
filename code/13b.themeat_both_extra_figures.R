##load packages
library(ggpubr)
library(tidyverse)

##load data##
resid_foliar_extra <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\resid_foliar_extra.rds")
resid_soil_extra <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\resid_soil_extra.rds")

##SOIL GRAPHS##
ed_ph <- ggplot(resid_soil_extra, aes(x = pH, y = e_resid_ph)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = substrate)) + ylab("emf diversity") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + theme(text = element_text(size=15)) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.06", size = 5, vjust = -0.5, hjust = -0.2) 
ed_ph

ed_sp <- ggplot(resid_soil_extra, aes(x = species, y = e_resid_sp)) + geom_boxplot(aes(color = substrate)) + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab("tree host species") + theme(text = element_text(size=15)) + ylab("emf diversity") + annotate("text", x=-Inf, y=-Inf, label= "p = 0.01", size = 5, vjust = -0.5, hjust = -0.2) + geom_signif(comparisons = list(c("psme", "tshe"), c("psme", "tabr")), 
                                                                                                                                                                                                                                                                                                                                                                                                                  map_signif_level=TRUE)
ed_sp

ad_sp <- ggplot(resid_soil_extra, aes(x = species, y = a_resid_sp)) + geom_boxplot(aes(color = substrate)) + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab("tree host species") + theme(text = element_text(size=15)) + ylab("amf diversity") + annotate("text", x=-Inf, y=-Inf, label= "p = 0.04", size = 5, vjust = -0.5, hjust = -0.2) + geom_signif(comparisons = list(c("psme", "tabr"), c("tabr", "tshe")), 
                                                                                                                                                                                                                                                                                                                                                                                                                 map_signif_level=TRUE)
ad_sp

acwm_sp <- ggplot(resid_soil_extra, aes(x = species, y = a_resid_cwm)) + geom_boxplot(aes(color = substrate)) + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab("tree host species") + theme(text = element_text(size=15)) + ylab("amf cwm") + annotate("text", x=-Inf, y=-Inf, label= "p = 0.04", size = 5, vjust = -0.5, hjust = -0.2) +   geom_signif(comparisons = list(c("tshe", "psme")), 
                                                                                                                                                                                                                                                                                                                                                                                                              map_signif_level=TRUE)
acwm_sp

sapcwm_ph <- ggplot(resid_soil_extra, aes(x = pH, y = sap_resid_cwm)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = substrate)) + ylab("saprobe cwm") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + theme(text = element_text(size=15)) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.01", size = 5, vjust = -0.5, hjust = -0.2) 
sapcwm_ph

ehs_ph <- ggplot(resid_soil_extra, aes(x = pH, y = ecm_resid_hs)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = substrate)) + ylab("emf host specificifty") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + theme(text = element_text(size=15)) + annotate("text", x=-Inf, y=-Inf, label= "p < 0.00", size = 5, vjust = -0.5, hjust = -0.2) 
ehs_ph

##FOLIAR GRAPHS##
fpd_ph <- ggplot(resid_foliar_extra, aes(x = pH, y = fpath_div)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = substrate)) + ylab("foliar pathogen diversity") + theme_classic(base_size=15) + scale_color_manual(values=c("#006D2C")) + guides(color = "none") + theme(text = element_text(size=15)) + annotate("text", x=-Inf, y=-Inf, label= "p < 0.00", size = 5, vjust = -0.5, hjust = -0.2) 
fpd_ph

fsd_sp <- ggplot(resid_foliar_extra, aes(x = species, y = fsap_div)) + geom_boxplot(aes(color = substrate)) + theme_classic(base_size=15) + scale_color_manual(values=c("#006D2C")) + guides(color = "none") + xlab("tree host species") + theme(text = element_text(size=15)) + ylab("foliar saprobe diversity") + annotate("text", x=-Inf, y=-Inf, label= "p < 0.00", size = 5, vjust = -0.5, hjust = -0.2) +  geom_signif(comparisons = list(c("psme", "tshe"), c("tabr", "tshe")), 
                                                                                                                                                                                                                                                                                                                                                                                                                              map_signif_level=TRUE)
fsd_sp

fpcwm_ph <- ggplot(resid_foliar_extra, aes(x = pH, y = fpath_cwm)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = substrate)) + ylab("foliar pathogen cwm") + theme_classic(base_size=15) + scale_color_manual(values=c("#006D2C")) + guides(color = "none") + theme(text = element_text(size=15)) + annotate("text", x=-Inf, y=-Inf, label= "p < 0.00", size = 5, vjust = -0.5, hjust = -0.2) 
fpcwm_ph

fscwm_sp <- ggplot(resid_foliar_extra, aes(x = species, y = fsap_cwm)) + geom_boxplot(aes(color = substrate)) + theme_classic(base_size=15) + scale_color_manual(values=c("#006D2C")) + guides(color = "none") + xlab("tree host species") + theme(text = element_text(size=15)) + ylab("foliar saprobe cwm") + annotate("text", x=-Inf, y=-Inf, label= "p < 0.00", size = 5, vjust = -0.5, hjust = -0.2) +   geom_signif(comparisons = list(c("psme", "tabr"), c("tabr", "tshe")), 
                                                                                                                                                                                                                                                                                                                                                                                                                          map_signif_level=TRUE)
fscwm_sp

fscwm_nut <- ggplot(resid_foliar_extra, aes(x = foliar, y = fsap_cwm_nut)) + geom_point(color = "darkgrey") + geom_smooth(method = "lm", size = 2, aes(color = substrate)) + xlab("foliar nutrients") + ylab("foliar saprobe cwm") + theme_classic(base_size=15) + scale_color_manual(values=c("#006D2C")) + guides(color = "none") + theme(text = element_text(size=15)) + annotate("text", x=-Inf, y=-Inf, label= "p = 0.03", size = 5, vjust = -0.5, hjust = -0.2) 
fscwm_nut

fpspe_sp <- ggplot(resid_foliar_extra, aes(x = species, y = fpath_spec)) + geom_boxplot(aes(color = substrate)) + theme_classic(base_size=15) + scale_color_manual(values=c("#006D2C")) + guides(color = "none") + xlab("tree host species") + theme(text = element_text(size=15)) + ylab("foliar pathogen specificity") + annotate("text", x=-Inf, y=-Inf, label= "p < 0.00", size = 5, vjust = -0.5, hjust = -0.2) +   geom_signif(comparisons = list(c("psme", "tabr"), c("tabr", "tshe")), 
                                                                                                                                                                                                                                                                                                                                                                                                                                     map_signif_level=TRUE)
fpspe_sp

fsspe_sp <- ggplot(resid_foliar_extra, aes(x = species, y = fsap_spec)) + geom_boxplot(aes(color = substrate)) + theme_classic(base_size=15) + scale_color_manual(values=c("#006D2C")) + guides(color = "none") + xlab("tree host species") + theme(text = element_text(size=15)) + ylab("foliar saprobe specificity") + annotate("text", x=-Inf, y=-Inf, label= "p < 0.00", size = 5, vjust = -0.5, hjust = -0.2) +   geom_signif(comparisons = list(c("psme", "tabr"), c("tabr", "tshe")), 
                                                                                                                                                                                                                                                                                                                                                                                                                                   map_signif_level=TRUE)
fsspe_sp

soil_plots <- ggarrange(ed_ph, ed_sp, ehs_ph, ad_sp, acwm_sp, sapcwm_ph,
                   ncol = 3, nrow = 2, align = c("hv"), labels = c("A", "B", "C", "D", "E", "F"))
soil_plots
foliar_plots <- ggarrange(fpd_ph, fpspe_sp, fsd_sp, fscwm_sp, fscwm_nut, fsspe_sp,
                        ncol = 3, nrow = 2, align = c("hv"), labels = c("A", "B", "C", "D", "E", "F"))
foliar_plots
