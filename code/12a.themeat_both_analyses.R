##load packages##
library(ggpubr)
library(gridExtra)
library(broom)
library(tidyverse)
library(AICcmodavg)
library(MuMIn)

##load soil data
soil.diversity <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\soil.diversity_24.rds")
soil.spec <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\soil.spec_24.rds")

##load foliar data
foliar.diversity <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\needle.div_24.rds")
foliar.spec <- readRDS("\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\foliar.spec_24.rds")

soil.diversity <- subset(soil.diversity, select = c(1,2,3,8,33,35:39,89:98))

merge_soil <- merge(soil.diversity, soil.spec, by = 0)
merge_soil <- column_to_rownames(merge_soil, "Row.names")
resid_soil <- subset(merge_soil, select = c(7,10))

foliar.diversity <- subset(foliar.diversity, select = c(2,3,27:33,83:88))

merge_foliar <- merge(foliar.diversity, foliar.spec, by = 0)
merge_foliar <- column_to_rownames(merge_foliar, "Row.names")
resid_foliar <- subset(merge_foliar, select = c(6,9))

##GENERAL DIVERSITY & RICHNESS##
sr <- lm(total.richness ~ climate, data = merge_soil)
sd <- lm(total.diversity ~ climate, data = merge_soil)
summary(sr)
summary(sd)
sr <- tidy(sr)
sd <- tidy(sd)


fr <- lm(total.richness ~ climate, data = merge_foliar)
fd <- lm(total.diversity ~ climate, data = merge_foliar)
summary(fr)
summary(fd)

fr <- tidy(fr)
fd <- tidy(fd)

write.csv(sr, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.richness.csv")
write.csv(sd, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.diversity.csv")

write.csv(fr, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar.richness.csv")
write.csv(fd, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar.diversity.csv")

##SOIL DIVERSITY##
##PATHOGENS##
p <- lm(path.diversity ~ species + pH + soil + climate, data = merge_soil)
p1 <- lm(path.diversity ~ species + climate + soil, data = merge_soil)
p2 <- lm(path.diversity ~ climate + soil, data = merge_soil)
p3 <- lm(path.diversity ~ climate, data = merge_soil)
models_p <- list(p, p1, p2, p3)
aic_p <- aictab(cand.set = models_p)
aic_p
summary(p3)
model_full <- tidy(p3)
model_full
write.csv(model_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.path.diversity.csv")
write.csv(aic_p, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.path.diversity.aic.csv")

p_full <- lm(path.diversity ~ species + pH + soil, data = merge_soil)
summary(p_full)
p_resid_d <- resid(p_full)
p_climate <- lm(p_resid_d ~ climate, data = merge_soil)
summary(p_climate)

p_resid_d <- as.data.frame(p_resid_d)
resid_soil <- cbind(resid_soil, p_resid_d)

pd <- ggplot(resid_soil, aes(x = climate, y = p_resid_d)) + geom_point() + geom_smooth(method = "lm", size = 3) + ylab("") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ")
pd

##ECM##
e <- lm(ecm.diversity ~ species + pH + soil + climate, data = merge_soil)
e1 <- lm(ecm.diversity ~ species + pH + soil, data = merge_soil)
e2 <- lm(ecm.diversity ~ species + pH, data = merge_soil)
e3 <- lm(ecm.diversity ~ pH, data = merge_soil)
models_e <- list(e, e1, e2, e3)
aic_e <- aictab(cand.set = models_e)
aic_e
summary(e2)
model_full <- tidy(e2)
write.csv(model_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.ecm.diversity.csv")
write.csv(aic_e, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.ecm.diversity.aic.csv")

e_full <- lm(ecm.diversity ~ species + pH + soil, data = merge_soil)
summary(e_full)
e_resid_d <- resid(e_full)
e_climate <- lm(e_resid_d ~ climate, data = merge_soil)
summary(e_climate)
e_resid_d <- as.data.frame(e_resid_d)
resid_soil <- cbind(resid_soil, e_resid_d)

ed <- ggplot(resid_soil, aes(x = climate, y = e_resid_d)) + geom_point() + geom_smooth(method = "lm", size = 3) + ylab("") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ")
ed

##AM##
a <- lm(am.diversity ~ species + pH + soil + climate, data = merge_soil)
a1 <- lm(am.diversity ~ species + pH + climate, data = merge_soil)
a2 <- lm(am.diversity ~ species + pH, data = merge_soil)
a3 <- lm(am.diversity ~ species, data = merge_soil)
models_a <- list(a, a1, a2, a3)
aic_a <- aictab(cand.set = models_a)
aic_a
summary(a)
model_avg <- summary(model.avg(models_a))
model_full <- model_avg$coefmat.full
model_full
write.csv(model_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.am.diversity.csv")
write.csv(aic_a, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.am.diversity.aic.csv")

a_full <- lm(am.diversity ~ species + pH + soil, data = merge_soil)
summary(a_full)
a_resid_d <- resid(a_full)
a_climate <- lm(a_resid_d ~ climate, data = merge_soil)
summary(a_climate)
a_resid_d <- as.data.frame(a_resid_d)
resid_soil <- cbind(resid_soil, a_resid_d)

ad <- ggplot(resid_soil, aes(x = climate, y = a_resid_d)) + geom_point() + geom_smooth(method = "lm", size = 3) + ylab("") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ")
ad


##SAPROBES##
s <- lm(sap.diversity ~ species + pH + soil + climate, data = merge_soil)
s1 <- lm(sap.diversity ~ species + pH + climate, data = merge_soil)
s2 <- lm(sap.diversity ~ species + pH, data = merge_soil)
s3 <- lm(sap.diversity ~ pH, data = merge_soil)
models_s <- list(s, s1, s2, s3)
aic_s <- aictab(cand.set = models_s)
aic_s
summary(s3)
model_full <- tidy(s3)
write.csv(model_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.sap.diversity.csv")
write.csv(aic_s, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.sap.diversity.aic.csv")

s_full <- lm(sap.diversity ~ species + pH + soil, data = merge_soil)
summary(s_full)
s_resid_d <- resid(s_full)
s_climate <- lm(s_resid_d ~ climate, data = merge_soil)
summary(s_climate)
s_resid_d <- as.data.frame(s_resid_d)
resid_soil <- cbind(resid_soil, s_resid_d)

sd <- ggplot(resid_soil, aes(x = climate, y = s_resid_d)) + geom_point() + geom_smooth(method = "lm", size = 3) + ylab("") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ")
sd
##THE END OF PT1: DIVERSITY##



##SOIL CWM##
##PATHOGENS##
p <- lm(pathogen ~ species + pH + soil + climate, data = merge_soil)
p1 <- lm(pathogen ~ species + climate + soil, data = merge_soil)
p2 <- lm(pathogen ~ climate + soil, data = merge_soil)
p3 <- lm(pathogen ~ climate, data = merge_soil)
models_p <- list(p, p1, p2, p3)
aic_p <- aictab(cand.set = models_p)
aic_p
summary(p)
models_p <- list(p2,p3)
model_avg <- summary(model.avg(models_p))
m_coefs <- coefTable(model_avg, full=TRUE, adjust.se = TRUE, df = TRUE)
coefnames <- row.names(m_coefs)
m_coefs <- as.data.frame(m_coefs)
m_coefs <- mutate(m_coefs,
                  coefficient = coefnames,
                  t = Estimate / `Std. Error`,
                  p = pt(abs(t), df = 65, lower.tail = FALSE),
                  lower95 = Estimate - 1.96 * `Std. Error`,
                  upper95 = Estimate + 1.96 * `Std. Error`) %>% select(coefficient, Estimate, `Std. Error`, t, p, lower95, upper95)
write.csv(m_coefs, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.path.cwm.csv")
write.csv(aic_p, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.path.cwm.aic.csv")

p_full <- lm(pathogen ~ species + pH + soil, data = merge_soil)
summary(p_full)
p_resid_cwm <- resid(p_full)
p_climate <- lm(p_resid_cwm ~ climate, data = merge_soil)
summary(p_climate)
p_resid_cwm <- as.data.frame(p_resid_cwm)
resid_soil <- cbind(resid_soil, p_resid_cwm)

p_cwm <- ggplot(resid_soil, aes(x = climate, y = p_resid_cwm)) + geom_point() + geom_smooth(method = "lm", size = 3) + ylab("") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ")
p_cwm

##ECM##
e <- lm(ecm ~ species + pH + soil + climate, data = merge_soil)
e1 <- lm(ecm ~ species + pH + soil, data = merge_soil)
e2 <- lm(ecm ~ species + pH, data = merge_soil)
e3 <- lm(ecm ~ soil, data = merge_soil)
models_e <- list(e, e1, e2, e3)
aic_e <- aictab(cand.set = models_e)
aic_e
summary(e3)
model_full <- tidy(e3)
write.csv(model_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.ecm.cwm.csv")
write.csv(aic_e, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.ecm.cwm.aic.csv")

e_full <- lm(ecm ~ species + pH + soil, data = merge_soil)
summary(e_full)
e_resid_cwm <- resid(e_full)
e_climate <- lm(e_resid_cwm ~ climate, data = merge_soil)
summary(e_climate)
e_resid_cwm <- as.data.frame(e_resid_cwm)
resid_soil <- cbind(resid_soil, e_resid_cwm)

e_cwm <- ggplot(resid_soil, aes(x = climate, y = e_resid_cwm)) + geom_point() + geom_smooth(method = "lm", size = 3) + ylab("") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ")
e_cwm

##AM##
a <- lm(am ~ species + pH + soil + climate, data = merge_soil)
a1 <- lm(am ~ species + pH + climate, data = merge_soil)
a2 <- lm(am ~ species + pH, data = merge_soil)
a3 <- lm(am ~ species, data = merge_soil)
models_a <- list(a, a1, a2, a3)
aic_a <- aictab(cand.set = models_a)
aic_a
summary(a)
models_a <- list(a1, a2, a3)
model_avg <- summary(model.avg(models_a))
m_coefs <- coefTable(model_avg, full=TRUE, adjust.se = TRUE, df = TRUE)
coefnames <- row.names(m_coefs)
m_coefs <- as.data.frame(m_coefs)
m_coefs <- mutate(m_coefs,
                  coefficient = coefnames,
                  t = Estimate / `Std. Error`,
                  p = pt(abs(t), df = 64, lower.tail = FALSE),
                  lower95 = Estimate - 1.96 * `Std. Error`,
                  upper95 = Estimate + 1.96 * `Std. Error`) %>% select(coefficient, Estimate, `Std. Error`, t, p, lower95, upper95)

write.csv(m_coefs, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.am.cwm.csv")
write.csv(aic_a, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.am.cwm.aic.csv")

a_full <- lm(am ~ species + pH + soil, data = merge_soil)
summary(a_full)
a_resid_cwm <- resid(a_full)
a_climate <- lm(a_resid_cwm ~ climate, data = merge_soil)
summary(a_climate)
a_resid <- as.data.frame(a_resid_cwm)
resid_soil <- cbind(resid_soil, a_resid_cwm)

a_cwm <- ggplot(resid_soil, aes(x = climate, y = a_resid_cwm)) + geom_point() + geom_smooth(method = "lm", size = 3) + ylab("") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ")
a_cwm


##SAPROBES##
s <- lm(saprobe ~ species + pH + soil + climate, data = merge_soil)
s1 <- lm(saprobe ~ species + pH + climate, data = merge_soil)
s2 <- lm(saprobe ~ species + pH, data = merge_soil)
s3 <- lm(saprobe ~ pH, data = merge_soil)
models_s <- list(s, s1, s2, s3)
aic_s <- aictab(cand.set = models_s)
aic_s
summary(s3)
model_full <- tidy(s3)
write.csv(model_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.sap.cwm.csv")
write.csv(aic_s, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.sap.cwm.aic.csv")

s_full <- lm(saprobe ~ species + pH + soil, data = merge_soil)
summary(s_full)
s_resid_cwm <- resid(s_full)
s_climate <- lm(s_resid_cwm ~ climate, data = merge_soil)
summary(s_climate)
s_resid_cwm <- as.data.frame(s_resid_cwm)
resid_soil <- cbind(resid_soil, s_resid_cwm)

s_cwm <- ggplot(resid_soil, aes(x = climate, y = s_resid_cwm)) + geom_point() + geom_smooth(method = "lm", size = 3) + ylab("") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ")
s_cwm

##THE END OF PT 2: CWM##



##SOIL SPECIFICITY##
##PATHOGENS##
p <- lm(path.host.spe ~ species + pH + soil + climate, data = merge_soil)
p1 <- lm(path.host.spe ~ species + climate + soil, data = merge_soil)
p2 <- lm(path.host.spe ~ climate + species, data = merge_soil)
p3 <- lm(path.host.spe ~ climate, data = merge_soil)
models_p <- list(p, p1, p2, p3)
aic_p <- aictab(cand.set = models_p)
aic_p
models_p <- list(p2, p3)
model_avg <- summary(model.avg(models_p))
model_full <- model_avg$coefmat.full
m_coefs <- coefTable(model_avg, full=TRUE, adjust.se = TRUE, df = TRUE)
coefnames <- row.names(m_coefs)
m_coefs <- as.data.frame(m_coefs)
m_coefs <- mutate(m_coefs,
                   coefficient = coefnames,
                   t = Estimate / `Std. Error`,
                   p = pt(abs(t), df = 65, lower.tail = FALSE),
                   lower95 = Estimate - 1.96 * `Std. Error`,
                   upper95 = Estimate + 1.96 * `Std. Error`) %>% select(coefficient, Estimate, `Std. Error`, t, p, lower95, upper95)
write.csv(m_coefs, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.path.spe.csv")
write.csv(aic_p, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.path.spe.aic.csv")

p_full <- lm(path.host.spe ~ species + pH + soil, data = merge_soil)
summary(p_full)
p_resid_spe <- resid(p_full)
p_climate <- lm(p_resid_spe ~ climate, data = merge_soil)
summary(p_climate)
p_resid_spe <- as.data.frame(p_resid_spe)
resid_soil <- cbind(resid_soil, p_resid_spe)

p_spe <- ggplot(resid_soil, aes(x = climate, y = p_resid_spe)) + geom_point() + geom_smooth(method = "lm", size = 3) + ylab("") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ")
p_spe

##ECM##
e <- lm(ecm.host.spe ~ species + pH + soil + climate, data = merge_soil)
e1 <- lm(ecm.host.spe ~ soil + pH + climate, data = merge_soil)
e2 <- lm(ecm.host.spe ~ climate + soil, data = merge_soil)
e3 <- lm(ecm.host.spe ~ soil, data = merge_soil)
models_e <- list(e, e1, e2, e3)
aic_e <- aictab(cand.set = models_e)
aic_e
models_e <- list(e1, e3)
model_avg <- summary(model.avg(models_e))
m_coefs <- coefTable(model_avg, full=TRUE, adjust.se = TRUE, df = TRUE)
coefnames <- row.names(m_coefs)
m_coefs <- as.data.frame(m_coefs)
m_coefs <- mutate(m_coefs,
                  coefficient = coefnames,
                  t = Estimate / `Std. Error`,
                  p = pt(abs(t), df = 64, lower.tail = FALSE),
                  lower95 = Estimate - 1.96 * `Std. Error`,
                  upper95 = Estimate + 1.96 * `Std. Error`) %>% select(coefficient, Estimate, `Std. Error`, t, p, lower95, upper95)
write.csv(m_coefs, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.ecm.spe.csv")
write.csv(aic_e, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.ecm.spe.aic.csv")

e_full <- lm(ecm.host.spe ~ species + pH + soil, data = merge_soil)
summary(e_full)
e_resid_spe <- resid(e_full)
e_climate <- lm(e_resid_spe ~ climate, data = merge_soil)
summary(e_climate)
e_resid_spe <- as.data.frame(e_resid_spe)
resid_soil <- cbind(resid_soil, e_resid_spe)

e_spe <- ggplot(resid_soil, aes(x = climate, y = e_resid_spe)) + geom_point() + geom_smooth(method = "lm", size = 3) + ylab("") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ")
e_spe

##AM##
a <- lm(am.host.spe ~ species + pH + soil + climate, data = merge_soil)
a1 <- lm(am.host.spe ~ species + pH + climate, data = merge_soil)
a2 <- lm(am.host.spe ~ species + pH, data = merge_soil)
a3 <- lm(am.host.spe ~ species, data = merge_soil)
models_a <- list(a, a1, a2, a3)
aic_a <- aictab(cand.set = models_a)
aic_a
summary(a)
model_avg <- summary(model.avg(models_a))
m_coefs <- coefTable(model_avg, full=TRUE, adjust.se = TRUE, df = TRUE)
coefnames <- row.names(m_coefs)
m_coefs <- as.data.frame(m_coefs)
m_coefs <- mutate(m_coefs,
                  coefficient = coefnames,
                  t = Estimate / `Std. Error`,
                  p = pt(abs(t), df = 63, lower.tail = FALSE),
                  lower95 = Estimate - 1.96 * `Std. Error`,
                  upper95 = Estimate + 1.96 * `Std. Error`) %>% select(coefficient, Estimate, `Std. Error`, t, p, lower95, upper95)
write.csv(m_coefs, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.am.spe.csv")
write.csv(aic_a, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.am.spe.aic.csv")

a_full <- lm(am.host.spe ~ species + pH + soil, data = merge_soil)
summary(a_full)
a_resid_spe <- resid(a_full)
a_resid_spe <- as.data.frame(a_resid_spe)
resid_soil <- merge(resid_soil, a_resid_spe, by = 0, all = TRUE)
a_climate <- lm(a_resid_spe ~ climate, data = resid_soil)
summary(a_climate)

a_spe <- ggplot(resid_soil, aes(x = climate, y = a_resid_spe)) + geom_point() + geom_smooth(method = "lm", size = 3) + ylab("") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ")
a_spe


##SAPROBES##
s <- lm(sap.host.spe ~ species + pH + soil + climate, data = merge_soil)
s1 <- lm(sap.host.spe ~ species + pH + climate, data = merge_soil)
s2 <- lm(sap.host.spe ~ species + climate, data = merge_soil)
s3 <- lm(sap.host.spe ~ climate, data = merge_soil)
models_s <- list(s, s1, s2, s3)
aic_s <- aictab(cand.set = models_s)
aic_s
summary(s2)
models_s <- list(s, s1)
model_avg <- summary(model.avg(models_s))
m_coefs <- coefTable(model_avg, full=TRUE, adjust.se = TRUE, df = TRUE)
coefnames <- row.names(m_coefs)
m_coefs <- as.data.frame(m_coefs)
m_coefs <- mutate(m_coefs,
                  coefficient = coefnames,
                  t = Estimate / `Std. Error`,
                  p = pt(abs(t), df = 63, lower.tail = FALSE),
                  lower95 = Estimate - 1.96 * `Std. Error`,
                  upper95 = Estimate + 1.96 * `Std. Error`) %>% select(coefficient, Estimate, `Std. Error`, t, p, lower95, upper95)
write.csv(m_coefs, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.sap.spe.csv")
write.csv(aic_s, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\soil.sap.spe.aic.csv")

s_full <- lm(sap.host.spe ~ species + pH + soil, data = merge_soil)
summary(s_full)
s_resid_spe <- resid(s_full)
s_climate <- lm(s_resid_spe ~ climate, data = merge_soil)
summary(s_climate)
s_resid_spe <- as.data.frame(s_resid_spe)
resid_soil <- cbind(resid_soil, s_resid_spe)

s_spe <- ggplot(resid_soil, aes(x = climate, y = s_resid_spe)) + geom_point() + geom_smooth(method = "lm", size = 3) + ylab("") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ")
s_spe

##THE END OF PT 3: SPECIFICITY##






##FOLIAR DIVERSITY##
##PATHOGENS##
p <- lm(path.diversity ~ species + pH + foliar + climate, data = merge_foliar)
p1 <- lm(path.diversity ~ species + pH + foliar, data = merge_foliar)
p2 <- lm(path.diversity ~ pH + foliar, data = merge_foliar)
p3 <- lm(path.diversity ~ pH, data = merge_foliar)
models_p <- list(p, p1, p2, p3)
aic_p <- aictab(cand.set = models_p)
aic_p
models_p <- list(p1, p3)
summary(p)
model_avg <- summary(model.avg(models_p))
m_coefs <- coefTable(model_avg, full=TRUE, adjust.se = TRUE, df = TRUE)
coefnames <- row.names(m_coefs)
m_coefs <- as.data.frame(m_coefs)
m_coefs <- mutate(m_coefs,
                  coefficient = coefnames,
                  t = Estimate / `Std. Error`,
                  p = pt(abs(t), df = 61, lower.tail = FALSE),
                  lower95 = Estimate - 1.96 * `Std. Error`,
                  upper95 = Estimate + 1.96 * `Std. Error`) %>% select(coefficient, Estimate, `Std. Error`, t, p, lower95, upper95)
write.csv(m_coefs, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar.path.diversity.csv")
write.csv(aic_p, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar.path.diversity.aic.csv")

p_full <- lm(path.diversity ~ species + pH + foliar, data = merge_foliar)
summary(p_full)
p_resid_d <- resid(p_full)
p_climate <- lm(p_resid_d ~ climate, data = merge_foliar)
summary(p_climate)
p_resid_d <- as.data.frame(p_resid_d)
resid_foliar <- cbind(resid_foliar, p_resid_d)

pd <- ggplot(resid_foliar, aes(x = climate, y = p_resid_d)) + geom_point() + geom_smooth(method = "lm", size = 3) + ylab("") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ")
pd


##SAPROBES##
s <- lm(sap.diversity ~ species + pH + foliar + climate, data = merge_foliar)
s1 <- lm(sap.diversity ~ species + pH + climate, data = merge_foliar)
s2 <- lm(sap.diversity ~ species + pH, data = merge_foliar)
s3 <- lm(sap.diversity ~ pH, data = merge_foliar)
models_s <- list(s, s1, s2, s3)
aic_s <- aictab(cand.set = models_s)
aic_s
summary(s2)
model_full <- tidy(s2)
write.csv(model_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar.sap.diversity.csv")
write.csv(aic_s, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar.sap.diversity.aic.csv")

s_full <- lm(sap.diversity ~ species + pH + foliar, data = merge_foliar)
summary(s_full)
s_resid_d <- resid(s_full)
s_climate <- lm(s_resid_d ~ climate, data = merge_foliar)
summary(s_climate)
s_resid_d <- as.data.frame(s_resid_d)
resid_foliar <- cbind(resid_foliar, s_resid_d)

sd <- ggplot(resid_foliar, aes(x = climate, y = s_resid_d)) + geom_point() + geom_smooth(method = "lm", size = 3) + ylab("") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ")
sd
##THE END OF PT1: DIVERSITY##



##FOLIAR CWM##
##PATHOGENS##
p <- lm(pathogen ~ species + pH + foliar + climate, data = merge_foliar)
p1 <- lm(pathogen ~ species + climate + pH, data = merge_foliar)
p2 <- lm(pathogen ~ pH + species, data = merge_foliar)
p3 <- lm(pathogen ~ pH, data = merge_foliar)
models_p <- list(p, p1, p2, p3)
aic_p <- aictab(cand.set = models_p)
aic_p
summary(p)
models_p <- list(p1,p2,p3)
model_avg <- summary(model.avg(models_p))
m_coefs <- coefTable(model_avg, full=TRUE, adjust.se = TRUE, df = TRUE)
coefnames <- row.names(m_coefs)
m_coefs <- as.data.frame(m_coefs)
m_coefs <- mutate(m_coefs,
                  coefficient = coefnames,
                  t = Estimate / `Std. Error`,
                  p = pt(abs(t), df = 61, lower.tail = FALSE),
                  lower95 = Estimate - 1.96 * `Std. Error`,
                  upper95 = Estimate + 1.96 * `Std. Error`) %>% select(coefficient, Estimate, `Std. Error`, t, p, lower95, upper95)
write.csv(m_coefs, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar.path.cwm.csv")
write.csv(aic_p, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar.path.cwm.aic.csv")

p_full <- lm(pathogen ~ species + pH + foliar, data = merge_foliar)
summary(p_full)
p_resid_cwm <- resid(p_full)
p_climate <- lm(p_resid_cwm ~ climate, data = merge_foliar)
summary(p_climate)
p_resid_cwm <- as.data.frame(p_resid_cwm)
resid_foliar <- cbind(resid_foliar, p_resid_cwm)

p_cwm <- ggplot(resid_foliar, aes(x = climate, y = p_resid_cwm)) + geom_point() + geom_smooth(method = "lm", size = 3) + ylab("") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ")
p_cwm


##SAPROBES##
s <- lm(saprobe ~ species + pH + foliar + climate, data = merge_foliar)
s1 <- lm(saprobe ~ species + pH + foliar, data = merge_foliar)
s2 <- lm(saprobe ~ foliar + pH, data = merge_foliar)
s3 <- lm(saprobe ~ foliar, data = merge_foliar)
models_s <- list(s, s1, s2, s3)
aic_s <- aictab(cand.set = models_s)
aic_s
summary(s1)
model_full <- tidy(s1)
write.csv(model_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar.sap.cwm.csv")
write.csv(aic_s, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar.sap.cwm.aic.csv")

s_full <- lm(saprobe ~ species + pH + foliar, data = merge_foliar)
summary(s_full)
s_resid_cwm <- resid(s_full)
s_climate <- lm(s_resid_cwm ~ climate, data = merge_foliar)
summary(s_climate)
s_resid_cwm <- as.data.frame(s_resid_cwm)
resid_foliar <- cbind(resid_foliar, s_resid_cwm)

s_cwm <- ggplot(resid_foliar, aes(x = climate, y = s_resid_cwm)) + geom_point() + geom_smooth(method = "lm", size = 3) + ylab("") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ")
s_cwm

##THE END OF PT 2: CWM##



##FOLIAR SPECIFICITY##
##PATHOGENS##
p <- lm(path.host.spe ~ species + pH + foliar + climate, data = merge_foliar)
p1 <- lm(path.host.spe ~ species + climate + pH, data = merge_foliar)
p2 <- lm(path.host.spe ~ climate + species, data = merge_foliar)
p3 <- lm(path.host.spe ~ climate, data = merge_foliar)
summary(p)
models_p <- list(p, p1, p2, p3)
aic_p <- aictab(cand.set = models_p)
aic_p
models_p <- list(p1,p2)
model_avg <- summary(model.avg(models_p))
m_coefs <- coefTable(model_avg, full=TRUE, adjust.se = TRUE, df = TRUE)
coefnames <- row.names(m_coefs)
m_coefs <- as.data.frame(m_coefs)
m_coefs <- mutate(m_coefs,
                  coefficient = coefnames,
                  t = Estimate / `Std. Error`,
                  p = pt(abs(t), df = 61, lower.tail = FALSE),
                  lower95 = Estimate - 1.96 * `Std. Error`,
                  upper95 = Estimate + 1.96 * `Std. Error`) %>% select(coefficient, Estimate, `Std. Error`, t, p, lower95, upper95)
write.csv(m_coefs, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar.path.spe.csv")
write.csv(aic_p, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar.path.spe.aic.csv")

p_full <- lm(path.host.spe ~ species + pH + foliar, data = merge_foliar)
summary(p_full)
p_resid_spe <- resid(p_full)
p_climate <- lm(p_resid_spe ~ climate, data = merge_foliar)
summary(p_climate)
p_resid_spe <- as.data.frame(p_resid_spe)
resid_foliar <- cbind(resid_foliar, p_resid_spe)

p_spe <- ggplot(resid_foliar, aes(x = climate, y = p_resid_spe)) + geom_point() + geom_smooth(method = "lm", size = 3) + ylab("") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ")
p_spe

##SAPROBES##
s <- lm(sap.host.spe ~ species + pH + foliar + climate, data = merge_foliar)
s1 <- lm(sap.host.spe ~ species + pH + climate, data = merge_foliar)
s2 <- lm(sap.host.spe ~ species + climate, data = merge_foliar)
s3 <- lm(sap.host.spe ~ climate, data = merge_foliar)
models_s <- list(s, s1, s2, s3)
aic_s <- aictab(cand.set = models_s)
aic_s
summary(s2)
model_full <- tidy(s2)
write.csv(model_full, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar.sap.spe.csv")
write.csv(aic_s, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\data\\diversity models\\foliar.sap.spe.aic.csv")

s_full <- lm(sap.host.spe ~ species + pH + foliar, data = merge_foliar)
summary(s_full)
s_resid_spe <- resid(s_full)
s_climate <- lm(s_resid_spe ~ climate, data = merge_foliar)
summary(s_climate)
s_resid_spe <- as.data.frame(s_resid_spe)
resid_foliar <- cbind(resid_foliar, s_resid_spe)

s_spe <- ggplot(resid_foliar, aes(x = climate, y = s_resid_spe)) + geom_point() + geom_smooth(method = "lm", size = 3) + ylab("") + theme_classic(base_size=15) + scale_color_manual(values=c("#993404")) + guides(color = "none") + xlab(" ")
s_spe

##THE END OF PT 3: SPECIFICITY##

##SAVE FILES FOR FIGURES##
saveRDS(resid_foliar, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\resid_foliar_24.rds")
saveRDS(resid_soil, "\\Users\\abbey\\Desktop\\2023-2024\\sgfh manuscript\\rds.files\\resid_soil_24.rds")

##THE END##

