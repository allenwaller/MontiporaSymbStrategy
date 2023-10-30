# Mcap 2019 physiology graphing and analysis
# Script for Allen-Waller & Barott 'Symbiotic dinoflagellates divert energy away from mutualism during coral bleaching recovery' Symbiosis 2023
# Luella Allen-Waller

# Load packages:
library(dplyr)
library(ggplot2)
library(csv)
library(tidyverse)
library(ggpubr)
library(readr)
library(Rmisc)
library(lme4)
library(rstatix)
library(gridExtra)
library(data.table)
library(LambertW)
library(car)
library(patchwork)
library(mgcv)
library(scales)

# Load data
data <- read.csv("Mcap 2019 Phys.csv")
data$delta13C <- data$symb_13C - data$host_13C
data$date.airbrushed <- as.factor(data$date.airbrushed)

# More data, using all the color/PAM data we have, not just from the timepoints
colorpam <- read.csv("Mcap 2019 FvFm and Color.csv")
key <- read.csv("MpHi frag key.csv")
colorpam <- merge(colorpam, key, by = "frag")

# subset for particular timepoints
junecolorpam <- colorpam %>% subset(day == "6/27/19")
julycolorpam <- colorpam %>% subset(day == "7/21/19")
junedata <- data %>% subset(timepoint == "Acute stress")
julydata <- data %>% subset(timepoint == "Initial recovery")

############# Temperature reactions for each timepoint
# FvFm
JunePAMRxn <- ggplot(junecolorpam, aes(x=treatment, y=FvFm, 
                                       color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("FvFm") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(legend.position = "none", panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("") + coord_cartesian(ylim = c(.45,0.75))
JunePAMRxn
JulyPAMRxn <- ggplot(julycolorpam, aes(x=treatment, y=FvFm, 
                                       color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("FvFm") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(legend.position = "none", panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("")+ coord_cartesian(ylim = c(0.45,0.75))
JulyPAMRxn
# Stats:
junePAM.lmm <- lmer(FvFm ~ treatment* hist_symb + (1|parent_colony), data = junecolorpam)
anova(junePAM.lmm)
plot(junePAM.lmm) # ***treatment
julyPAM.lmm <-lmer(FvFm ~ treatment* hist_symb + (1|parent_colony), data = julycolorpam)
anova(julyPAM.lmm)
plot(julyPAM.lmm) # *treatment
# qqs:
qqnorm(resid(junePAM.lmm))
qqline(resid(junePAM.lmm))
qqnorm(resid(julyPAM.lmm))
qqline(resid(julyPAM.lmm))

# chl
JuneChlRxn <- ggplot(junedata, aes(x=treatment, y=totalchl_cm2, 
                                   color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Chl (µg/cm2)") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(legend.position = "none", panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("") + coord_cartesian(ylim = c(4,17.5))
JuneChlRxn
JulyChlRxn <- ggplot(julydata, aes(x=treatment, y=totalchl_cm2, 
                                   color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Chl (µg/cm2)") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(legend.position = "none", panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("") + coord_cartesian(ylim = c(4,17.5))
JulyChlRxn
ChlRxns <- ggarrange(JuneChlRxn, JulyChlRxn)
ChlRxns
# Stats:
junechl.lmm <- lmer(totalchl_cm2 ~ treatment* hist_symb + (1|colony) + (1|date.airbrushed), data = junedata)
summary(junechl.lmm)
anova(junechl.lmm)
plot(junechl.lmm)
julychl.lmm <- lmer(totalchl_cm2 ~ treatment* hist_symb + (1|colony) + (1|date.airbrushed), data = julydata)
anova(julychl.lmm)
plot(julychl.lmm) # *treatment, **treatmentX hist_symb
# qqs:
qqnorm(resid(junechl.lmm))
qqline(resid(junechl.lmm))
qqnorm(resid(julychl.lmm))
qqline(resid(julychl.lmm))

# Symb density
SymbdensJuneRxn <- ggplot(junedata, aes(x=treatment, y=symb_cm2 / 1000 / 1000, 
                                        color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Symbionts (10^3 cells µm-2)") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw() +
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("")
SymbdensJuneRxn
SymbdensJulyRxn <- ggplot(julydata, aes(x=treatment, y=symb_cm2 / 1000 / 1000, 
                                        color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Symbionts (10^3 cells µm-2") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("")
SymbdensJulyRxn
SymbdensRxns <- ggarrange(SymbdensJuneRxn, SymbdensJulyRxn, common.legend = T)
SymbdensRxns
junesymb.lmm <- lmer(symb_cm2 ~ treatment* hist_symb + (1|colony), data = junedata)
anova(junesymb.lmm)
plot(junesymb.lmm)
julysymb.lmm <- lmer(symb_cm2 ~ treatment* hist_symb + (1|colony) + (1|date.airbrushed), data = julydata)
anova(julysymb.lmm)
plot(julysymb.lmm) # ** treatmentX symb
# qqs:
qqnorm(resid(junesymb.lmm))
qqline(resid(junesymb.lmm))
qqnorm(resid(julysymb.lmm))
qqline(resid(julysymb.lmm))

# GP
GPJuneRxn <- ggplot(junedata, aes(x=treatment, y=GP_cm2 *60, 
                                  color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Gross photosynthesis (µmol/hr/cm2)") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(legend.position = "none",panel.grid = element_blank(), axis.title = element_text(size = 12),text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("")
GPJuneRxn
GPJulyRxn <- ggplot(julydata, aes(x=treatment, y=GP_cm2 *60, 
                                  color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Gross photosynthesis (µmol/hr/cm2)") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(legend.position = "none",panel.grid = element_blank(), axis.title = element_text(size = 12),text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("")
GPJulyRxn
GPRxns <- ggarrange(GPJuneRxn, GPJulyRxn)
GPRxns
juneGP.lmm <- lmer(GP_cm2 ~ treatment* hist_symb + (1|colony) + (1|date.airbrushed), data = junedata)
anova(juneGP.lmm)
plot(juneGP.lmm)
julyGP.lmm <- lmer(GP_cm2 ~ treatment* hist_symb + (1|colony) + (1|date.airbrushed), data = julydata)
anova(julyGP.lmm)
plot(julyGP.lmm) # hmm big gap wrt fitted()
# qqs:
qqnorm(resid(juneGP.lmm))
qqline(resid(juneGP.lmm))
qqnorm(resid(julyGP.lmm))
qqline(resid(julyGP.lmm))

# alpha
AlphaJuneRxn <- ggplot(junedata, aes(x=treatment, y=alpha*1000, 
                                     color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("alpha * 1000") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("") + coord_cartesian(ylim = c(0, 0.42))
AlphaJuneRxn
AlphaJulyRxn <- ggplot(julydata, aes(x=treatment, y=alpha *1000, 
                                     color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("alpha*1000") + xlab("") + scale_y_continuous(labels = label_number(accuracy = 0.1))+
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("") + coord_cartesian(ylim = c(0, 0.42))
AlphaJulyRxn
AlphaRxns <- ggarrange(AlphaJuneRxn, AlphaJulyRxn, common.legend = T)
AlphaRxns
junealpha.lmm <- lmer(alpha ~ treatment* hist_symb + (1|colony) + (1|date.airbrushed), data = junedata)
anova(junealpha.lmm)
plot(junealpha.lmm)
julyalpha.lmm <- lmer(alpha ~ treatment* hist_symb + (1|colony)+ (1|date.airbrushed), data = julydata)
anova(julyalpha.lmm)
plot(julyalpha.lmm)
# qqs:
qqnorm(resid(junealpha.lmm))
qqline(resid(junealpha.lmm))
qqnorm(resid(julyalpha.lmm))
qqline(resid(julyalpha.lmm))

# Symbiont 13C
Symb13CJuneRxn <- ggplot(junedata, aes(x=treatment, y=symb_13C, 
                                       color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Symbiont 13C at-%") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(legend.position = "none",panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("") + coord_cartesian(ylim = c(1.105,1.14))
Symb13CJuneRxn
Symb13CJulyRxn <- ggplot(julydata, aes(x=treatment, y=symb_13C, 
                                       color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Symbiont 13C at-%") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(legend.position = "none",panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("") + coord_cartesian(ylim = c(1.105,1.14))
Symb13CJulyRxn
Symb13CRxns <- ggarrange(Symb13CJuneRxn, Symb13CJulyRxn)
Symb13CRxns
junesymb13C.lmm <- lmer(symb_13C ~ treatment* hist_symb + (1|colony) + (1|date.airbrushed), data = junedata)
anova(junesymb13C.lmm)
plot(junesymb13C.lmm)
julysymb13C.lmm <- lmer(symb_13C ~ treatment* hist_symb + (1|colony) + (1|date.airbrushed), data = julydata)
anova(julysymb13C.lmm)
plot(julysymb13C.lmm)
# qqs:
qqnorm(resid(junesymb13C.lmm))
qqline(resid(junesymb13C.lmm))
qqnorm(resid(julysymb13C.lmm))
qqline(resid(julysymb13C.lmm))

# Host protein
ProtJuneRxn <- ggplot(junedata, aes(x=treatment, y=prot_cm2/1000, 
                                    color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Host protein (mg/cm^2)") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("") + coord_cartesian(ylim = c(0,0.75))
ProtJuneRxn
ProtJulyRxn <- ggplot(julydata, aes(x=treatment, y=prot_cm2 /1000, 
                                    color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Host protein (mg/cm^2)") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("") + coord_cartesian(ylim = c(0,0.7))
ProtJulyRxn
ProtRxns <- ggarrange(ProtJuneRxn, ProtJulyRxn, common.legend = T)
ProtRxns
juneprot.lmm <- lmer(prot_cm2 ~ treatment* hist_symb + (1|colony) + (1|date.airbrushed), data = junedata)
anova(juneprot.lmm)
plot(juneprot.lmm)
julyprot.lmm <- lmer(prot_cm2 ~ treatment* hist_symb + (1|colony) + (1|date.airbrushed), data = julydata)
anova(julyprot.lmm)
plot(julyprot.lmm)
# qqs:
qqnorm(resid(juneprot.lmm))
qqline(resid(juneprot.lmm))
qqnorm(resid(julyprot.lmm))
qqline(resid(julyprot.lmm))

# Lipids
LipJuneRxn <- ggplot(junedata, aes(x=treatment, y=lipids.mg.cm2, 
                                   color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Host lipids (mg/cm^2)") + xlab("") + scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("") + coord_cartesian(ylim = c(0,1))
LipJuneRxn
LipJulyRxn <- ggplot(julydata, aes(x=treatment, y=lipids.mg.cm2, 
                                   color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Host lipids (mg/cm^2)") + xlab("") + scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("") + coord_cartesian(ylim =c(0,.5))
LipJulyRxn
LipRxns <- ggarrange(LipJuneRxn, LipJulyRxn, common.legend = T)
LipRxns
junelip.lmm <- lmer(lipids.ug.cm2 ~ treatment* hist_symb + (1|colony) + (1|date.airbrushed), data = junedata)
anova(junelip.lmm)
plot(junelip.lmm)
julylip.lmm <- lmer(lipids.ug.cm2 ~ treatment* hist_symb + (1|colony) + (1|date.airbrushed), data = julydata)
anova(julylip.lmm) # nothing
plot(julylip.lmm)
# qqs:
qqnorm(resid(junelip.lmm))
qqline(resid(junelip.lmm))
qqnorm(resid(julylip.lmm))
qqline(resid(julylip.lmm))

# calcif.pg.d.cm^2ication
CalcJuneRxn <- ggplot(junedata, aes(x=treatment, y=calcif.pg.d.cm^2*1000, 
                                    color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("calcif.pg.d.cm^2ication (ng/day/cm^2)") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("")+ coord_cartesian(ylim =c(1.5,4.6))
CalcJuneRxn
CalcJulyRxn <- ggplot(julydata, aes(x=treatment, y=calcif.pg.d.cm^2*1000, 
                                    color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("calcif.pg.d.cm^2ication (ng/day/cm^2)") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("") + coord_cartesian(ylim =c(1.5,3))
CalcJulyRxn
CalcRxns <- ggarrange(CalcJuneRxn, CalcJulyRxn, common.legend = T)
CalcRxns
junecalc.lmm <- lmer(calcif.pg.d.cm^2 ~ treatment* hist_symb + (1|colony), data = junedata)
anova(junecalc.lmm)
plot(junecalc.lmm) # treatment
julycalc.lmm <- lmer(calcif.pg.d.cm^2 ~ treatment* hist_symb + (1|colony), data = julydata)
anova(julycalc.lmm)
plot(julycalc.lmm)
# qqs:
qqnorm(resid(junecalc.lmm))
qqline(resid(junecalc.lmm))
qqnorm(resid(julycalc.lmm))
qqline(resid(julycalc.lmm))

# P:R
JunePR <- ggplot(junedata, aes(x=treatment, y=P.R, 
                               color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("P:R") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) + ggtitle("")+
  coord_cartesian(ylim = c(4,8))
JunePR
JulyPR <- ggplot(julydata, aes(x=treatment, y=P.R, 
                               color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("P:R") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("") +coord_cartesian(ylim = c(4,8))
JulyPR
PRRxns <- ggarrange(JunePR, JulyPR, common.legend = T)
PRRxns
# stats
junePR.lmm <- lmer(P.R ~ treatment* hist_symb + (1|colony), data = junedata)
anova(junePR.lmm) # borderline effect of treatment0.05472
plot(junePR.lmm)
julyPR.lmm <- lmer(P.R ~ treatment* hist_symb + (1|colony), data = julydata)
anova(julyPR.lmm)
plot(julyPR.lmm)
# qqs:
qqnorm(resid(junePR.lmm))
qqline(resid(junePR.lmm))
qqnorm(resid(julyPR.lmm))
qqline(resid(julyPR.lmm))

# Host 13C
Host13CJuneRxn <- ggplot(junedata, aes(x=treatment, y=host_13C, 
                                       color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Host 13C at-%") + xlab("") + scale_y_continuous(labels(label_number(accuracy = 0.01))) +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("") + coord_cartesian(ylim = c(1.092,1.107))
Host13CJuneRxn
Host13CJulyRxn <- ggplot(julydata, aes(x=treatment, y=host_13C, 
                                       color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Host 13C at-%") + xlab("") + scale_y_continuous(labels(label_number(accuracy = 0.01))) +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("") + coord_cartesian(ylim = c(1.092,1.107))
Host13CJulyRxn
Host13CRxns <- ggarrange(Host13CJuneRxn, Host13CJulyRxn, common.legend = T)
Host13CRxns
junehost13C.lmm <- lmer(host_13C ~ treatment* hist_symb + (1|colony) + (1|date.airbrushed), data = junedata)
anova(junehost13C.lmm)
plot(junehost13C.lmm)
julyhost13C.lmm <- lmer(host_13C ~ treatment* hist_symb + (1|colony)  + (1|date.airbrushed), data = julydata)
anova(julyhost13C.lmm)
plot(julyhost13C.lmm)
# qqs:
qqnorm(resid(junehost13C.lmm))
qqline(resid(junehost13C.lmm))
qqnorm(resid(julyhost13C.lmm))
qqline(resid(julyhost13C.lmm))

### Compilations (aspect ratio for .pdf export 5.5:10)
# Acute organismal response
AcuteOrganismal <- ggarrange(JunePAMRxn,SymbdensJuneRxn, JuneChlRxn,GPJuneRxn,AlphaJuneRxn, Symb13CJuneRxn, JunePR,ProtJuneRxn,LipJuneRxn,CalcJuneRxn, Host13CJuneRxn,
                             ncol = 6, nrow = 2, widths = 1.1, common.legend = T)
AcuteOrganismal
# Initial recovery organismal response
InitRecovOrganismal <- ggarrange(JulyPAMRxn,SymbdensJulyRxn, JulyChlRxn, GPJulyRxn,AlphaJulyRxn,Symb13CJulyRxn,JulyPR,ProtJulyRxn,LipJulyRxn,CalcJulyRxn, Host13CJulyRxn,
                                  ncol = 6, nrow = 2, widths = 1.1, common.legend = T)
InitRecovOrganismal


########## Recovery trajectories for different temperatures
# Acute to initial-recovery trajectories for each variable
hot <- data %>% subset(treatment== "31°")
ambient <- data %>% subset(treatment== "28°")
hotcolorpam <- colorpam %>% subset(treatment == "31°") %>% 
  subset(day != "7/8/19")
ambcolorpam <- colorpam %>% subset(treatment == "28°") %>%
  subset(day != "7/8/19")

# FvFm
FvFmRecovery <- ggplot(hotcolorpam, aes(x=day, y=FvFm, 
                                        color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Photosynthetic efficiency (Fv/Fm)") + xlab("") + scale_y_continuous(labels = label_number(accuracy = 0.1))+
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  scale_x_discrete(labels = c("acute", "subacute")) +
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("") + coord_cartesian(ylim = c(.4,.7))
FvFmRecovery
# lmm
fvfmrecov.lmm <- lmer(FvFm ~ hist_symb * day + (1|parent_colony), data = hotcolorpam)
qqnorm(resid(fvfmrecov.lmm))
qqline(resid(fvfmrecov.lmm))
plot(fvfmrecov.lmm) # heteroscedastic; heavy-tailed
# Try LambertW Guassianize:
# first subset out frags for which we don't have FvFm data that time:
hotcolorpam.lam <- subset(hotcolorpam, frag!="52" & frag!="57")
hotcolorpam.lam$lambPAM <- hotcolorpam.lam$FvFm %>% Gaussianize()
hotlambPAM.lmm <- lmer(lambPAM ~ day * hist_symb + (1|parent_colony), data = hotcolorpam.lam)
plot(hotlambPAM.lmm) # much better residual plot!
qqnorm(resid(hotlambPAM.lmm))
qqline(resid(hotlambPAM.lmm))
anova(hotlambPAM.lmm)
# significant effect of day X hist_symb (*)

# Symbiont counts
SymbRecovery <- ggplot(hot, aes(x=timepoint, y=symb_cm2 / 1000 / 1000, 
                                color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Symbiont density (10^3 cells/µm2)") + xlab("") +scale_y_continuous(labels = label_number(accuracy = 0.1))+
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("")
SymbRecovery
# lm
symbrecov.lmm <- lmer(symb_cm2 ~ hist_symb * timepoint + (1|colony) + (1|date.airbrushed), data = hot)
anova(symbrecov.lmm)
qqnorm(resid(symbrecov.lmm))
plot(symbrecov.lmm)

# Chl
ChlRecovery <- ggplot(hot, aes(x=timepoint, y=totalchl_cm2, 
                               color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Chl (µg/cm2)") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("")
ChlRecovery
# lm
chlrecov.lmm <- lmer(totalchl_cm2 ~ hist_symb * timepoint + (1|colony) + (1|date.airbrushed), data = hot)
anova(chlrecov.lmm)
qqnorm(resid(chlrecov.lmm))
qqline(resid(chlrecov.lmm))
plot(chlrecov.lmm)
# Try LambertW Guassianize:
hot$lambchl <- hot$totalchl_cm2 %>% Gaussianize()
hotlambchl.lmm <- lmer(lambchl ~ timepoint * hist_symb + (1|colony), data = hot)
plot(hotlambchl.lmm) # better residual plot
qqnorm(resid(hotlambchl.lmm)) # closer to normal
qqline(resid(hotlambchl.lmm))
anova(hotlambPAM.lmm)

# alpha
hot$alphax1000 <- hot$alpha * 1000
AlphaRecovery <- ggplot(hot, aes(x=timepoint, y=alphax1000, 
                                 color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("alpha x 10^3") + xlab("") + scale_y_continuous(labels = label_number(accuracy = 0.1))+
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("") + coord_cartesian(ylim = c(0,0.45))
AlphaRecovery
# lm
AlphaRecovery.lmm <- lmer(alpha ~ hist_symb * timepoint + (1|colony) + (1|date.airbrushed), data = hot)
anova(AlphaRecovery.lmm)
qqnorm(resid(AlphaRecovery.lmm))
qqline(resid(AlphaRecovery.lmm))
plot(AlphaRecovery.lmm)

# GP
GPRecovery <- ggplot(hot, aes(x=timepoint, y=GP_cm2 * 60, 
                              color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("GP (µmol O2/hr/cm2)") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("")
GPRecovery
# lm
GPRecovery.lmm <- lmer(GP_cm2 ~ hist_symb * timepoint + (1|colony) + (1|date.airbrushed), data = hot)
anova(GPRecovery.lmm)
qqnorm(resid(GPRecovery.lmm))
qqline(resid(GPRecovery.lmm))
plot(GPRecovery.lmm)

# Symbiont 13C
Symb13CRecovery <- ggplot(hot, aes(x=timepoint, y=symb_13C, 
                                   color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Symbiont 13C at-%") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("")
Symb13CRecovery
# lm
Symb13CRecovery.lmm <- lmer(symb_13C ~ hist_symb * timepoint + (1|colony) + (1|date.airbrushed), data = hot)
anova(Symb13CRecovery.lmm)
qqnorm(resid(Symb13CRecovery.lmm))
qqline(resid(Symb13CRecovery.lmm))
plot(Symb13CRecovery.lmm)

# Host 13C
Host13CRecovery <- ggplot(hot, aes(x=timepoint, y=host_13C, 
                                   color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Host 13C at-%") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("")
Host13CRecovery
# lm
Host13CRecovery.lmm <- lmer(host_13C ~ hist_symb * timepoint + (1|colony) + (1|date.airbrushed), data = ambient)
anova(Host13CRecovery.lmm)
qqnorm(resid(Host13CRecovery.lmm))
qqline(resid(Host13CRecovery.lmm))
plot(Host13CRecovery.lmm)

# Protein
HostProteinRecovery <- ggplot(hot, aes(x=timepoint, y=prot_cm2/1000, 
                                       color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Protein (mg/cm2)") + xlab("") + scale_y_continuous(labels = label_number(accuracy = 0.1))+
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("") + coord_cartesian(ylim = c(0,0.6))
HostProteinRecovery
# lmm
HostProteinRecovery.lmm <- lmer(prot_cm2 ~ hist_symb * timepoint + (1|colony) + (1|date.airbrushed), data = hot)
anova(HostProteinRecovery.lmm)
qqnorm(resid(HostProteinRecovery.lmm))
qqline(resid(HostProteinRecovery.lmm))
plot(HostProteinRecovery.lmm)

# calcif.pg.d.cm^2ication
CalcRecovery <- ggplot(hot, aes(x=timepoint, y=calcif.pg.d.cm^2 *1000, 
                                color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("calcif.pg.d.cm^2ication (ng/day/cm2)") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("")
CalcRecovery
# lm
CalcRecovery.lmm <- lmer(calcif.pg.d.cm^2 ~ hist_symb * timepoint + (1|colony), data = hot)
anova(CalcRecovery.lmm)
qqnorm(resid(CalcRecovery.lmm))
qqline(resid(CalcRecovery.lmm))
plot(CalcRecovery.lmm)

# P:R
PRRecov <- ggplot(hot, aes(x=timepoint, y=P.R, 
                           color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("P:R") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("") + coord_cartesian(ylim = c(4,8))
PRRecov
# lmm
PRRecov.lmm <- lmer(P.R ~ hist_symb * timepoint + (1|colony) + (1|date.airbrushed), data = hot)
anova(PRRecov.lmm)
qqnorm(resid(PRRecov.lmm))
qqline(resid(PRRecov.lmm))
plot(PRRecov.lmm)

# Lipid recovery
LipidRecov <- ggplot(hot, aes(x=timepoint, y=lipids.mg.cm2, 
                              color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Host lipids (mg/cm2)") + xlab("") + scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("")+ coord_cartesian(ylim = c(0,1))
LipidRecov
LipidRecov.lmm <- lmer(lipids.mg.cm2 ~ hist_symb * timepoint + (1|colony), data = hot)
anova(LipidRecov.lmm)
qqnorm(resid(LipidRecov.lmm))
qqline(resid(LipidRecov.lmm))
plot(LipidRecov.lmm)

# Heated recovery compilation
RecoveryOrganismalPhys <- ggarrange(FvFmRecovery,SymbRecovery,ChlRecovery,GPRecovery,AlphaRecovery,Symb13CRecovery,
                                    PRRecov,HostProteinRecovery, LipidRecov, CalcRecovery,Host13CRecovery,
                                    ncol = 6, nrow = 2, widths = 1.1, common.legend = T)
RecoveryOrganismalPhys

# Ambient only
ambient <- data %>% subset(treatment== "28º")
ambientcolorpam <- colorpam %>% subset(treatment == "28º") %>% 
  subset(day != "7/8/19")
ambient$group <- as.factor(paste(ambient$hist_symb, ambient$timepoint))

# FvFm
AmbFvFmRecovery <- ggplot(ambientcolorpam, aes(x=day, y=FvFm, 
                                               color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Photosynthetic efficiency (Fv/Fm)") + xlab("") + scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  scale_x_discrete(labels = c("acute", "lag")) +
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("") + coord_cartesian(ylim = c(0.5,0.8))
AmbFvFmRecovery
# lm
ambfvfmrecov.lmm <- lmer(FvFm ~ hist_symb * day + (1|parent_colony) + (1|initial_tank), data = ambientcolorpam)
anova(ambfvfmrecov.lmm) # ***timepoint
qqnorm(resid(ambfvfmrecov.lmm))
qqline(resid(ambfvfmrecov.lmm))
plot(ambfvfmrecov.lmm)

# Symbiont counts
AmbSymbRecovery <- ggplot(ambient, aes(x=timepoint, y=symb_cm2 / 1000 / 1000, 
                                       color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Symbiont density (10^3 cells/µm2)") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("")
AmbSymbRecovery
# lmm
ambsymbrecov.lmm <- lmer(symb_cm2 ~ hist_symb * timepoint + (1|colony), data = ambient)
anova(ambsymbrecov.lmm)
qqnorm(resid(ambsymbrecov.lmm))
qqline(resid(ambsymbrecov.lmm))
plot(ambsymbrecov.lmm)

# Chl
AmbChlRecovery <- ggplot(ambient, aes(x=timepoint, y=totalchl_cm2, 
                                      color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Chl (µg/cm2)") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("")
AmbChlRecovery
# lm
chlrecov.lmm <- lmer(totalchl_cm2 ~ hist_symb * timepoint + (1|colony), data = ambient)
anova(chlrecov.lmm)
qqnorm(resid(chlrecov.lmm))
qqline(resid(chlrecov.lmm))
plot(chlrecov.lmm)

# alpha
AmbAlphaRecovery <- ggplot(ambient, aes(x=timepoint, y=alpha*1000, 
                                        color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("alpha x 10^3") + xlab("") + scale_y_continuous(labels = label_number(accuracy = 0.1))+
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("") + coord_cartesian(ylim = c(0,0.5))
AmbAlphaRecovery
# lm
AlphaRecovery.lmm <- lmer(alpha ~ hist_symb * timepoint + (1|colony), data = ambient)
anova(AlphaRecovery.lmm)
qqnorm(resid(AlphaRecovery.lmm))
qqline(resid(AlphaRecovery.lmm))
plot(AlphaRecovery.lmm)

# GP
AmbGPRecovery <- ggplot(ambient, aes(x=timepoint, y=GP_cm2 *60, 
                                     color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("GP (µmol O2/min/cm2)") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("")
AmbGPRecovery
# lm
GPRecovery.lmm <- lmer(GP_cm2 ~ hist_symb * timepoint + (1|colony), data = ambient)
anova(GPRecovery.lmm)
qqnorm(resid(GPRecovery.lmm))
qqline(resid(GPRecovery.lmm))
plot(GPRecovery.lmm)

# Symbiont 13C
AmbSymb13CRecovery <- ggplot(ambient, aes(x=timepoint, y=symb_13C, 
                                          color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Symbiont 13C at-%") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("")
AmbSymb13CRecovery
# lm
AmbSymb13CRecovery.lmm <- lmer(symb_13C ~ hist_symb * timepoint + (1|colony), data = ambient)
anova(AmbSymb13CRecovery.lmm)
qqnorm(resid(AmbSymb13CRecovery.lmm))
qqline(resid(AmbSymb13CRecovery.lmm))
plot(AmbSymb13CRecovery.lmm)

# Host 13C
AmbHost13CRecovery <- ggplot(ambient, aes(x=timepoint, y=host_13C, 
                                          color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Host 13C at-%") + xlab("") + scale_y_continuous(labels = label_number(accuracy = 0.001))+
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("") + coord_cartesian(ylim = c(1.09,1.110))
AmbHost13CRecovery
# lm
AmbHost13CRecovery.lmm <- lmer(host_13C ~ hist_symb * timepoint + (1|colony), data = ambient)
anova(AmbHost13CRecovery.lmm)
qqnorm(resid(AmbHost13CRecovery.lmm))
qqline(resid(AmbHost13CRecovery.lmm))

# Protein
AmbHostProteinRecovery <- ggplot(ambient, aes(x=timepoint, y=prot_cm2 /1000, 
                                              color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Protein (mg/cm2)") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("") + coord_cartesian(ylim = c(0,0.8))
AmbHostProteinRecovery
# lm
HostProteinRecovery.lmm <- lmer(prot_cm2 ~ hist_symb * timepoint + (1|colony), data = ambient)
anova(HostProteinRecovery.lmm)
qqnorm(resid(HostProteinRecovery.lmm))
qqline(resid(HostProteinRecovery.lmm))
plot(HostProteinRecovery.lmm)

# calcif.pg.d.cm^2ication
AmbCalcRecovery <- ggplot(ambient, aes(x=timepoint, y=calcif.pg.d.cm^2 * 1000, 
                                       color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("calcif.pg.d.cm^2ication (ng/day/cm2)") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("")
AmbCalcRecovery
# lm
CalcRecovery.lmm <- lmer(calcif.pg.d.cm^2 ~ hist_symb * timepoint + (1|colony), data = ambient)
anova(CalcRecovery.lmm)
qqnorm(resid(CalcRecovery.lmm))
qqline(resid(CalcRecovery.lmm))
plot(CalcRecovery.lmm)

# PR
AmbPRRecov <- ggplot(ambient, aes(x=timepoint, y=P.R, 
                                  color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("P:R") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("")
AmbPRRecov
# lmm
AmbPRRecov.lmm <- lmer(P.R ~ hist_symb * timepoint + (1|colony), data = ambient)
anova(AmbPRRecov.lmm)
qqnorm(resid(AmbPRRecov.lmm))
qqline(resid(AmbPRRecov.lmm))
plot(AmbPRRecov.lmm)

# Lipids
AmbLipRecov <- ggplot(ambient, aes(x=timepoint, y=lipids.mg.cm2, 
                                   color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Lipids (mg/cm2)") + xlab("") + scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("")+ coord_cartesian(ylim = c(0,0.7))
AmbLipRecov
# lmm
AmbLipRecov.lmm <- lmer(lipids.mg.cm2 ~ hist_symb * timepoint + (1|colony), data = ambient)
qqnorm(resid(AmbLipRecov.lmm))
qqline(resid(AmbLipRecov.lmm)) # very heavy upper tail
plot(AmbLipRecov.lmm)
# Lambert's W:
lambamb <- ambient %>% drop_na(lipids.mg.cm2) 
lambamb$LambLip <- lambamb$lipids.mg.cm2 %>% Gaussianize()
LambAmbLipRecov.lmm <- lmer(LambLip ~ hist_symb * timepoint + (1|colony), data = lambamb)
qqnorm(resid(LambAmbLipRecov.lmm))
qqline(resid(LambAmbLipRecov.lmm)) # a bit better
plot(resid(LambAmbLipRecov.lmm))
anova(LambAmbLipRecov.lmm)

# Ambient compilation
AmbOrganismalResponse <- ggarrange(AmbFvFmRecovery,AmbSymbRecovery, AmbChlRecovery, AmbGPRecovery,AmbAlphaRecovery,AmbSymb13CRecovery,
                                   AmbPRRecov,AmbHostProteinRecovery,AmbLipRecov,AmbCalcRecovery, AmbHost13CRecovery,
                                   ncol = 6, nrow = 2, widths = 1.1, common.legend = T)
AmbOrganismalResponse

######## Effect of photosynthate provisioning on intracellular pH and calcif.pg.d.cm^2ication
# basal pHi
SymbBasalpHiJulyRxn <- ggplot(julydata, aes(x=treatment, y=symb_pHi_basal, 
                                            color=hist_symb)) +
  stat_summary(fun = mean, geom = "path", aes(group = hist_symb)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(width=.1, group=hist_symb)) +
  ylab("Symbiocyte intracellular pH") + xlab("") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14), legend.key = element_blank(), plot.title = element_text(hjust = 0.5)) +
  ggtitle("") +coord_cartesian(ylim = c(7.3,8.2))
SymbBasalpHiJulyRxn
# stats:
SymbBasalpHiJuly.lmm <- lmer(symb_pHi_basal ~ treatment* hist_symb + (1|colony), data = julydata)
anova(SymbBasalpHiJuly.lmm)
qqline(resid(SymbBasalpHiJuly.lmm))
qqnorm(resid(SymbBasalpHiJuly.lmm))

# symbiocyte basal pHi by host carbon assimilation:
host13C.vs.symbpHi.lm <- lm(symb_pHi_basal ~ host_13C, data = data)
summary(host13C.vs.symbpHi.lm)
coef(host13C.vs.symbpHi.lm)
# p = 0.0339
host13C.vs.symbpHi <- ggplot(data, aes(x = host_13C, y = symb_pHi_basal)) +
  geom_point(size = 2, aes(shape = treatment, color = hist_symb, group = interaction(hist_symb, treatment))) +
  ylab("Intracellular pH") + 
  xlab("Host 13C at-%") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  scale_shape_manual(values = c(2, 1, 17, 16))+
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14)) +
  labs(title = "") + coord_cartesian(xlim = c(1.0915,1.104))+
  geom_smooth(color = "black", size = 0.7,method = lm) +
  labs(color = "Natural history", shape = "treatment")
host13C.vs.symbpHi
host13C.vs.symbpHi.lm <- cor.test(data$host_13C, data$symb_pHi_basal)
host13C.vs.symbpHi.lm
# p=0.03392, R2 = 0.2830

# calcif.pg.d.cm^2ication by host carbon assimilation:
host13C.vs.calcif.pg.d.cm^2 <- ggplot(data, aes(x = host_13C, y = calcif.pg.d.cm^2 * 1000)) +
  geom_point(size = 2, aes(color = hist_symb, 
                           shape = treatment)) +
  ylab("calcif.pg.d.cm^2ication (ng cm-2 day-1)") + 
  xlab("Host 13C at-%") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  scale_shape_manual(values = c(2, 1, 17, 16))+
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14)) +
  labs(title = "") +
  labs(color = "Natural history", shape = "treatment") +
  geom_smooth(color = "black", size = 0.7, method=lm)
host13C.vs.calcif.pg.d.cm^2
host13C.vs.calcif.pg.d.cm^2.lm <- cor.test(data$host_13C, data$calcif.pg.d.cm^2)
host13C.vs.calcif.pg.d.cm^2.lm
# t = 2.4268, df = 38, p-value = 0.02008, cor = 0.3663127

# Host carbon assim. effects on physiology compilation figure (Fig 3):
pHi.Calc.Figure <- ggarrange(SymbBasalpHiJulyRxn,host13C.vs.symbpHi,host13C.vs.calcif.pg.d.cm^2,
                             nrow = 1, ncol = 3,
                             widths = c(1.5,1.5,1.5), common.legend = T)
pHi.Calc.Figure

######## Effect of symbiont density on carbon translocation
# host 13C by #symbs
symbperprot.vs.host13c <- ggplot(data, aes(x = symb_nprot, y = host_13C)) +
  geom_point(size = 2, aes(color = hist_symb, 
                           shape = treatment)) +
  ylab("Host 13C at-%") + 
  xlab("Symbiont density (cells/ng host protein)") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  scale_shape_manual(values = c(2, 1, 17, 16))+
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14)) +
  labs(title = "") + 
  labs(color = "Natural history", shape = "treatment") +
  geom_smooth(color = "black", method=lm, size = 0.7, fullrange=T)
symbperprot.vs.host13c
symbperprot.vs.host13c.lm <- cor.test(data$symb_nprot, data$host_13C)
symbperprot.vs.host13c.lm
# simple linear corr: t = 1.075, df = 38, R2 = 0.0295 (p=0.2891)
symbperprot.vs.host13c.lmm <- lmer(host_13C ~ symb_nprot * treatment* hist_symb + (1|colony), data = data)
anova(symbperprot.vs.host13c.lmm)
plot(symbperprot.vs.host13c.lmm)
# **treatmentand *symb species but no density effects or interactions

# delta13C by #symbs
symbperprot.vs.delta13c.lm <- lm(delta13C ~ symb_nprot, data = data)
summary(symbperprot.vs.delta13c.lm)
# *p = 0.016781
symbperprot.vs.delta13c <- ggplot(data, aes(x = symb_nprot, y = delta13C)) +
  geom_point(size = 2, aes(color = hist_symb, 
                           shape = treatment)) +
  ylab("∆13C at% (symb-host)") + 
  xlab("Symbiont density (cells/ng host protein)") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  scale_shape_manual(values = c(2, 1, 17, 16))+
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14)) +
  labs(title = "") + 
  labs(color = "Natural history", shape = "treatment") +
  coord_cartesian(ylim=c(0.002,0.045)) +
  geom_smooth(method=lm, fullrange=T, size = 0.7, color = "black")
symbperprot.vs.delta13c
symbperprot.vs.delta13c.lm <- cor.test(data$symb_ugprot, data$delta13C)
symbperprot.vs.delta13c.lm
# simple linear corr: t = 2.5018, df = 38, R2 = 0.1414, p-value = 0.01678

# Symbiont population and host 13C compilation figure:
SymbPopSelfishness <- ggarrange(symbperprot.vs.host13c,symbperprot.vs.delta13c,
                                common.legend = T)
SymbPopSelfishness


########## Correlation plots for supplemental figures

# GP/symb vs. pHi
GPpersymb.vs.pHi <- ggplot(data, aes(x = GP_persymb * 60 * 1000, y = symb_pHi_basal)) +
  geom_point(size = 2, aes(color = hist_symb, 
                           shape = treatment)) +
  ylab("Intracellular pH (pHi)") + 
  xlab("Photosynthesis per symbiont (umol/10^3 cells/hr)") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  scale_shape_manual(values = c(2, 1, 17, 16))+
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14)) +
  labs(title = "") + 
  labs(color = "Natural history", shape = "treatment") +
  coord_cartesian(xlim = c(0,0.02)) +
  geom_smooth(color = "black", method=lm, size = 0.7)
GPpersymb.vs.pHi
GPpersymb.vs.pHi.lm <- cor.test(data$GP_persymb, data$symb_pHi_basal)
GPpersymb.vs.pHi.lm
# t = 0.84295, df = 14, p-value = 0.4134, cor = 0.21977

# calcif.pg.d.cm^2ication vs. pHi
calcif.pg.d.cm^2.vs.pHi <- ggplot(data, aes(x = calcif.pg.d.cm^2 * 1000, y = symb_pHi_basal)) +
  geom_point(size = 2, aes(color = hist_symb, 
                           shape = treatment)) +
  ylab("Intracellular pH (pHi)") + 
  xlab("calcif.pg.d.cm^2ication (ng cm-2 day-1)") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  scale_shape_manual(values = c(2, 1, 17, 16))+
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14)) +
  labs(title = "") + coord_cartesian(xlim = c(0,4)) +
  labs(color = "Natural history", shape = "treatment") +
  geom_smooth(color = "black", size = 0.7,  method=lm)
calcif.pg.d.cm^2.vs.pHi
calcif.pg.d.cm^2.vs.pHi.lm <- cor.test(data$calcif.pg.d.cm^2, data$symb_pHi_basal)
calcif.pg.d.cm^2.vs.pHi.lm
# t = 0.26914, df = 14, p-value = 0.7917, cor = 0.0717

# Symbiont selfishness, density, and respiration (Fig S8)
# symb per cm2 vs selfishness:
symbpercm2.vs.delta13C <- ggplot(data, aes(x = symb_cm2 / 1000 / 1000, y = delta13C)) +
  geom_point(size = 2, aes(color = hist_symb, 
                           shape = treatment)) +
  ylab("Symbiont selfishness (Delta13C)") + 
  xlab("Symbiont density (10^3 cells µm-2)") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  scale_shape_manual(values = c(2, 1, 17, 16))+
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14)) +
  labs(title = "") +
  labs(color = "Natural history", shape = "treatment") +
  geom_smooth(color = "black", method=lm, size=0.7)
symbpercm2.vs.delta13C
symbpercm2.vs.delta13C.lm <- cor.test(data$symb_cm2, data$delta13C)
symbpercm2.vs.delta13C.lm
# t = 2.6087, df = 38, p-value = 0.01292, R2 = 0.15188
# respiration vs. selfishness:
data$R_perngprot <- data$R_cm2 * 1000 / data$prot_cm2
resp.vs.delta13C <- ggplot(data, aes(x = R_perngprot * 60, y = delta13C)) +
  geom_point(size = 2, aes(color = hist_symb, 
                           shape = treatment)) +
  ylab("Symbiont selfishness (Delta13C)") + 
  xlab("Respiration (umol O2/ng host protein/hr)") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  scale_shape_manual(values = c(2, 1, 17, 16))+
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14)) +
  labs(title = "") +
  labs(color = "Natural history", shape = "treatment") +
  geom_smooth(color = "black", method=lm, size=0.7)
resp.vs.delta13C
resp.vs.delta13C.lm <- cor.test(data$R_perngprot, data$delta13C)
resp.vs.delta13C.lm
# t = 0.65453, df = 37, p-value = 0.5168, R2 = 0.01145
# both selfishness panels:
Selfishness.supp <- ggarrange(symbpercm2.vs.delta13C, resp.vs.delta13C, common.legend = T)
Selfishness.supp

# Symbiont density is related to respiration and chlorophyll but not to lipids or calcif.pg.d.cm^2ication - Fig S9
# symbdens vs. calcif.pg.d.cm^2ication:
symbpercm2.vs.calc <- ggplot(data, aes(x = symb_cm2 / 1000, y = calcif.pg.d.cm^2 * 1000)) +
  geom_point(size = 2, aes(color = hist_symb, 
                           shape = treatment)) +
  ylab("calcif.pg.d.cm^2ication (ng cm-2 day-1)") + 
  xlab("Symbiont density (cells µm-2)") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  scale_shape_manual(values = c(2, 1, 17, 16))+
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14)) +
  labs(title = "") +
  labs(color = "Natural history", shape = "treatment") +
  geom_smooth(color = "black", method=lm, size=0.7)
symbpercm2.vs.calc
symbpercm2.vs.calc.lm <- cor.test(data$symb_cm2, data$calcif.pg.d.cm^2)
symbpercm2.vs.calc.lm
# t = 1.9113, df = 38, p-value = 0.06353, cor = 0.2961487
# symbdens vs. chl/symb:
symbpercm2.vs.chlsymb <- ggplot(data, aes(x = symb_cm2 / 1000, y = totalchl_symb * 1000)) +
  geom_point(size = 2, aes(color = hist_symb, 
                           shape = treatment)) +
  ylab("Chlorophyll (mg/symbiont)") + 
  xlab("Symbiont density (cells µm-2)") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  scale_shape_manual(values = c(2, 1, 17, 16))+
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14)) +
  labs(title = "") +
  labs(color = "Natural history", shape = "treatment") +
  geom_smooth(color = "black", method=loess, size=0.7)
symbpercm2.vs.chlsymb
symbpercm2.vs.chlsymb.gam <- gam(symb_cm2 ~ s(symb_cm2), method = "REML", data = data)
symbpercm2.vs.chlsymb.gam
anova(symbpercm2.vs.chlsymb.gam)
# edf = 1.06, F = 1.001e+32, p<0.0001, REML score = -777.8559
# symbdens vs. respiration/cm2:
symbpercm2.vs.resp <- ggplot(data, aes(x = symb_cm2 / 1000, y = R_cm2 *60)) +
  geom_point(size = 2, aes(color = hist_symb, 
                           shape = treatment)) +
  ylab("Respiration (umol O2 cm-2 hr-1)") + 
  xlab("Symbiont density (cells µm-2)") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  scale_shape_manual(values = c(2, 1, 17, 16))+
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14)) +
  labs(title = "") +
  labs(color = "Natural history", shape = "treatment") +
  geom_smooth(color = "black", method=lm, size=0.7)
symbpercm2.vs.resp
symbpercm2.vs.resp.lm <- cor.test(data$symb_cm2, data$R_cm2)
symbpercm2.vs.resp.lm
# t = 2.1064, df = 37, p-value = 0.04201, cor = 0.3272279
# symbdens vs. lipids/cm2 (note: Loess function/GAM with method = REML, not linear):
symbpercm2.vs.lip <- ggplot(data, aes(x = symb_cm2 / 1000, y = lipids.mg.cm2)) +
  geom_point(size = 2, aes(color = hist_symb, 
                           shape = treatment)) +
  ylab("Lipids (mg/cm2)") + 
  xlab("Symbiont density (cells µm-2)") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  scale_shape_manual(values = c(2, 1, 17, 16))+
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14)) +
  labs(title = "") +
  labs(color = "Natural history", shape = "treatment") +
  geom_smooth(color = "black", method = loess, size=0.7)
symbpercm2.vs.lip
symbpercm2.vs.lip.gam <- gam(lipids.mg.cm2 ~ s(symb_cm2), method = "REML", data = data)
symbpercm2.vs.lip.gam
anova(symbpercm2.vs.lip.gam) # no effect; edf = 2.571, F = 1.193, p = 0.292

# Symbiont density supplemental figure (Fig. S9)
SymdensSupp <- ggarrange(symbpercm2.vs.lip, symbpercm2.vs.calc, symbpercm2.vs.chlsymb, symbpercm2.vs.resp,
                         nrow = 2, ncol = 2, common.legend = T)
SymdensSupp

# Gross photosynthesis vs. host carbon assimilation (Fig S10)
GPprot.vs.host13C <- ggplot(data, aes(x = gp_ugprot * 60, y = host_13C)) +
  geom_point(size = 2, aes(color = hist_symb, 
                           shape = treatment)) +
  ylab("Host C assimilation (13C at-%)") + 
  xlab("Photosynthesis (umol O2/ug host protein/hr)") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  scale_shape_manual(values = c(2, 1, 17, 16))+
  theme_bw()  + 
  theme(panel.grid = element_blank(), text = element_text(size = 14)) +
  labs(title = "") +
  labs(color = "Natural history", shape = "treatment") +
  geom_smooth(color = "black", method=lm, size=0.7)
GPprot.vs.host13C
GPprot.vs.host13C.lm <- cor.test(data$gp_ugprot, data$host_13C)
GPprot.vs.host13C.lm
# t = 0.33062, df = 38, p-value = 0.7427, R2 = 0.00287

