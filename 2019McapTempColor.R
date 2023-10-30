# Mcap 2019 temperature and color graphing and analysis
# Script for Allen-Waller & Barott 'Symbiotic dinoflagellates divert energy away from mutualism during coral bleaching recovery' Symbiosis 2023
# Figure 1
# Luella Allen-Waller

# load packages
library(dplyr)
library(ggplot2)
library(csv)
library(tidyverse)
library(ggpubr)
library(readr)
library(Rmisc)
library(lme4)
library(lmerTest)
library(rstatix)
library(patchwork)
library(ggpmisc)
library(lubridate)
library(MuMIn)
library(scales)

# load data
color <- read.csv("Mcap 2019 FvFm and Color.csv")
# fix date parsing
color$date <- color$day %>% mdy()
# tank temperatures and eDHW:
eDHW <- read.csv("MpHi temp eDHW.csv")
eDHW$Date <- eDHW$Date %>% mdy()

# Color over time
colortime <- ggplot(color, aes(x = day, y = avg_red_intensity, color = hist_symb, linetype = treatment)) +
  stat_summary(fun = mean, geom = "path", size = 0.7, aes(group = interaction(treatment,hist_symb))) +
  stat_summary(fun.data = mean_se, size = 0.7, geom = "errorbar", aes(width=1, group=interaction(treatment,hist_symb))) +
  ylab(expression(paste("Color Score")))+
  xlab("Date") +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  scale_linetype_manual(values = c('solid','dotted'))+
  scale_x_date(limit=c(as.Date("2019-06-26"),as.Date("2019-07-22")), 
               breaks = c(as.Date("2019-06-27"),as.Date("2019-07-08"),as.Date("2019-07-21")),
               date_labels = "%b %d")+
  theme_bw()  + 
  theme(axis.title=element_text(size=12, color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_rect(fill = NA), legend.title = element_text(size = 10, color = "black"), legend.key = element_blank(), 
        legend.text = element_text(size = 10, color = "black"), legend.position = "none", legend.key.size = unit(0.3, "cm")) +
  labs(title = "") +
  labs(color = "Bleaching Susceptibility", linetype = "Treatment")+
  coord_cartesian(ylim = c(0.85,1.15))
colortime
# lmem: did color change with time in response to treatment?
# various models with temp, bleach hist, plus initial tank as random
color$factor_date <- as.factor(color$day)
color.lme <- lmer(avg_red_intensity ~ factor_date*temp*bleach_hist + (1|initial_tank), data = color)
color.lm <- lm(avg_red_intensity ~ factor_date*temp*bleach_hist, data = color)
color.lme2 <- lmer(avg_red_intensity ~ factor_date*temp + (1|initial_tank), data = color)
color.lm2 <-  lm(avg_red_intensity ~ factor_date*temp, data = color)
AICc(color.lme,color.lm,color.lme2,color.lm2)
# simplest is best - color.lm2
plot(color.lm2)
# no points with excessive leverage, and pretty normal
summary(color.lme)
# significant effect of date *** and temp*** and interaction**
# tukey on full LM
color.tukey <- emmeans(color.lm, list(pairwise ~ temp*factor_date*bleach_hist), simple = "temp", adjust = "tukey")
color.tukey
### differences between dates
# 31 B: first TP ***different from the rest
# 28 NB: first TP *different from final
# 31 NB: first TP **different from second and *different from third
# second and third timepoints do not differ from each other for any groups
# no timepoints differ for 28 B
### differences between temps:
# color differs by temp for every bleaching history at every date ()
### differences b/w bleach histories: NONE
# so now a version with no bleach hist:
color.tukey2 <- emmeans(color.lm2, list(pairwise ~ temp*factor_date), simple = "temp", adjust = "tukey")
color.tukey2
# ***at each date for temp
# for 28, just T1 and final
# for 31, both July different from T1

# 2-way repeated measures ANOVA to make use of repeated genet sampling
# generate summary stats for each temp-genet category at each day
colorsummary <-color %>%
  group_by(temp, parent, day, bleach_hist) %>%
  get_summary_stats(avg_red_intensity, type = "mean_sd")
data.frame(colorsummary)
# pairwise:
pair<-color %>% 
  pairwise_t_test(avg_red_intensity~date,paired=TRUE, p.adjust.method = "bonferroni" ) 
data.frame(pair) 
# based on temperature alone:
pwc <- color %>%
  group_by(temp, parent, bleach_hist) %>%
  t.test(
    avg_red_intensity ~ day, paired = T,
    p.adjust.method = "bonferroni"
  )
data.frame(pwc)

# Temperatures: mean with min/max ribbon
MpHiTreatmentTemps <- ggplot(eDHW, aes(x = Date, color = treatment)) + 
  geom_line(aes(y = Mean), size = 0.5) +
  geom_ribbon(aes(ymin = Min, ymax = Max, fill = treatment), alpha = 0.3, linetype = "blank") +  
  scale_color_manual(values = c("heat" = "#990000", "ambient" = "#FF9999")) +
  theme(aspect.ratio = .3, axis.title=element_text(size=12, color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_rect(fill = NA), legend.title = element_text(size = 10, color = "black"), legend.key = element_blank(), 
        legend.text = element_text(size = 10, color = "black"), legend.position = "none", legend.key.size = unit(0.3, "cm")) +
  scale_fill_manual(values=alpha(c("#FF9999", "#990000")), name = "") + 
  labs(y = "temperature (°C)", x = "", color = "", title = "") +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  geom_vline(xintercept = 2019-06-27, color = "black", size = 0.5, linetype = "dotted") +
  geom_vline(xintercept = 2019-02-21, color = "black", size = 0.5, linetype = "dotted")
MpHiTreatmentTemps
# eDHW
# Make one df for each treatment
HT.eDHW <- subset(eDHW, treatment == "heat")
ambient.eDHW <- subset(eDHW, treatment == "ambient")
# bind everything together
eDHW.bytreatment <- rbind(HT.eDHW,ambient.eDHW)
# Graph eDHW by treatment
MpHi.eDHW <- ggplot(eDHW.bytreatment, aes(x = Date, color = treatment)) + 
  geom_line(aes(y = cumul_eDHW_meanover29), size = 0.5) + 
  scale_color_manual(values = c("heat" = "#990000", "ambient" = "#FF9999")) +
  theme(aspect.ratio = .3, axis.title=element_text(size=12, color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_rect(fill = NA), legend.title = element_text(size = 10, color = "black"), legend.key = element_blank(), 
        legend.text = element_text(size = 10, color = "black"), legend.background = element_blank(), legend.position = c(0.8, 0.55), legend.key.size = unit(0.3, "cm")) +
  scale_fill_manual(values=alpha(c("#FF9999", "#990000")), name = "") + 
  labs(y = "eDHW (°C)", x = "Date", color = "", title = "") +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  geom_vline(xintercept = 14, color = "black", size = 0.4, linetype = "dotted") +
  geom_vline(xintercept = 46.5, color = "black", size = 0.4, linetype = "longdash")
MpHi.eDHW

# stack temp and eDHW
MpHiTempsStack <- ggarrange(MpHiTreatmentTemps, MpHi.eDHW.Days,
  ncol = 1, nrow = 2)
MpHiTempsStack

# Higher-res temp data
tanktemps <- read.csv('DetailMpHiTankTemps.csv')
tanktemps$Date <- tanktemps$Date %>% mdy()
tanktemps$date_time <-  ymd_hms(paste(tanktemps$Date,tanktemps$Time))
# Clean so that I have one variable per column
tanktemps <-  tanktemps%>% gather('Tmp.9','Tmp.10','Tmp.11','Tmp.12','Tmp.13','Tmp.14',key='tank',value='temp')
# Add a column to raw data output denoting treatment group with ambient as default
tanktemps <-  mutate(tanktemps,treatment = 'ambient')
# assign HT tanks
tanktemps$treatment <- with(tanktemps, ifelse(tank=='Tmp.10', 'heat', treatment))
tanktemps$treatment <- with(tanktemps, ifelse(tank=='Tmp.12', 'heat', treatment))
tanktemps$treatment <- with(tanktemps, ifelse(tank=='Tmp.14', 'heat', treatment))
# Check date format, fix, and arrange according to tank and then date within tank
str(tanktemps)
tanktemps <-  tanktemps %>% arrange(tanktemps$tank, tanktemps$Date)

# Plot each tank over the whole experiment
# use black and white aesthetics
# add rectangle for the treatment period
MpHiTankTempsDetailed <- ggplot(tanktemps,aes(x=date_time)) +
  geom_rect(aes(xmin=as.POSIXct("2019-06-14"), xmax=as.POSIXct("2019-06-27"),ymin=26,ymax=33), color="#ffffff", fill = "#ffcccc", alpha=0.05)+
  geom_line(aes(y=temp, linetype=treatment, group = tank), size = 0.7, alpha = 0.4) +
  scale_linetype_manual(values = c("heat" = "dotted","ambient" = "solid")) +
  theme(axis.title=element_text(size=12, color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_rect(fill = NA), legend.title = element_text(size = 10, color = "black"), legend.key = element_blank(), 
        legend.text = element_text(size = 10, color = "black"), legend.background = element_blank(), legend.position = c(0.8, 0.55), legend.key.size = unit(0.3, "cm")) +
  labs(y="Temperature (ºC)", x="Date",linetype="Treatment",title="")+
  scale_x_datetime(breaks = c(as.POSIXct("2019-06-14"),as.POSIXct("2019-06-27"),as.POSIXct("2019-07-08"),as.POSIXct("2019-07-21"), as.POSIXct("2019-07-31")),
                   labels = date_format(format = "%b %d")) +
  coord_cartesian(ylim = c(26.9,32.2))
MpHiTankTempsDetailed
# remake eDHW with BW aesthetics
BW.MpHi.eDHW.Days <- ggplot(eDHW.bytreatment, aes(x = Date, linetype = treatment)) + 
  geom_line(aes(y = cumul_eDHW_meanover29), size = 0.7) + 
  scale_linetype_manual(values = c("heat" = "dotted","ambient" = "solid")) +
  theme(axis.title=element_text(size=12, color = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
        panel.border = element_rect(fill = NA), legend.title = element_text(size = 10, color = "black"), legend.key = element_blank(), 
        legend.text = element_text(size = 10, color = "black"), legend.background = element_blank(), legend.position = c(0.8, 0.55), legend.key.size = unit(0.3, "cm")) +
  labs(y="eDHW (ºC)", x="Date",linetype="Treatment",title="")
BW.MpHi.eDHW.Days

# stack temp and eDHW
MpHiBWTempsStack <- ggarrange(MpHiTankTempsDetailed, BW.MpHi.eDHW.Days,
                            ncol = 1, nrow = 2)
MpHiBWTempsStack

# add color, with room for frag pics
blankbox <-  ggplot() + theme(axis.title=element_text(size=12, color = "black"),
                              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                              panel.border = element_blank(), legend.title = element_text(size = 10, color = "black"), legend.key = element_blank(), 
                              legend.text = element_text(size = 10, color = "black"), legend.background = element_blank(), legend.position = c(0.8, 0.55), legend.key.size = unit(0.3, "cm"))
blankbox
MpHiColorTempSquare <- ggarrange(MpHiTankTempsDetailed, blankbox, BW.MpHi.eDHW.Days, colortime,
                              ncol = 2, nrow = 2, 
                              widths = c(1,0.7), heights = c(1,0.55))
MpHiColorTempSquare

# Make a version all on one row without eDHW
MpHiRow <-ggarrange(MpHiTankTempsDetailed, blankbox, colortime,blankbox,
                      ncol = 4, nrow = 1, 
                      widths = c(1.2,1,0.7,0.1))
MpHiRow
