# Plotting PI models for Allen-Waller & Barott 2023 supplement
# Luella Allen-Waller, as adapted from Teegan Innis (Innis et al. 2021 GCB)

# load libraries
library(dplyr)
library(ggplot2)
library(csv)
library(tidyverse)
library(ggpubr)
library(readr)
library(Rmisc)
library(rstatix)
library(gridExtra)
library(data.table)


####--------------------------------------------------------
#### Plot: 4 panels: column 1 = june, column 2 = july
# load raw data and ensure frag ID is a factor
acute <- read.csv("2019 Mcap Acute Respirometry Raw.csv")
acute$ID <- as.factor(acute$ID)
init.recv <- read.csv("2019 Mcap Initial Recov Respirometry Raw.csv.csv")
init.recv$ID <- as.factor(init.recv$ID)

# normalize to surface area, account for chamber volume, and multiply by 60 so it's an hourly rate instead of by minute
SA <- read.csv("2019Mcap_surface_areas.csv")
SA$ID <- as.factor(SA$ID)
acute <- merge(x = acute, y = SA, by = "ID", all.x = T)
acute$umolcm2hr <- (acute$coefficients.x1 * 60 * acute$chamber_vol) / acute$SA
init.recv <- merge(x= init.recv, y = SA, by = "ID", all.x = T)
init.recv$umolcm2hr <- (init.recv$coefficients.x1 * 60 * init.recv$chamber_vol) / init.recv$SA

# add colony factors
key <- read.csv("2019 Mcap frag key.csv")
key$ID <- as.factor(key$ID)
acute <- merge(x = acute, y = key, by = "ID", all.x = T)
init.recv <- merge(x = init.recv, y = key, by = "ID", all.x = T)

# combine into one df
PI <- bind_rows(acute, init.recv)
# eliminate NBC
PI$histsymb <- as.factor(paste(PI$hist, PI$symb))
PI <- subset(PI, histsymb != "NB C")

# split by temp
PIHot <- subset(PI, temp == "31º")
PIAmb <- subset(PI, temp == "28º")

# Plot heat-treated
PICurvesHot <- ggplot(PIHot, aes(x = rad, y = umolcm2hr, color = symb, group = ID)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 2), se = F, size = 0.5) +
  scale_x_continuous(breaks = c(0, 200, 400, 600, 800)) +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  scale_fill_manual(values = c("#92C5DE", "#D6604D")) +
  facet_grid(. ~ timepoint)  + coord_cartesian(ylim = c(-5,8)) +
  labs(y = expression(atop("Oxygen Evolution", paste((µmol~O["2"]~cm^-2~hr^-1)))), x = expression(atop("Irradiance", paste((umol~m^-2~s^-1)), color = "Symbiont"))) +
  theme(aspect.ratio = 1, 
        axis.title = element_text(size = 14), axis.text = element_text(size = 14, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.key = element_rect(fill = NA),
        strip.background = element_blank(), strip.text = element_blank()) 
PICurvesHot

# Plot ambient
PICurvesAmbient <- ggplot(PIAmb, aes(x = rad, y = umolcm2hr, color = symb, group = ID)) +
  stat_smooth(method = "lm", formula = y ~ poly(x, 2), se = F, size = 0.5) +
  scale_x_continuous(breaks = c(0, 200, 400, 600, 800)) +
  scale_color_manual(values = c("#92C5DE", "#D6604D")) +
  scale_fill_manual(values = c("#92C5DE", "#D6604D")) +
  facet_grid(. ~ timepoint) + coord_cartesian(ylim = c(-5,8)) +
  labs(y = expression(atop("Oxygen Evolution", paste((µmol~O["2"]~cm^-2~hr^-1)))), x = expression(atop("Irradiance", paste((umol~m^-2~s^-1)), color = "Symbiont"))) +
  theme(aspect.ratio = 1, 
        axis.title = element_text(size = 14), axis.text = element_text(size = 14, color = "black"),
        panel.background = element_blank(), panel.border = element_rect(fill = NA),
        legend.position = "right", legend.key = element_rect(fill = NA),
        strip.background = element_blank(), strip.text = element_blank()) 
PICurvesAmbient

# Plot both together
PICurves <- ggarrange(PICurvesHot, PICurvesAmbient,
                      ncol = 1, nrow = 2, common.legend = T)
PICurves
