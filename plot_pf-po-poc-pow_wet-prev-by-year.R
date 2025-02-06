################################################
# Plot Plasmodium prevalences within either wet season by year
# Author: Kelly Carey-Ewend 
# Date 10 November 2024
# Modified: 12 January 2025
################################################

#load necessary packages
suppressPackageStartupMessages({
  library(dplyr)
  library(plyr)
  library(tidyverse)
  library(ggtext)
  library(forcats)
  library(ggpubr)
})

#set working directory
setwd("~/Documents/P. ovale Epi Bagamoyo")

#load in table of data values
year_prev <- read.csv("TranSMIT_prev_by_year_inwetseason.csv")

#flip to long-format data
year_prev_long <- year_prev %>% pivot_longer(!Year.of.Study,names_to = "species", values_to = "prevalence")

#divide dataset into Pf and Po prevalences
year_prev_long_pfpo <- year_prev_long[(year_prev_long$species == 'Pf' | year_prev_long$species == 'Po'),]

#relevel species in Pf, Po order
year_prev_long_pfpo$species <-fct_relevel(year_prev_long_pfpo$species, 'Pf','Po')

#assign new year name including sample size on next row for plotting
year_prev_long_pfpo$Year.of.Study <- revalue(year_prev_long_pfpo$Year.of.Study, c(
  "2018 (n=159)" = "2018\n(n=159)", 
  "2019 (n=1637)" = "2019\n(n=1637)", 
  "2020 (n=326)" = "2020\n(n=326)", 
  "2021 (n=1088)" = "2021\n(n=1088)",
  "2022 (n=766)" = "2022\n(n=766)"
))

#create plot
prevalence_by_year_pfpo_plot <- ggplot(year_prev_long_pfpo, aes(fill=species, y=prevalence, x=Year.of.Study)) + 
  geom_bar(position="dodge", stat="identity",show.legend = T) +
  #manually assign colors and italicization to each species; Poc and Pow not used in this plot
  scale_fill_manual (values = c('Pf'= 'darkblue','Po'='#D95F02',"Poc"='#1B9E77', 'Pow' = '#7570B3'),
                     labels = c('Pf'= expression(paste(italic("Pf"))),
                                'Po'= expression(paste(italic("Po"))),
                                'Poc'= expression(paste(italic("Poc"))),
                                'Pow'= expression(paste(italic("Pow"))))) +
  labs(y="Wet Season Prevalence", x="Year of Study") +
  theme_bw() +
  theme(legend.title = element_blank())

#print plot
prevalence_by_year_pfpo_plot

#divide larger dataset into Poc and Pow prevalences
year_prev_long_pocpow <- year_prev_long[(year_prev_long$species == 'Poc' | year_prev_long$species == 'Pow'),]

#assign correct sample size (now excluding Po-positive samples in which species-identification assay was not run) from precalence denominator, 
#and put sample size on new row
year_prev_long_pocpow$Year.of.Study <- revalue(year_prev_long_pocpow$Year.of.Study, c(
  "2018 (n=159)" = "2018\n(n=158)", 
  "2019 (n=1637)" = "2019\n(n=1557)", 
  "2020 (n=326)" = "2020\n(n=319)", 
  "2021 (n=1088)" = "2021\n(n=1046)",
  "2022 (n=766)" = "2022\n(n=760)"
))

#create plot
prevalence_by_year_pocpow_plot <- ggplot(year_prev_long_pocpow, aes(fill=species, y=prevalence, x=Year.of.Study)) + 
  geom_bar(position="dodge", stat="identity") +
  #manually assign colors and italicization to each species; Po and Pf not used in this plot
  scale_fill_manual (values = c('Pf'= 'darkblue','Po'='#D95F02',"Poc"='#1B9E77', 'Pow' = '#7570B3'),
                     labels = c('Pf'= expression(paste(italic("Pf"))),
                                'Po'= expression(paste(italic("Po"))),
                                'Poc'= expression(paste(italic("Poc"))),
                                'Pow'= expression(paste(italic("Pow"))))) +
  labs(y="Wet Season Prevalence", x="Year of Study") +
  theme_bw() +
  theme(legend.title = element_blank())

#print plot
prevalence_by_year_pocpow_plot

#arrange pf-po and poc-pow plots on top of each other
prev_by_year_plot <- ggarrange(prevalence_by_year_pfpo_plot, prevalence_by_year_pocpow_plot, nrow = 2)

#print dual plot
prev_by_year_plot

#save plot
ggsave('Figures/prev_by_year.png', prev_by_year_plot, height = 5, width = 6, dpi = 600)
