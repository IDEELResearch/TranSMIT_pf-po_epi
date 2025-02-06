################################################
# Plot P. ovale and falciparum interaction by season
# Author: Kelly Carey-Ewend
# Date: 15 November 2024
# Last modified: January 2025
################################################


#set working directory with input (and output) files
setwd("~/Documents/P. ovale Epi Bagamoyo")

#load necessary R packages (some may not be needed). Install new packages with install.packages("newpackage")
suppressPackageStartupMessages({
  library(survey)
  library(srvyr)
  library(tableone)
  library(EpiStats)
  library(lmtest)
  library(gplots)
  library(ggplot2)
  library(ggtext)
  library(tidyverse)
  library(fastDummies)
  library(gmodels)
  library(readxl)
  library(mgcv)
  library(splines)
  library(plotrix)
  library(nloptr)
  library(lme4)
  library(flextable)
  library(RColorBrewer)
})

#Determine colorblind-friendly color palettes
library("RColorBrewer")
display.brewer.all(colorblindFriendly = TRUE)
display.brewer.pal(n=4,"Dark2")
brewer.pal(n=4,"Dark2")


##import csv file containing season, ID number (arbitrary for ordering of groups), log(estimate), log(lower CI), log(upper CI) and p-value
data <- read.csv("TranSMIT_PfPo_PR.csv", header = FALSE)

#assign names to variables in table
#ID is given to each observation and is used for ordering and assigning colors in plot; is arbitrary
colnames(data) <- c("season","ID","point","lci","uci","p-value")
for (n in 1:length(data)){
  if (data$season[n] == "Long wet season"){
    data$season.name[n] <- "*Masika* rains"
  } else if (data$season[n] == "Short wet season"){
    data$season.name[n] <- "*Vuli* rains"
  } else if (data$season[n] == "Dry season") {
    data$season.name[n] <- "Dry season"
  } else {
    data$season.name[n] <- "Overall"
  }
}

#generate text expression with italicized portions for use in graphing
#first phrase places arrows showing the significances of prevalence ratios on either side of the null
y_title <- "Fewer co-infections \u2194  More co-infections <br><br>*Po* PR (95% CI) between *Pf*+ and *Pf*-"

#generate graph
pfpo_plot <- ggplot(data, aes(x=ID, y=point)) +
  #add null value line
  geom_hline(yintercept=1, linetype='dashed') +
  #plot confidence intervals
  geom_pointrange(aes(x=ID, y=point, ymin=lci, ymax=uci), shape=15, size=0.8, color="black", show.legend=F, fatten=0.2) + 
  #plot point estimate
  geom_point(shape=15, size=3, color="black", show.legend=F, alpha=0.7) + 
  #plot point estimate with color (if desired)
  #geom_point(shape=15, size=5, aes(color=as.factor(ID)), show.legend=F, alpha=0.7) + 
  #assign color to specific observations, if desired
  #scale_color_manual(values = c("1" = "darkblue", "2" = "#D95F02", "3" = 'grey50', "4" = '#E6AB02')) +
  #flip coordinates if desired
  coord_flip() + 
  #change theme to black-white
  theme_bw() + 
  #determine x-axis breaks
  scale_x_continuous(breaks=data$ID, labels=data$season.name, trans = "reverse") + 
  #determine y axis breaks
  scale_y_log10(breaks = c(0.13, 0.25, 0.5, 1, 2, 4, 8), lim = c(0.12,8.1)) +
  #adjust labels, font sizes
  labs(x="Season", y=y_title) + 
  theme(axis.text.y = element_markdown(size = 16), axis.text.x = element_markdown(size=16, angle = 0), 
        axis.title.y = element_markdown(size = 18), axis.title.x = element_markdown(size = 18),
        #axis.ticks.y=element_blank(),
        panel.grid.minor=element_blank()) 

#print plot
pfpo_plot
#save plot to png file in figures directory (within working directory)
ggsave("Figures/pfpo_antagonism.png", plot = pfpo_plot, device = "png", dpi = 600, width = 6.5, height = 5)



#save version of pfpo_plot for stacking with the poc and pow-specific associations
y_title <- expression(paste(italic("Po")," PR (95% CI) between ", italic("Pf"),"+ and ",italic("Pf")))

stacked_pfpo_plot <- ggplot(data, aes(x=ID, y=point)) +
  #add null value line
  geom_hline(yintercept=1, linetype='dashed') +
  #plot confidence intervals
  geom_pointrange(aes(x=ID, y=point, ymin=lci, ymax=uci), shape=15, size=0.8, color="black", show.legend=F, fatten=0.2) + 
  #plot point estimate
  geom_point(shape=15, size=3, color="black", show.legend=F, alpha=0.7) + 
  #plot point estimate with color (if desired)
  #geom_point(shape=15, size=5, aes(color=as.factor(ID)), show.legend=F, alpha=0.7) + 
  #assign color to specific observations, if desired
  #scale_color_manual(values = c("1" = "darkblue", "2" = "#D95F02", "3" = 'grey50', "4" = '#E6AB02')) +
  #flip coordinates if desired
  coord_flip() + 
  #change theme to black-white
  theme_bw() + 
  #determine x-axis breaks
  scale_x_continuous(breaks=data$ID, labels=data$season.name, trans = "reverse") + 
  #determine y axis breaks
  scale_y_log10(breaks = c(0.13, 0.25, 0.5, 1, 2, 4), lim = c(0.12,4.2)) +
  #adjust labels, font sizes
  labs(x="Season", y=y_title) + 
  theme(axis.text.y = element_text(size = 11), axis.text.x = element_text(size=12), 
        axis.title.y = element_text(size = 11), axis.title.x = element_text(size = 12),
        #axis.ticks.y=element_blank(),
        panel.grid.minor=element_blank()) 
ggsave("Figures/pfpo_antagonism.png", plot = stacked_pfpo_plot, device = "png", dpi = 600, width = 6.5, height = 3.4)


#Use following code to plot additional graphs of interaction of Pf to specific Po species
#repeat above for additional dataset showing species-specific interactions (Poc and Pow)
pf_poc_pow <- read.csv("TranSMIT_Pf-vs-Poc-Pow_PR.csv", header = FALSE)
colnames(pf_poc_pow) <- c("species","ID","point","lci","uci","p-value")

#adjust titles and labels for italicization
y_title <- expression(paste("Adjusted ",italic("Po")," species prevalence ratio (95% CI) between ", italic("Pf"),"+ and ",italic("Pf"),"-"))
x_title <- expression(paste(italic("Po")," species"))
pf_poc_pow$label <- c(expression(paste(italic("Poc"))),expression(paste(italic("Pow"))))

#create plot of species-specific interactions
pf_poc_pow_plot <- ggplot(pf_poc_pow, aes(x=as.factor(ID), y=point)) +
  geom_hline(yintercept=1, linetype='dashed') +
  geom_pointrange(aes(x=as.factor(ID), y=point, ymin=lci, ymax=uci), shape=15, size=0.8, color="black", show.legend=F, fatten=0.2) + 
  #plot points without color
  geom_point(shape=15, size=3, color="black", show.legend=F, alpha=0.7) + 
  #plot points with color, if desired
  #geom_point(shape=15, size=3, aes(color=as.factor(ID)), show.legend=F, alpha=0.7) + 
  coord_flip() + theme_bw() + 
  #scale_color_manual(values = c("1" = "#1B9E77", "2" = "#7570B3")) +
  scale_x_discrete(label=c(expression(paste(italic("Poc"))),expression(paste(italic("Pow"))))) + 
  scale_y_log10(breaks = c(0.13, 0.25, 0.5, 1, 2, 4), lim = c(0.12,4.2)) +
  labs(x=x_title, y=y_title) + 
  theme(axis.text.y = element_text(size = 11), axis.text.x = element_text(size=12), 
        axis.title.y = element_text(size = 11), axis.title.x = element_text(size = 12),
        #axis.ticks.y=element_blank(),
        panel.grid.minor=element_blank()) 

#print plot
pf_poc_pow_plot
#save plot
ggsave("Figures/pf-poc-pow_antagonism.png", plot = pf_poc_pow_plot, device = "png", dpi = 600, width = 6.5, height = 2.5)

#load ggpubr package for organizing plots in same image
library(ggpubr)

#arrange two plots vertically, with extra space so the axes line up
#first, create empty plot
blank <- ggplot() + theme_void()

#create lower row of blank plot, then pf-poc/pow plot. adjust widths manually so y-axes match up
bottom_row <- ggarrange(blank, pf_poc_pow_plot, ncol = 2, widths = c(1.07,6.93))
#create plot of pf-po plot on top, pf-poc/pow plot on bottom
dual_plot <- ggarrange(pfpo_plot, bottom_row, ncol=1, nrow=2, widths = c(8,8), heights = c(3.4,2))
#print plot
dual_plot

#save plot
ggsave("Figures/pf-po-poc-pow_antagonism.png", plot = dual_plot, device = "png", dpi = 600, width = 6.5, height = 5.5)

