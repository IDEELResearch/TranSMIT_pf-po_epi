################################################
# Plot Pf and Po prevalence by villages
# Author: Kelly Carey-Ewend 
# Date 5 November 2024
# Modified: 15 January 2025
################################################

#load needed packages
suppressPackageStartupMessages({
  library("ggplot2")
  library("reshape2")
})

#set working directory
setwd("~/Documents/P. ovale Epi Bagamoyo")

#load in manually-exported table of Pf and Po prevalences by village
prevs <- read.table("TranSMIT Epi Data Tables - Village Prevalence.tsv", header = TRUE, sep = "\t")

#assign column names
names(prevs) <- c("Village", "Po","Pf")

#test correlation (Pearson's) of village-wise Pf and Po prevalence
cor.test(prevs$Po,prevs$Pf)

#sort data by ascending Pf prevalence
prevs.sort <- prevs[order(prevs$Pf),]

#create long dataset, with a new observation for each Pf and Po prevalence
prevs.long <- melt(prevs.sort)

#create sorted village factor variable based on villages present in earlier dataset
prevs.long$sort.village <- factor(prevs.long$Village, levels = prevs.sort$Village)

#create barplot
village_prev <- ggplot(prevs.long, aes(x = sort.village,value,fill=variable)) + 
  geom_bar(position = "dodge", stat = "identity") +
  theme(plot.title = element_text(hjust = 0.5, size = 20),
        legend.position = "right",
        legend.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_text(size=18),
        axis.title.x = element_text(size=18),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(angle = 60, size = 16, vjust = 0.5)) +
  labs(x = "Village", 
       y = "Prevalence (%)")

#print plot
village_prev

#save plot
ggsave('Figures/transmit_pf-po_village-prev.png', village_prev, height = 4, width = 7, dpi = 600)
