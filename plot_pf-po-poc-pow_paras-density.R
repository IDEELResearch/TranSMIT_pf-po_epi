################################################
# Plot Pf-Po and Poc-Pow Parasite Density Histograms
# Author: Kelly Carey-Ewend 
# Date 25 November 2024
# Modified: 12 January 2025
################################################

#load necessary packages
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
})


#########generate Pf and Po parasite density / Ct histograms
#set working directory
setwd("~/Documents/P. ovale Epi Bagamoyo")

#load in data
para_dens <- read.csv("transmit_parasite_density.csv")

#assign infection status to all positive samples based on positivity for either Po or Pf
para_dens$status <- NA
for (i in 1:nrow(para_dens)) {
  if (para_dens$pf[i] == 1 & para_dens$po[i] == 1) {
    para_dens$status[i] <- "Pf-Po mixed infection"
  }else if (para_dens$pf[i] == 1 & para_dens$po[i] == 0) {
    para_dens$status[i] <- "Pf mono-infection"
  }else if (para_dens$pf[i] == 0 & para_dens$po[i] == 1) {
    para_dens$status[i] <- "Po mono-infection"
  }}

#reformat status so it appears in common legend
para_dens$status <- factor(para_dens$status, levels = c("Pf mono-infection","Pf-Po mixed infection","Po mono-infection"))

#subset to Pf-positives and Po-positives
pf_para_dens <- para_dens[para_dens$pf == 1,]
po_para_dens <- para_dens[para_dens$po == 1,]

#plot histogram of Pf parasite density
options(scipen = 999)
pf_density_plot <- ggplot(pf_para_dens, aes(x=qpcr_pfdens_screening, fill=status)) + 
  #geom_histogram(aes(y=after_stat(density)), show.legend = T) + 
  geom_density(aes(y=after_stat(density)), show.legend = T, alpha = 0.5) +
  #add line showing submicroscopic threshold (<100p/uL)
  geom_vline(xintercept = 100, linetype = "dashed", color = 'grey50') +
  theme_bw() + 
  scale_fill_manual(values = c("Pf mono-infection" = "skyblue2","Pf-Po mixed infection" = 'grey50',"Po mono-infection" = '#D95F02'), drop = FALSE)+
  labs(x="Pf parasite density per microliter", y = "Proportion",fill=element_blank()) + 
  theme(plot.title = element_text(hjust =0.5)) + 
  scale_x_log10(breaks=c(0.1,10,1000,100000)) +
  theme_bw()+ 
  theme(axis.title.y = element_text(size = 15),axis.title.x = element_text(size=15),title = element_text(size = 20), legend.text = element_text(size = 15),axis.text= element_text(size = 12))

#print plot
pf_density_plot

#plot histogram of Po Ct values
po_density_plot <- ggplot(po_para_dens, aes(x=po18s_ct, fill=status)) + 
  # geom_histogram(aes(y=after_stat(density)),show.legend = T) + 
  geom_density(aes(y=after_stat(density)), show.legend = T, alpha = 0.5) +
  #add line showing submicroscopic threshold (<100p/uL), corresponding to Po18S Ct of 37.0598 (as calculated from a representative standard run, though not run in each Po18S assay)
  geom_vline(xintercept = 37.0598, linetype = "dashed", color = 'grey50') +
  theme_bw() + 
  scale_fill_manual(values = c("Pf mono-infection" = "skyblue2","Pf-Po mixed infection" = 'grey50',"Po mono-infection" = '#D95F02'), drop = FALSE)+
  labs(x="po18S cyclic threshold", y = "Proportion",fill=element_blank()) + 
  theme(plot.title = element_text(hjust =0.5)) + 
  scale_x_reverse() +
  theme_bw()+ 
  theme(axis.title.y = element_text(size = 15),axis.title.x = element_text(size=15),title = element_text(size = 20), legend.text = element_text(size = 15),axis.text= element_text(size = 12))

#print plot
po_density_plot

#separate densities between mono and mixed infections of Pf-positives
pf_mono_para_dens <- pf_para_dens$qpcr_pfdens_screening[pf_para_dens$status == "Pf mono-infection"]
pf_mixed_para_dens <- pf_para_dens$qpcr_pfdens_screening[pf_para_dens$status == "Pf-Po mixed infection"]
#test for difference in Pf density distributions between Pf and Pf-Po infections
ks.test(pf_mono_para_dens,pf_mixed_para_dens, alternative = c("two.sided"),exact = NULL, simulate.p.value = FALSE, B = 2000)

#separate densities between mono and mixed infections of Po-positives
po_mono_para_dens <- po_para_dens$po18s_ct[po_para_dens$status == "Po mono-infection"]
po_mixed_para_dens <- po_para_dens$po18s_ct[po_para_dens$status == "Pf-Po mixed infection"]
#test for difference in Pf density distributions between Pf and Pf-Po infections
ks.test(po_mono_para_dens,po_mixed_para_dens, alternative = c("two.sided"),exact = NULL, simulate.p.value = FALSE, B = 2000)

#arrange Pf and Po density plots side-by-side
dual_density_plot <- ggarrange(pf_density_plot,po_density_plot,ncol =2,common.legend = TRUE,legend="right")
dual_density_plot

#save Pf-Po density plot
ggsave('Figures/pf-po_parasite_density_plots.png', dual_density_plot, height = 6, width = 12, dpi = 600)


###generate plots of Poc and Pow parasite density

#set working directory
setwd("~/Documents/P. ovale Epi Bagamoyo")

#load in Poc-Pow parasite densities
po_dens <- read.csv("transmit_parasite_density_final.csv")

#assign category for subsetting and labeling
for (i in 1:nrow(po_dens)) {
  if (po_dens$mixed[i] == 1) {
    po_dens$status[i] <- "Poc-Pow mixed infection"
  }else if (po_dens$poc[i] == 1) {
    po_dens$status[i] <- "Poc mono-infection"
  }else if (po_dens$pow[i] == 1) {
    po_dens$status[i] <- "Pow mono-infection"
  }}
po_dens$status <- factor(po_dens$status, levels = c("Poc-Pow mixed infection","Poc mono-infection","Pow mono-infection"))

#subset by species
poc_dens <- po_dens[po_dens$poc==1,]
pow_dens <- po_dens[po_dens$pow==1,]

#remove scientific notation from axis labels
options(scipen = 999)

#plot histogram of Poc parasite density
poc_density_plot <- ggplot(poc_dens, aes(x=poc_parasite_density, fill=status)) + 
  #geom_histogram(aes(y=after_stat(density)), show.legend = T) + 
  geom_density(aes(y=after_stat(density)), show.legend = T, alpha = 0.5) +
  theme_bw() + 
  scale_fill_manual(values = c("Poc-Pow mixed infection" = "#E6AB02","Poc mono-infection" = '#1B9E77', "Pow mono-infection" = '#7570B3'), drop = FALSE)+
  labs(x="Poc parasite density per microliter", y = "Proportion",fill=element_blank()) + 
  theme(plot.title = element_text(hjust =0.5)) + 
  scale_x_log10() +
  theme_bw()+ 
  theme(axis.title.y = element_text(size = 15),axis.title.x = element_text(size=15),title = element_text(size = 20), legend.text = element_text(size = 15),axis.text= element_text(size = 12))
#print plot
poc_density_plot

#plot histogram of Pow parasite density
pow_density_plot <- ggplot(pow_dens, aes(x=pow_parasite_density, fill=status)) + 
  # geom_histogram(aes(y=after_stat(density)),show.legend = T) + 
  geom_density(aes(y=after_stat(density)), show.legend = T, alpha = 0.5) +
  theme_bw() + 
  scale_fill_manual(values = c("Poc-Pow mixed infection" = "#E6AB02","Poc mono-infection" = '#1B9E77', "Pow mono-infection" = '#7570B3'), drop = FALSE)+
  labs(x="Pow parasite density per microliter", y = "Proportion",fill=element_blank()) + 
  theme(plot.title = element_text(hjust =0.5)) + 
  scale_x_log10() +
  theme_bw()+ 
  theme(axis.title.y = element_text(size = 15),axis.title.x = element_text(size=15),title = element_text(size = 20), legend.text = element_text(size = 15),axis.text= element_text(size = 12))
#print plot
pow_density_plot


#arrange Poc and Pow plots side-by-side
dual_po_density_plot <- ggarrange(poc_density_plot,pow_density_plot,ncol =2,common.legend = TRUE,legend="right")
dual_po_density_plot

#save plot
ggsave('Figures/poc-pow_parasite_density_plots.png', dual_po_density_plot, height = 6, width = 12, dpi = 600)
