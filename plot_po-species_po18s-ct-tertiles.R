################################################
# Plot Po species-identification status by Po18S Cyclic threshold tertiles
# Author: Kelly Carey-Ewend 
# Date 10 November 2024
# Modified: 12 January 2025
################################################

#load necessary packages
suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(tidyverse)
  library(reshape2)
})

#set working directory
setwd("~/Documents/P. ovale Epi Bagamoyo")

#load in data table
status_density <- read.csv("Po-status_by_Po-density.csv", header = TRUE)

#create long version of dataset, with a new observation per Ct tertile and species combination
status_density_long <- status_density %>% pivot_longer(!Po.species, names_to = "density", values_to = "prop")

#reorder density for plotting
status_density_long$species <-fct_relevel(status_density_long$Po.species, 'Poc mono-inf. (n=198)','Pow mono-inf. (n=104)','Poc-Pow mixed inf. (n=91)', 'Po species unknown (n=268)')

#create long version of ordered density
status_density_long$dens_ordered <-fct_relevel(status_density_long$density, 'High.Po.density','Medium.Po.density',"Low.Po.density")

#create plot
density_plot <- ggplot(status_density_long, aes(fill=dens_ordered, y=prop, x=species)) + 
  geom_bar(position="stack", stat="identity") +
  #manually assign color to each density group, and change label to be italicized
  scale_fill_manual (values = c('High.Po.density'= '#E7298A','Medium.Po.density'='#E6AB02',"Low.Po.density"='grey50'),
                     labels = c('High.Po.density'= expression(paste("High ",italic("Po")," density")),
                                'Medium.Po.density'=expression(paste('Medium ',italic('Po'),' density')),
                                "Low.Po.density"=expression(paste("Low ",italic("Po")," density")))) +
  labs(y="Proportion", x=expression(paste(italic("Po")," species composition"))) +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 40, vjust = 0.5)) +
  #rename x axis species labels to feature sample size on next row
  scale_x_discrete(label = c(
    "Poc mono-inf. (n=195)"=expression(atop(paste(italic("Poc")," mono-inf."),"(n=195)")),
    "Pow mono-inf. (n=103)"=expression(atop(paste(italic("Pow")," mono-inf."),"(n=103)")),
    "Poc-Pow mixed inf. (n=90)"=expression(atop(paste(italic("Poc"),"-",italic("Pow")," mixed inf."),"(n=90)")),
    'Po species unknown (n=263)' = expression(atop(paste(italic('Po'),' species unknown'),'(n=263)'))))

#print plot
density_plot

#save plot
ggsave('Figures/species-status_by_po-density.png', density_plot, height = 4, width = 5, dpi = 600)
