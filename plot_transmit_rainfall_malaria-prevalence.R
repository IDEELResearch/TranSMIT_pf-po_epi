################################################
# Plot TranSMIT parasite prevalence by month
# Author: Julia Muller and Kelly Carey-Ewend 
# Date 8 November 2024
# Modified: 15 January 2025
################################################

# Load packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(gridExtra)
  library(ggpp)
  library(zoo)
  library(dplyr)
})


#set working directory
setwd('~/Documents/P. ovale Epi Bagamoyo')

# Import data
prev <- read.csv('transmit_rainfall_malaria-prev_by-month.csv')
rainfall <- read.csv('bagamoyo_rainfall.csv')

#rename columns for plotting
prev$month <- prev$study_month

prev$total_observations <- prev$N


#---------------#
# Functions ----
#---------------#
generate_observations_table <- function(data, cols, base_size = 8) {
  tableGrob(subset(data, select = c(month, total_observations)), 
            rows = NULL, 
            cols = c('Month', 'Total Screened'),
            ttheme_gtbw(base_size = base_size,
                        colhead = list(bg_params = list(fill = '#f2f2f2')),
                        core = list(bg_params = list(fill = c('white', '#fafafa')))))
}

add_observations_table <- function(plot, table) { 
  grid.arrange(grobs = list(plot, table), ncol = 2, widths = list(2, 0.35)) 
}

add_annotations <- function(label_height) {
  # Helper function to create annotation elements for a given range
  annotations <- function(label_height, month_start, month_end) {
    list(ggplot2::annotate('text', x = (month_start+month_end)/2, y = label_height, label = 'No Screening', color = 'black', size = 3),
         ggplot2::annotate('segment', x = month_start, xend = month_end, y = label_height * 0.95, yend = label_height * 0.95, color = 'black', size = 0.3),
         ggplot2::annotate('segment', x = month_start, xend = month_start, y = label_height * 0.95, yend = label_height * 0.90, color = 'black', size = 0.3),
         ggplot2::annotate('segment', x = month_end, xend = month_end, y = label_height * 0.95, yend = label_height * 0.90, color = 'black', size = 0.3),
         ggplot2::annotate('segment', x = month_start, xend = month_start, y = 0, yend = label_height * 0.90, color = 'black', size = 0.3, linetype = 'dashed'),
         ggplot2::annotate('segment', x = month_end, xend = month_end, y = 0, yend = label_height * 0.90, color = 'black', size = 0.3, linetype = 'dashed'),
         ggplot2::annotate('segment', x = (month_start+month_end)/2, xend = (month_start+month_end)/2, y = label_height * 0.95, yend = label_height * 0.97, color = 'black', size = 0.3)
    )
  }
  # Create annotations for all specified date ranges
  list(annotations(label_height, 3, 6),         # Dec 2018 - Jan 2019
       annotations(label_height, 11, 12),       # Aug 2019 - Sep 2019
       annotations(label_height, 15.75, 16.25), # Jan 2020
       annotations(label_height, 19, 24),       # Apr 2020 - Sep 2020
       annotations(label_height, 40, 41)        # Jan 2022 - Feb 2022
  )
}


#------------------#
# Prepare data ----
#------------------#
# Merge prevalence and rainfall data 
prev_rain <- full_join(prev, rainfall, by = 'month')

# Generate a list of all study months (October 2018 to June 2022)
study_months <- as.character(seq(as.yearmon(as.character(201810), '%Y%m'),
                                as.yearmon(as.character(202206), '%Y%m'),
                                1/12))

# Replace months without screening with '...' to indicate absent data
observations <- prev %>%
  select(month, total_observations) %>%
  mutate(
    month = study_months,
    total_observations = as.character(total_observations),
    across(everything(), ~ if_else(is.na(total_observations), '...', .)), # Replace NA with '...'
    flag = c(TRUE, (total_observations[-1] != '...' | total_observations[-length(total_observations)] != '...'))) %>%
  filter(flag) %>% # Keep only the first occurrence of '...'
  select(-flag)

# Create a table of number of monthly observations for the entire study period
obs_table <- generate_observations_table(observations)


###Calculate Pearson's correlation coefficient between monthly Po and Pf prevalences
prev_nomissing <- na.omit(prev)
cor.test(prev_nomissing$pf_prev, prev_nomissing$po_prev, method = "pearson")
cor.test(prev_nomissing$poc_prev, prev_nomissing$pow_prev, method = "pearson")

#--------------------#
# Generate plots ----
#--------------------#

#determine color palette
library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)
display.brewer.pal(n=6,"Dark2")
brewer.pal(n=6,"Dark2")

#load package for variable reordering
library("forcats")

#create set of columns
cols = c('pf_prev','po_prev','pfpo_prev')

# P. falciparum and P. ovale plot
plot_pf_po <- prev_rain %>%
  pivot_longer(cols = c('pf_prev','po_prev','pfpo_prev'), 
               names_to = 'variable', values_to = 'value') %>%
  mutate(variable = fct_relevel(variable, 
                                'pf_prev','po_prev','pfpo_prev')) %>%
  ggplot() +
    geom_area(aes(x = month, y = CHIRPS_rainfall/10, group = 1), color = 'lightblue3', fill = 'lightblue', alpha = 0.75) + 
    geom_bar(aes(x = month, y = value, fill = variable), stat = 'identity', position = 'dodge', color = 'black', linewidth = 0.1) +
    scale_fill_manual(name = NULL,
                      labels = c('pf_prev'='P. falciparum', 'po_prev'='P. ovale', 'pfpo_prev'='Pf-Po co-infection', 'poc_prev'='Poc mono-inf.', 'pow_prev'='Pow mono-inf.','mixed_prev'='Poc-Pow mixed inf.'), 
                      values = c('pf_prev'='darkblue', 'po_prev'='#D95F02', 'pfpo_prev'='grey50', 'poc_prev'='#1B9E77', 'pow_prev'='#7570B3','mixed_prev'='#E6AB02')) +
    theme_light() + 
    theme(legend.position = 'top',
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    guides(fill = guide_legend(nrow = 1)) +
    labs(x = NULL, y = 'Percent positive (%)') +
    scale_x_discrete(labels = study_months) + 
    scale_y_continuous(sec.axis = sec_axis(~.*10, name = 'Rainfall (mm)'), breaks = seq(0, 60, 10)) +
    add_annotations(label_height = 52)

#Poc and Pow plot
plot_poc_pow <- prev_rain %>%
  pivot_longer(cols = c('poc_prev','pow_prev','mixed_prev'), 
               names_to = 'variable', values_to = 'value') %>%
  mutate(variable = fct_relevel(variable, 
                                'poc_prev','pow_prev','mixed_prev')) %>%
  ggplot() +
  geom_area(aes(x = month, y = CHIRPS_rainfall/35, group = 1), color = 'lightblue3', fill = 'lightblue', alpha = 0.75) + 
  geom_bar(aes(x = month, y = value, fill = variable), stat = 'identity', position = 'dodge', color = 'black', linewidth = 0.1) +
  scale_fill_manual(name = NULL,
                    labels = c('pf_prev'='P. falciparum', 'po_prev'='P. ovale', 'pfpo_prev'='Pf-Po co-infection', 'poc_prev'='Poc mono-inf.', 'pow_prev'='Pow mono-inf.','mixed_prev'='Poc-Pow mixed inf.'), 
                    values = c('pf_prev'='darkblue', 'po_prev'='#D95F02', 'pfpo_prev'='grey50', 'poc_prev'='#1B9E77', 'pow_prev'='#7570B3','mixed_prev'='#E6AB02')) +
  theme_light() + 
  theme(legend.position = 'top',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  guides(fill = guide_legend(nrow = 1)) +
  labs(x = NULL, y = 'Percent positive (%)') +
  scale_x_discrete(labels = study_months) + 
  scale_y_continuous(sec.axis = sec_axis(~.*35, name = 'Rainfall (mm)'), breaks = seq(0, 15, 3)) +
  add_annotations(label_height = 14)

#add sampling table to pf-po plot
plot_pf_po_table <- add_observations_table(plot_pf_po, obs_table)

#save pf-po plot
ggsave('Figures/prev_by_month_pf_po.png', plot_pf_po_table, height = 4, width = 11, dpi = 300)

#add sampling table to pf-po plot
plot_poc_pow_table <- add_observations_table(plot_poc_pow, obs_table)

#save pf-po plot
ggsave('Figures/prev_by_month_poc_pow.png', plot_poc_pow_table, height = 4, width = 11, dpi = 300)

#load package for arranging plots within same image
library(ggpubr)

#arrange prior plots together
pfpo_pocpow <- ggarrange(plot_pf_po, plot_poc_pow, nrow=2)

#add single sampling table
pfpopocpow_plot <- add_observations_table(pfpo_pocpow, obs_table)


#save plot
ggsave('Figures/prev_by_month_pf_po_poc_pow.png', pfpopocpow_plot, height = 8, width = 11, dpi = 300)






# P. falciparum, P. ovale, and P. malariae plot
plot_pf_po_pm <- prev_rain %>%
  pivot_longer(cols = c('pf_prevalence', 'po_prevalence', 'po_no_pf_prevalence', 'pm_prevalence'), 
               names_to = 'variable', values_to = 'value') %>%
  ggplot() +
    geom_area(aes(x = month, y = rainfall/6.67, group = 1), color = 'lightblue3', fill = 'lightblue', alpha = 0.75) + 
    geom_bar(aes(x = month, y = value*100, fill = variable), stat = 'identity', position = 'dodge', color = 'black', linewidth = 0.1) +
    scale_fill_manual(name = NULL,
                      labels = c('P. falciparum', 'P. malariae', 'P. ovale (mono-infection)', 'P. ovale (all)'),
                      values = c('darkblue', 'seagreen3', 'grey50', 'orange')) +
    theme_light() + 
    theme(legend.position = 'top',
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(x = NULL, y = 'Percent positive (%)') +
    scale_x_discrete(labels = study_months) + 
    scale_y_continuous(sec.axis = sec_axis(~.*6.67, name = 'Rainfall (mm)'), breaks = seq(0, 60, 10)) +
    add_annotations(label_height = 52)
plot_pf_po_pm_table <- add_observations_table(plot_pf_po_pm, obs_table)

# Export figures
ggsave('figures/prev_by_month_pf_po.png', plot_pf_po_table, height = 8, width = 12, dpi = 300)
ggsave('figures/prev_by_month_pf_po_pm.png', plot_pf_po_pm_table, height = 8, width = 12, dpi = 300)
