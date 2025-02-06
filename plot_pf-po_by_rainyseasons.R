################################################
# Plot P. ovale mono- vs. co-infection by season
# Author: Julia Muller, modified by Kelly Carey-Ewend
# Date: 13 January 2025
# Last modified: 27 January 2025
################################################

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggtext)
  library(ggpubr)
})

# Import data
setwd('~/Documents/P. ovale Epi Bagamoyo')
data <- read_csv('transmit_rainfall_malaria-prev_studied-months_lag6.csv', show_col_types = F)

# Classify short/long wet season by calendar month, filter out dry season
data <- data %>%
  separate(study_month_lag6, into = c('year', 'month'), sep = '-', remove = F) %>%
  mutate(across(year:month, as.numeric)) %>%
  mutate(season = case_when(
    month %in% c(3, 4, 5) ~ 'Long Wet',
    month %in% c(10, 11, 12) ~ 'Short Wet',
    TRUE ~ 'Dry')) %>%
  filter(season != 'Dry')
  


#----------------------------------------#
# Plot yearly short/long wet seasons ----
#----------------------------------------#
# Summarize by year and season
data_yearly <- data %>%
  group_by(year, season) %>%
  summarize(
    po_mono_prev = sum(po_count, na.rm = T)/sum(N),
    po_co_prev = sum(pfpo_count, na.rm = T)/sum(N),
    total_screened = sum(N)) %>%
  ungroup()

# Pivot data to long format for plotting
data_lf <- data_yearly %>%
  pivot_longer(cols = c('po_mono_prev', 'po_co_prev'), names_to = 'variable', values_to = 'value')

# Select long wet season; add in blank data for 2018 and 2020 for cross-season comparisons
data_long_wet <- data_lf %>% 
  filter(season == 'Long Wet') %>%
  bind_rows(expand.grid(
    year = 2018,
    season = c('Short wet', 'Long wet'),
    variable = c('po_mono_prev', 'po_co_prev'),
    value = 0)) %>%
  bind_rows(expand.grid(
    year = 2020,
    season = c('Short Wet', 'Long Wet'),
    variable = c('po_mono_prev', 'po_co_prev'),
    value = 0))

# Select short wet season; add in blank data for 2022 for cross-season comparisons
data_short_wet <- data_lf %>% 
  filter(season == 'Short Wet') %>%
  bind_rows(expand.grid(
    year = 2022,
    season = c('Short Wet', 'Long Wet'),
    variable = c('po_mono_prev', 'po_co_prev'),
    value = 0))

# Create new year+sample size variable for plotting
data_long_wet$year.n <- paste(data_long_wet$year,"\n(n=",data_long_wet$total_screened,")", sep = "")
data_short_wet$year.n <- paste(data_short_wet$year,"\n(n=",data_short_wet$total_screened,")", sep = "")


# Plot long wet seasons, by year
plot_long_wet <- ggplot(data = data_long_wet) +
  geom_bar(aes(x = year.n, y = value, fill = variable), stat = 'identity', position = 'dodge') +
  labs(x = 'Year of study', y = expression(paste(italic('Masika'),' prevalence')), fill = NULL) + 
  scale_fill_manual(values = c('grey50', '#D95F02'), labels = c('*Pf-Po* co-infection', '*Po* mono-infection')) +
  theme_bw(base_size = 14) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.7,0.8),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.text.x = element_text(face="bold"), 
        axis.text.y = element_text(size=14),
        legend.text = element_markdown()) +
  scale_y_continuous(limits = c(0, 0.25))
    

# Plot short wet seasons, by year
plot_short_wet <- ggplot(data = data_short_wet) +
  geom_bar(aes(x = year.n, y = value, fill = variable), stat = 'identity', position = 'dodge') +
  labs(x = 'Year of study', y = expression(paste(italic('Vuli'),' prevalence')), fill = NULL) + 
  scale_fill_manual(values = c('grey50', '#D95F02'), labels = c('*Po-Pf* co-infection', '*Po* mono-infection')) +
  theme_bw(base_size = 14) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.7,0.8),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(face="bold"), axis.text.y = element_text(size=14), legend.text = element_markdown()) +
  scale_y_continuous(limits = c(0, 0.25))

# Arrange short and long wet season plots 
combined_plot <- ggarrange(plot_long_wet, plot_short_wet, nrow = 2)
ggsave('Figures/po_mono_coinfection_yearly.png', combined_plot, height = 6.5, width = 4.5, dpi = 600)


#---------------------------------------------#
# Plot mono-infection/co-infection ratios ----
#---------------------------------------------#
# Summarize by month
data_monthly <- data %>%
  group_by(study_month_lag6, season) %>%
  summarize(po_ratio = sum(po_count, na.rm = T)/sum(pfpo_count, na.rm = T)) %>%
  filter(!is.na(po_ratio), po_ratio != Inf)

# Plot mono-infection to co-infection ratio
plot_ratio <- ggplot(data = data_monthly, aes(x = season, y = po_ratio)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.2) + 
  labs(x = 'Season', y = "*Po* mono-infection to *Po*-*Pf* co-infection ratio") +
  theme_bw(base_size = 16) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  scale_x_discrete(label=c("Short Wet" = "*Vuli* rains", "Long Wet" = "*Masika* rains")) +
  scale_y_log10(breaks = c(0,1,2,4,8,16,32)) +
  #coord_flip() +
  theme(axis.title = element_markdown(size = 16), axis.text.y = element_markdown(size=16), axis.text.x= element_markdown(size=16))
ggsave('Figures/po_mono_coinfection_ratio.png', plot_ratio, height = 5.5, width = 4, dpi = 600)
library(devtools)
devtools::install_github("kkeenan02/MsatAllele")
library(MsatAllele)
