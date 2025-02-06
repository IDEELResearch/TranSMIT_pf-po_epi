################################################
# Plot Po species-identification status by Po18S Cyclic threshold tertiles
# Author: Kelly Carey-Ewend 
# Date 10 November 2024
# Modified: 12 January 2025
################################################

#load necessary packages
suppressPackageStartupMessages({
  library(eulerr)
})

#set working directory
setwd("~/Documents/P. ovale Epi Bagamoyo")

#assign group names for all venn diagrams
#Pf-positive, Po-positive (and received species-identification assay), Poc-positive (containing the species per algorithm), Pow-positive (containing the species per algorithm), U (all Po-positives, including those that were not tested for species composition)
#the group combinations ensure that all Poc and Pow samples fall within the Po (Po species-tested) umbrella, and all Po samples that were tested for species fall within the U umbrella (which also includes untested Po-positives)
group <- c("Poc&Po&U","Pow&Po&U", "Poc&Pow&Po&U", "Pf", "Pf&Poc&Po&U","Pf&Pow&Po&U","Pf&Poc&Pow&Po&U","Pf&Po&U","Po&U", "U", "U&Pf")

#manually enter observed counts in all groups
counts <- c(
  #Poc&Po
  150, 
  #Pow&Po
  83,
  #Poc&Pow&Po
  74,
  #Pf
  1812,
  #Pf&Poc&Po
  39,
  #Pf&Pow&Po
  17,
  #Pf&Poc&Pow&Po
  13,
  #Pf&Po
  40,
  #Po
  215,
  #Po no assay run
  149,
  #Po no assay run & Pf
  47)

#assign names to each count for plotting
names(counts) <- group

#set random seed
set.seed(1)

#create observed counts venn object. Ellipse reshaping circles to match data distribution
venn<- euler(counts, shape = "ellipse")

#create and save plot
png('Figures/co-infection_venn.png', width = 5.5, height = 5.5, units = "in", res = 700)
obs_plot <- plot(venn, fills = list(fill = c('lightgreen','#D95F02','grey','orchid1','skyblue')), labels = list(font =3))
obs_plot
dev.off()


#Determine expected counts assuming independent assortment

#Marginal probabilities (except Pf) among 6977 (excluding 196 w/o species-identification assays).
#Pf probability was calculated across total dataset, but only applied to 6977 for expected counts
#Po probability refers to Po-positive with assay run; untested samples are treated as a separate group

#P(Pf)
pf = 1968/7173
#P(Poc)
poc = 276/6977
#P(Pow)
pow = 187/6977
#P(Po) is calculated only including samples potentially tested for species composition
po = 631/6977

#enter total sample size (not including untested Po-positives, which are not assumed to be independently assorted; testing status was based on laboratory logistics)
total = 6977

#calculate expected counts for each group (except untested Po and untested Po&Pf)
count_pf_poc_pow = total*pf*poc*pow
count_pf_pow = total*pf*pow - count_pf_poc_pow
count_pf_poc = total*pf*poc - count_pf_poc_pow
count_pf_po = total*pf*po - count_pf_poc_pow - count_pf_poc - count_pf_pow
count_poc_pow = total*poc*pow - count_pf_poc_pow
count_poc = total*poc - count_poc_pow - count_pf_poc - count_pf_poc_pow
count_pow = total*pow - count_poc_pow - count_pf_pow - count_pf_poc_pow
count_po = total*po - count_poc - count_pow - count_poc_pow - count_pf_po - count_pf_poc - count_pf_pow - count_pf_poc_pow 
count_pf = total*pf -count_pf_po - count_pf_poc - count_pf_pow - count_pf_poc_pow

#create expected data counts; note that Po and untested Po&Pf are the same values from observed dataset
counts_expected <- c(
  #Poc&Po
  count_poc,
  #Pow&Po
  count_pow,
  #Poc&Pow&Po
  count_poc_pow,
  #Pf
  count_pf,
  #Pf&Poc&Po
  count_pf_poc,
  #Pf&Pow&Po
  count_pf_pow,
  #Pf&Poc&Pow&Po
  count_pf_poc_pow,
  #Pf&Po
  count_pf_po,
  #Po
  count_po,
  #Po untested
  149,
  #Po untested&Pf
  47)

#assign names to each count for plotting
names(counts_expected) <- group

#set random seed
set.seed(1)

#create expected counts venn object
venn_expected <- euler(counts_expected, shape = "ellipse")

#save plot
png('Figures/co-infection_venn_expected.png', width = 5.5, height = 5.5, units = "in", res = 700)
plot(venn_expected, fills = list(fill = c('lightgreen','#D95F02','grey','orchid1','skyblue')), labels = list(font =3))
dev.off()
