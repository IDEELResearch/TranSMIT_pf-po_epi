The de-identified dataset in supplemental materials can be processed into a cleaned SAS
dataset using the script contained in P._ovale_epi_bagamoyo_datacleaning_KellyCE.sas.
Of note, the published dataset does not contain the enrollment_village variable for anonymity.
The script should be adjusted to ignore this input variable, or a full dataset can be requested from the authors.

The P._ovale_epi_bagamoyo_analysis_KellyCE.sas script can then be run on this cleaned dataset
to do additional variable processing and generate final analyses.

Some data visualizations require output of data from the SAS environment to an R
environment followed by processing with the included R scripts:

plot_pf-po_antagonism_forestplot.R
plot_pf-po_by_rainyseasons.R
plot_pf-po-poc-pow_paras-density.R
transmit_plasmodiuminteraction_NiDP.R
plot_poc-pow_obs-exp_venn.R
plot_pf-po-poc-pow_wet-prev-by-year.R
plot_po-species_po18s-ct-tertiles.R
plot_pf-po_village-prev.R
plot_transmit_rainfall_malaria-prevalence.R
NiDP_modelFitting.R