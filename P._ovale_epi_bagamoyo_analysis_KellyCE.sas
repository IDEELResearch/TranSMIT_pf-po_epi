DM "log; clear; odsresults; clear;";
/*****************************************************************************
	Name: Kelly Carey-Ewend
	Program: P._ovale_epi_bagamoyo_analysis_KellyCE.sas
   	Date: 06-05-2024
	Description: SAS program file for analyzing P. ovale prevalence and distribution in Bagamoyo in the TranSMIT study
*****************************************************************************/
OPTIONS MERGENOBY=warn NODATE NONUMBER FORMCHAR="|----|+|---+=|-/\<>*";
FOOTNOTE "P._ovale_epi_bagamoyo_analysis_KellyCE.sas run at %SYSFUNC(DATETIME(), DATETIME.) by Kelly Carey-Ewend";
/******************************* begin program ******************************/
TITLE;
LIBNAME dir '\\Mac\Home\Documents\P. ovale Epi Bagamoyo';

*format variables;
proc format;
	value enrollment_gender_ 0='Male' 1='Female';
	value enrollment_antimalarial_yn_ 0='No' 1='Yes' 
		98='Don''t know' 99='Refused';
	value season_ 0='Wet' 1='Dry';
	value region_ 0='Northwest' 1='South';
	value po18s_tert_ 2='Low Density' 1='Medium Density' 0='High Density';
	value pf18s_tert_ 2='Low Density' 1='Medium Density' 0='High Density';
	run;


*Import cleaned dataset;
DATA transmit_cleaned;
	SET dir.transmit_cleaned;
RUN;

PROC CONTENTS DATA=transmit_cleaned;
RUN;

*** Functional Form Assessment;
*evaluate distributions of continuous variables;
PROC MEANS DATA=transmit_cleaned P25 P50 P75 MAX MIN STDDEV MEAN;
	VAR enrollment_age;
	VAR rainfall_avg_lag0_past_1mo;
	VAR pf18s_ct;
	VAR po18s_ct;
RUN;

*Determine tertiles of rainfall distribution;
PROC UNIVARIATE DATA=transmit_cleaned;
	VAR rainfall_avg_lag0_past_1mo;
	OUTPUT out=rainfall_tertiles pctlpts=33.33 66.67 pctlpre=rainfall pctlname = pct33 pct66;
RUN;

PROC PRINT DATA=rainfall_tertiles;
RUN;

*Determine tertiles of po18S Ct values;
PROC UNIVARIATE DATA=transmit_cleaned;
	VAR po18s_ct;
	OUTPUT out=po18s_ct_tertiles pctlpts=33.33 66.67 pctlpre=po_ct_ pctlname = pct33 pct66;
RUN;

PROC PRINT DATA=po18s_ct_tertiles;
RUN;

*Determine tertiles of pf18S Ct values;
PROC UNIVARIATE DATA=transmit_cleaned;
	VAR pf18s_ct;
	OUTPUT out=pf18s_ct_tertiles pctlpts=33.33 66.67 pctlpre=pf_ct_ pctlname = pct33 pct66;
RUN;

PROC PRINT DATA=pf18s_ct_tertiles;
RUN;

*Determine how many samples were at risk for Pf-Po cross-reactivity;
*Pf+ samples with >100 parasites/uL had chance of causing false-positive Po;
PROC FREQ DATA=transmit_cleaned;
	WHERE qpcr_pfdens_screening > 100 & po18s_ct > 42;
	TABLES pf*po;
RUN;

PROC PRINT DATA=transmit_cleaned;
	WHERE qpcr_pfdens_screening > 100 & po18s_ct > 42;
RUN;

*Code new forms of key variables to be investigated;
DATA transmit_cleaned;
	SET transmit_cleaned;
	*implement new form of po positivity with cutoff at 40 cycles instead of 45;
	po_40 = 0;
	IF . < po18s_ct < 40 THEN po_40 = 1;
	*impute suspected Po parasite density from a po18S standard dilution (see supplemental data);
	po_parasite_density = 2.71828**((po18s_ct - 44.15)/(-1.532));
	*Define submicroscopic Pf infection as <100p/uL;
	pf_submicroscopic = .;
	IF pf = 1 THEN DO;
		pf_submicroscopic = 0;
		IF . < qpcr_pfdens_screening < 100 THEN pf_submicroscopic = 1;
		END;
	*Define submicroscopic Po infection as <100p/uL;
	po_submicroscopic = .;
	IF po = 1 THEN DO;
		po_submicroscopic = 0;
		IF . < po_parasite_density < 100 THEN po_submicroscopic = 1;
		END;
	*identify Po mono-infections, Pf mono-infections, and Po-Pf mixed infections;
	IF po = 1 & pf = 0 THEN po_mono = 1;
		ELSE po_mono = 0;
	IF pf = 1 & po = 0 THEN pf_mono = 1;
		ELSE pf_mono = 0;
	IF pf = 1 & po = 1 THEN popf_co = 1;
		ELSE popf_co = 0;
	*Determine appropriate forms for age in analysis:
	*Consider age quartiles, binary child-adult, restricted quadratic splines, restircted cubic splines;
	*first, encode squared and cubed versions of age in years;
	sqage = enrollment_age**2;
	cubage = enrollment_age**3;
	*Age quartiles are 11, 18, 30;
	*encode age in years over each quartile;
	age11 = enrollment_age - 11;
		IF age11 < 0 THEN age11 = 0;
	age18 = enrollment_age - 18;
		IF age18 < 0 THEN age18 = 0;
	age30 = enrollment_age - 30;
		IF age30 < 0 THEN age30 = 0;
	*create squared age in excess of each quartile;
	sqage11 = age11**2;
	sqage18 = age18**2;
	sqage30 = age30**2;
	*create restricted quadratic splines at each breakpoint;
	rqspline11 = sqage11 - sqage30;
	rqspline18 = sqage18 - sqage30;
	*Literature age breakpoints at 10, 15, above;
	age10 = enrollment_age - 10;
		IF age10 < 0 THEN age10 = 0;
	age15 = enrollment_age - 15;
		IF age15 < 0 THEN age15 = 0;
	sqage10 = age10**2;
	sqage15 = age15**2;
	rqspline10 = sqage10 - sqage15;
	*Code age categorically by quartiles;
	IF enrollment_age ^= . THEN DO; agecat11 = 0; agecat18 = 0; agecat30= 0; agecatover= 0; END;
	IF 0 < enrollment_age <= transmit_quartiles THEN agecat11 = 1;
		ELSE IF 11 < enrollment_age <= 18 THEN agecat18 = 1;
		ELSE IF 18 < enrollment_age <= 30 THEN agecat30 = 1;
		ELSE IF enrollment_age > 30 THEN agecatover = 1;
	*Add categorical age divisions at literature breakpoints (ages 10, 15, above);
	IF enrollment_age ^= . THEN DO; agecat10 = 0; agecat15 = 0; agecat16up= 0; END;
	IF 0 < enrollment_age <= 10 THEN agecat10 = 1;
		ELSE IF 10 < enrollment_age <= 15 THEN agecat15 = 1;
		ELSE IF 15 < enrollment_age THEN agecat16up = 1;
	*Add single categorical variable for age quartiles (DONT USE IN REGRESSIONS);
	IF agecat10 = 1 THEN agecat = 0;
		ELSE IF agecat15 = 1 THEN agecat = 1;
		ELSE IF agecat16up = 1 THEN agecat = 2;
		ELSE agecat = .;
	*create interaction terms for age categories and pf status;
	*one of these interactions terms doesn't need to be included in regressions because will be colinear with the other terms being 0;
	agecat10pf = agecat10*pf;
	agecat15pf = agecat15*pf;
	agecat16uppf = agecat16up*pf;
	*add age bins by years of 5 for use in indirect standardization;
	age5bin = floor(enrollment_age/5);
	*Add enrollment_gender and pf interaction term, which is 1 for females with positive pf;
	sexpf = enrollment_gender*pf;
	*code rainfall quadratic and categorical terms;
	rainsq = rainfall_avg_lag0_past_1mo**2;
	*rainfall quartiles: breakpoints must be coded manually;
	IF rainfall_avg_lag0_past_1mo ^= . THEN DO; raincat0 = 0; raincat1 = 0; raincat2= 0; raincatover= 0; END;
	IF 0 <= rainfall_avg_lag0_past_1mo <= 1.5261 THEN raincat0 = 1;
		ELSE IF 1.5261 < rainfall_avg_lag0_past_1mo <= 3.2954000 THEN raincat1 = 1;
		ELSE IF 3.2954000 < rainfall_avg_lag0_past_1mo <= 6.2858330 THEN raincat2 = 1;
		ELSE IF 6.2858330 < rainfall_avg_lag0_past_1mo THEN raincatover = 1;
	*rainfall tertiles;
	IF 0 <= rainfall_avg_lag0_past_1mo <= 2.08753 THEN raintert0 = 1;
		ELSE IF 2.08753 < rainfall_avg_lag0_past_1mo <= 5.11943 THEN raintert1 = 1;
		ELSE IF 5.11943 < rainfall_avg_lag0_past_1mo THEN raintert2 = 1;
	*Add single categorical variable for rain quartiles;
	IF raincat0 = 1 THEN raincat = 0;
		ELSE IF raincat1 = 1 THEN raincat = 1;
		ELSE IF raincat2 = 1 THEN raincat = 2;
		ELSE IF raincatover = 1 THEN raincat = 3;
		ELSE raincat = .;
	*add categorical variables for rainfall tertiles;
	IF raintert0 = 1 THEN raintert = 0;
		ELSE IF raintert1 = 1 THEN raintert = 1;
		ELSE IF raintert2 = 1 THEN raintert = 2;
		ELSE raintert = .;
	*code rainfall with 6 week lag (see above);
	rainlag6sq = rainfall_avg_lag6_past_1mo**2;
	IF rainfall_avg_lag6_past_1mo ^= . THEN DO; rainlag6cat0 = 0; rainlag6cat1 = 0; rainlag6cat2= 0; rainlag6catover= 0; END;
	IF 0 <= rainfall_avg_lag6_past_1mo <= 1.5261 THEN rainlag6cat0 = 1;
		ELSE IF 1.5261 < rainfall_avg_lag6_past_1mo <= 3.2954 THEN rainlag6cat1 = 1;
		ELSE IF 3.2954 < rainfall_avg_lag6_past_1mo <= 6.285833 THEN rainlag6cat2 = 1;
		ELSE IF 6.285833 < rainfall_avg_lag6_past_1mo THEN rainlag6catover = 1;
	IF 0 <= rainfall_avg_lag6_past_1mo <= 2.08753 THEN rainlag6tert0 = 1;
		ELSE IF 2.08753 < rainfall_avg_lag6_past_1mo <= 5.11943 THEN rainlag6tert1 = 1;
		ELSE IF 5.11943 < rainfall_avg_lag6_past_1mo THEN rainlag6tert2 = 1;
	*Add single categorical variable for rain quartiles and tertiles;
	IF rainlag6cat0 = 1 THEN rainlag6cat = 0;
		ELSE IF rainlag6cat1 = 1 THEN rainlag6cat = 1;
		ELSE IF rainlag6cat2 = 1 THEN rainlag6cat = 2;
		ELSE IF rainlag6catover = 1 THEN rainlag6cat = 3;
		ELSE rainlag6cat = .;
	IF rainlag6tert0 = 1 THEN rainlag6tert = 0;
		ELSE IF rainlag6tert1 = 1 THEN rainlag6tert = 1;
		ELSE IF rainlag6tert2 = 1 THEN rainlag6tert = 2;
		ELSE rainlag6tert = .;
	*code rainfall with 9 week lag (see above);
	rainlag9sq = rainfall_avg_lag9_past_1mo**2;
	IF rainfall_avg_lag9_past_1mo ^= . THEN DO; rainlag9cat0 = 0; rainlag9cat1 = 0; rainlag9cat2= 0; rainlag9catover= 0; END;
	IF 0 <= rainfall_avg_lag9_past_1mo <= 1.5261 THEN rainlag9cat0 = 1;
		ELSE IF 1.5261 < rainfall_avg_lag9_past_1mo <= 3.2954 THEN rainlag9cat1 = 1;
		ELSE IF 3.2954 < rainfall_avg_lag9_past_1mo <= 6.285833 THEN rainlag9cat2 = 1;
		ELSE IF 6.285833 < rainfall_avg_lag9_past_1mo THEN rainlag9catover = 1;
	IF 0 <= rainfall_avg_lag9_past_1mo <= 2.08753 THEN rainlag9tert0 = 1;
		ELSE IF 2.08753 < rainfall_avg_lag9_past_1mo <= 5.11943 THEN rainlag9tert1 = 1;
		ELSE IF 5.11943 < rainfall_avg_lag9_past_1mo THEN rainlag9tert2 = 1;
	*Add single categorical variable for rain quartiles;
	IF rainlag9cat0 = 1 THEN rainlag9cat = 0;
		ELSE IF rainlag9cat1 = 1 THEN rainlag9cat = 1;
		ELSE IF rainlag9cat2 = 1 THEN rainlag9cat = 2;
		ELSE IF rainlag9catover = 1 THEN rainlag9cat = 3;
		ELSE rainlag9cat = .;
	*Add single categorical variables for rain tertiles;
	IF rainlag9tert0 = 1 THEN rainlag9tert = 0;
		ELSE IF rainlag9tert1 = 1 THEN rainlag9tert = 1;
		ELSE IF rainlag9tert2 = 1 THEN rainlag9tert = 2;
		ELSE rainlag9tert = .;
	*code rainfall with 12 week lag (see above);
	rainlag12sq = rainfall_avg_lag12_past_1mo**2;
	IF rainfall_avg_lag12_past_1mo ^= . THEN DO; rainlag12cat0 = 0; rainlag12cat1 = 0; rainlag12cat2= 0; rainlag12catover= 0; END;
	IF 0 <= rainfall_avg_lag12_past_1mo <= 1.5261 THEN rainlag12cat0 = 1;
		ELSE IF 1.5261 < rainfall_avg_lag12_past_1mo <= 3.2954 THEN rainlag12cat1 = 1;
		ELSE IF 3.2954 < rainfall_avg_lag12_past_1mo <= 6.285833 THEN rainlag12cat2 = 1;
		ELSE IF 6.285833 < rainfall_avg_lag12_past_1mo THEN rainlag12catover = 1;
	IF 0 <= rainfall_avg_lag12_past_1mo <= 2.08753 THEN rainlag12tert0 = 1;
		ELSE IF 2.08753 < rainfall_avg_lag12_past_1mo <= 5.11943 THEN rainlag12tert1 = 1;
		ELSE IF 5.11943 < rainfall_avg_lag12_past_1mo THEN rainlag12tert2 = 1;
	*Add single categorical variable for rain quartiles;
	IF rainlag12cat0 = 1 THEN rainlag12cat = 0;
		ELSE IF rainlag12cat1 = 1 THEN rainlag12cat = 1;
		ELSE IF rainlag12cat2 = 1 THEN rainlag12cat = 2;
		ELSE IF rainlag12catover = 1 THEN rainlag12cat = 3;
		ELSE rainlag12cat = .;
	*Add single categorical variables for rain tertiles;
	IF rainlag12tert0 = 1 THEN rainlag12tert = 0;
		ELSE IF rainlag12tert1 = 1 THEN rainlag12tert = 1;
		ELSE IF rainlag12tert2 = 1 THEN rainlag12tert = 2;
		ELSE rainlag12tert = .;
	*Code rainfall anomaly in quadratic and cubic terms;
	rainanm_lag6_1mo_vs_3mo_cubic = rainanm_lag6_1mo_vs_3mo**3;
	rainanm_lag6_1mo_vs_3mo_sq = rainanm_lag6_1mo_vs_3mo**2;
	*add restricted quadratic splines at quartile breakpoints (-1.3770700, 0.0817440, 1.9877110);
	rainanm1v3_lag6_q1 = rainanm_lag6_1mo_vs_3mo + 1.37707;
		IF rainanm1v3_lag6_q1 < 0  THEN rainanm1v3_lag6_q1 = 0;
	rainanm1v3_lag6_q2 = rainanm_lag6_1mo_vs_3mo - 0.0817440;
		IF rainanm1v3_lag6_q2 < 0 THEN rainanm1v3_lag6_q2 = 0;
	rainanm1v3_lag6_q3 = rainanm_lag6_1mo_vs_3mo - 1.9877110;
		IF rainanm1v3_lag6_q3 < 0 THEN rainanm1v3_lag6_q3 = 0;
	rainanm1v3_lag6_sq1 = rainanm1v3_lag6_q1**2;
	rainanm1v3_lag6_sq2 = rainanm1v3_lag6_q2**2;
	rainanm1v3_lag6_sq3 = rainanm1v3_lag6_q3**2;
	rainanm1v3_lag6_spline1 = rainanm1v3_lag6_sq1 - rainanm1v3_lag6_sq3;
	rainanm1v3_lag6_spline2 = rainanm1v3_lag6_sq2 - rainanm1v3_lag6_sq3;
*Code year of study as binary indicator variables;
year18 = 0;
year19 = 0;
year20 = 0;
year21 = 0;
year22 = 0;
IF year_nolag = 2018 THEN year18 = 1;
	ELSE IF year_nolag = 2019 THEN year19 = 1;
	ELSE IF year_nolag = 2020 THEN year20 = 1;
	ELSE IF year_nolag = 2021 THEN year21 = 1;
	ELSE IF year_nolag = 2022 THEN year22 = 1;
*Add binary variable representing those screened in short wet season (October-December);
IF month_nolag = 10 | month_nolag = 11 | month_nolag = 12 THEN wetshort = 1;
ELSE wetshort = 0;
IF month_nolag = 3 | month_nolag = 4 | month_nolag = 5 THEN wetlong = 1;
ELSE wetlong = 0;
*create single categorical season variables for bivariate analyses (NOT for regressions);
IF wetshort = 1 THEN seasoncat = 2;
ELSE IF wetlong THEN seasoncat = 1;
ELSE seasoncat = 0;
*include interaction terms for short and long wet season and Pf positivity;
IF pf = 1 & wetlong = 1 THEN pf_wetlong = 1;
ELSE pf_wetlong = 0;
IF pf = 1 & wetshort = 1 THEN pf_wetshort = 1;
ELSE pf_wetshort = 0;
*Add binary variable representing those screened in (6w-lagged)short wet season (October-December);
IF month_lag6 = 10 | month_lag6 = 11 | month_lag6 = 12 THEN wetshort_lag6 = 1;
ELSE wetshort_lag6 = 0;
IF month_lag6 = 3 | month_lag6 = 4 | month_lag6 = 5 THEN wetlong_lag6 = 1;
ELSE wetlong_lag6 = 0;
*create single categorical season variables for bivariate analyses (NOT for regressions);
IF wetshort_lag6 = 1 THEN seasoncat_lag6 = 2;
ELSE IF wetlong_lag6 THEN seasoncat_lag6 = 1;
ELSE seasoncat_lag6 = 0;
*include interaction terms for short and long wet season and Pf positivity;
IF pf = 1 & wetlong_lag6 = 1 THEN pf_wetlong_lag6 = 1;
ELSE pf_wetlong_lag6 = 0;
IF pf = 1 & wetshort_lag6 = 1 THEN pf_wetshort_lag6 = 1;
ELSE pf_wetshort_lag6 = 0;
*add binary variable representing those screened in wet season (either) or dry season;
*dry = 0, wet = 1;
IF wetlong_lag6 = 1 | wetshort_lag6 = 1 THEN wet_either_lag6 = 1;
ELSE wet_either_lag6 = 0;
*Code po18S Ct value tertile;
IF 0 < po18s_ct THEN DO;
	IF po18s_ct <= 39.58 THEN po18s_tert0 = 1;
		ELSE po18s_tert0 = 0;
	IF (39.58 < po18s_ct & po18s_ct <= 41.66) THEN po18s_tert1 = 1;
		ELSE po18s_tert1=0;
	IF po18s_ct > 41.66 THEN po18s_tert2 = 1;
		ELSE po18s_tert2 = 0;
END;
*Code single categorical variable for tabulating all levels of po18S ct;
IF po18s_tert0 = 1 THEN po18s_tert = 0;
	ELSE IF po18s_tert1 = 1 THEN po18s_tert = 1;
	ELSE IF po18s_tert2 = 1 THEN po18s_tert = 2;
*Code pf18S Ct value tertile;
IF 0 < pf18s_ct THEN DO;
	IF pf18s_ct <= 32.4 THEN pf18s_tert0 = 1;
	ELSE pf18s_tert0 = 0;
	IF (32.4 < pf18s_ct & pf18s_ct <= 37.41) THEN pf18s_tert1 = 1;
		ELSE pf18s_tert1=0;
	IF pf18s_ct > 37.41 THEN pf18s_tert2 = 1;
		ELSE pf18s_tert2 = 0;
END;
*Code single categorical variable for tabulating all levels of pf18S ct;
IF pf18s_tert0 = 1 THEN pf18s_tert = 0;
	ELSE IF pf18s_tert1 = 1 THEN pf18s_tert = 1;
	ELSE IF pf18s_tert2 = 1 THEN pf18s_tert = 2;
*format and label Ct values;
FORMAT po18s_tert po18s_tert_.;
LABEL po18s_tert='Po18S Ct/Density';
FORMAT pf18s_tert pf18s_tert_.;
LABEL pf18s_tert='Pf18S Ct/Density';
*Reclassify cross-reactive Pf-Po double positives;
IF qpcr_pfdens_screening > 100 & po18s_ct > 42 THEN DO;
		po = 0;
		po18s_ct = .;
		qpcr_po_speciation___1 = 0;
		qpcr_po_speciation___2 = 0;
		qpcr_po_speciation___3 = 0;
		qpcr_po_speciation_curtisi = .;
		qpcr_po_speciation_wallikeri = .;
		END;
RUN;

*confirm no cross-reactive samples kept in dataset;
PROC FREQ DATA=transmit_cleaned;
	WHERE qpcr_pfdens_screening > 100 & po18s_ct > 42;
	TABLES pf*po;
RUN;

*Assess distirbution of binary variables and confirm correct recoding of other variables;
PROC FREQ DATA=transmit_cleaned;
	TABLES pf18s_tert*pf;
	TABLES po18s_ct*po_40;
	TABLES pf*po;
	*proportion positive for pf;
	TABLES pf / BINOMIAL (LEVEL = 2);
	*proportion positive for po;
	TABLES po / BINOMIAL (LEVEL = 2);
	*proportion female;
	TABLES enrollment_gender / BINOMIAL (LEVEL = 2);
	*proportion positive in North of Bagamoyo;
	TABLES region / BINOMIAL (LEVEL = 1);
	TABLES month_nolag*wetlong;
	TABLES pf*pf_wetlong;
	TABLES pf*pf_wetshort;
	TABLES month_lag6*seasoncat_lag6;
RUN;

*identify any samples with low Poc Ct and high Pow Ct to ensure are called as Poc mono-infections;
PROC PRINT DATA=transmit_cleaned;
	WHERE (qpcr_po_speciation_curtisi < 39 & qpcr_po_speciation_curtisi > 0 & qpcr_po_speciation_wallikeri >= 44);
RUN;

*identify any samples with low Pow Ct and high Poc Ct to ensure are called as Pow mono-infections;
PROC PRINT DATA=transmit_cleaned;
	WHERE (qpcr_po_speciation_wallikeri < 39 & qpcr_po_speciation_wallikeri > 0 & qpcr_po_speciation_curtisi >= 44);
RUN;

*identify any samples with erroneously inflated po18s Ct values (only run to 45 cycles);
PROC PRINT DATA=transmit_cleaned;
	WHERE (po18s_ct >= 45);
RUN;

*determine final counts and distributions for key variables;
*And examine marginal associations of specific variables;
PROC FREQ DATA=transmit_cleaned;
	TABLES pf*po;
	TABLES po*enrollment_gender / NOROW CHISQ;
	TABLES pf*enrollment_gender / NOROW CHISQ;
	TABLES po*agecat / NOROW CHISQ TREND;
	TABLES pf*agecat / NOROW CHISQ TREND;
	TABLES po*region / NOROW CHISQ;
	TABLES pf*region / NOROW CHISQ;
	TABLES po*seasonlag6 / NOROW CHISQ;
	TABLES pf*seasonlag6 / NOROW CHISQ;
	TABLES po*year_nolag /NOROW TREND CHISQ;
	TABLES pf*year_nolag / NOROW TREND CHISQ;
	TABLES po18s_tert*enrollment_gender / NOROW TREND;
	TABLES pf18s_tert*enrollment_gender / NOROW TREND;
	TABLES po18s_tert*pf / NOROW TREND;
	TABLES po18s_tert*agecat / NOROW TREND;
	TABLES pf18s_tert*year_nolag / CHISQ NOROW;
	TABLES po18s_tert*year_nolag / CHISQ NOROW;
	TABLES pf_submicroscopic;
	TABLES po_submicroscopic;
RUN;

*Calculate power to detect prevalence ratios observed in Pf for Po;
*First, by sex;
PROC POWER;
	TWOSAMPLEFREQ test=PCHI
	REFPROPORTION = 0.115
	RR=1.2520
	GROUPNS= 4853 | 2320
	POWER = .;
RUN;

*Next, between age groups;
*Important note: SAS' PROC POWER does not enable power calculations of three proportions;
*So we will calculate power to detect enrichment of adolescents (ages 12-15yo) vs. children (5-11yo);
PROC POWER;
	TWOSAMPLEFREQ test=PCHI
	REFPROPORTION = 0.123
	RR=1.07395
	GROUPNS= 3743 | 3429
	POWER = .;
RUN;

*Determine distributions of age and average rainfall over last 1mo with no lag;
PROC MEANS DATA=transmit_cleaned P25 P50 P75 MAX MIN MEAN STDDEV NMISS;
	VAR rainfall_avg_lag0_past_1mo;
	VAR enrollment_age;
RUN;

*Determine distribution of po18S rRNA qPCR values among those Po-positive;
PROC MEANS DATA=transmit_cleaned P25 P50 P75 MAX MIN MEAN STDDEV NMISS;
	WHERE po=1;
	VAR po18s_ct;
	VAR po_parasite_density;
RUN;

*Determine distribution of pf18S rRNA qPCR values among those Pf-positive;
PROC MEANS DATA=transmit_cleaned P25 P50 P75 MAX MIN MEAN STDDEV NMISS;
	WHERE pf=1;
	VAR pf18s_ct;
	VAR qpcr_pfdens_screening;
RUN;

*Prepare output dataset for creating histograms of parasite density by Pf;
DATA transmit_parasite_density;
	SET transmit_cleaned;
	KEEP pf po qpcr_pfdens_screening po18s_ct;
RUN;

PROC FREQ DATA=transmit_cleaned;
	TABLES pf*qpcr_pfdens_screening / MISSING;
RUN;

PROC EXPORT 
  DATA=transmit_parasite_density
  dbms=csv
  OUTFILE="\\Mac\Home\Documents\P. ovale Epi Bagamoyo\transmit_parasite_density.csv" 
  replace;
RUN;

*Perform wilcoxon rank sum of po18s and pf18s ct values;
DATA po18s_cts (rename=(po18s_ct=ct));
	SET transmit_cleaned;
	assay = "po18s";
	KEEP po18s_ct assay;
RUN;

DATA pf18s_cts (rename=(pf18s_ct=ct));
	SET transmit_cleaned;
	assay = "pf18s";
	KEEP pf18s_ct assay;
RUN;

DATA qpcr_ct_comparison;
	SET po18s_cts pf18s_cts;
RUN;

PROC PRINT DATA=pf18s_cts;
	WHERE ct >40.99;
RUN;

PROC NPAR1WAY DATA=qpcr_ct_comparison wilcoxon;
    class assay;
    var ct;
run;


******Graph LOESS lots of Pf and Po prevalence by various demographic variables;
ODS GRAPHICS / LOESSMAXOBS=70000;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=raincat X=rainanm_lag0_1mo_vs_3mo / MARKERATTRS=(size=0) SMOOTH=0.25;
	YAXIS LABELATTRS=(Size=14) VALUEATTRS=(size=12);
	XAXIS LABELATTRS=(Size=14) VALUEATTRS=(size=12);
RUN;

ODS GRAPHICS / LOESSMAXOBS=70000;
PROC SGPLOT DATA=transmit_cleaned;
	SCATTER Y=raincat X=rainanm_lag0_1mo_vs_3mo / MARKERATTRS=(size=1);
	YAXIS LABELATTRS=(Size=14) VALUEATTRS=(size=12);
	XAXIS LABELATTRS=(Size=14) VALUEATTRS=(size=12);
RUN;

ODS LISTING GPATH="\\Mac\Home\Documents\P. ovale Epi Bagamoyo\Figures";
ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="rainavg1mo-lag0_date_loess" imagefmt=jpeg height=6in width=8in noborder;
ODS GRAPHICS / LOESSMAXOBS=70000;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=rainfall_avg_lag0_past_1mo X=sasdate / MARKERATTRS=(size=1) SMOOTH = 0.05;
	YAXIS LABELATTRS=(Size=14) VALUEATTRS=(size=12);
	XAXIS LABELATTRS=(Size=14) VALUEATTRS=(size=12);
RUN;
ODS GRAPHICS OFF;

ODS LISTING GPATH="\\Mac\Home\Documents\P. ovale Epi Bagamoyo\Figures";
ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="rainanm1v3-lag0_date_loess" imagefmt=jpeg height=6in width=8in noborder;
ODS GRAPHICS / LOESSMAXOBS=70000;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=rainanm_lag0_1mo_vs_3mo X=sasdate / MARKERATTRS=(size=1) SMOOTH = 0.05;
	YAXIS LABELATTRS=(Size=14) VALUEATTRS=(size=12);
	XAXIS LABELATTRS=(Size=14) VALUEATTRS=(size=12);
RUN;
ODS GRAPHICS OFF;

****Investigte functional form of continuous variables like age and rainfall by loess plots;
ODS GRAPHICS / LOESSMAXOBS=70000;
ODS LISTING GPATH="\\Mac\Home\Documents\P. ovale Epi Bagamoyo\Figures";
goptions GSFMODE =REPLACE;

*po-date of screening;
ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po_date_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=po X=sasdate / MARKERATTRS=(size=0) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.3 by 0.10) LABELATTRS=(Size=14) VALUEATTRS=(size=12);
	XAXIS LABELATTRS=(Size=14) VALUEATTRS=(size=12);
RUN;
ODS GRAPHICS OFF;

*po-age at screening;
ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po_age_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=po X=enrollment_age / MARKERATTRS=(size=0) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.3 by 0.10) LABELATTRS=(size=14) VALUEATTRS=(size=12) LABEL="P. ovale prevalence";
	XAXIS VALUES = (0 to 100 by 10) LABELATTRS=(size=14) VALUEATTRS=(size=12) LABEL = "Participant age";
RUN;
ODS GRAPHICS OFF;

*po-avg rainfall over past 1mo;
ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po_rainfall-1mo-lag0_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=po X=rainfall_avg_lag0_past_1mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
RUN;
ODS GRAPHICS OFF;

*po-avg rainfall over past 1mo (6w lag);
ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po_rainfall-1mo-lag6_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=po X=rainfall_avg_lag6_past_1mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
RUN;
ODS GRAPHICS OFF;

*po-avg rainfall over past 1mo (9w lag);
ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po_rainfall-1mo-lag9_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=po X=rainfall_avg_lag9_past_1mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
RUN;
ODS GRAPHICS OFF;

*po-avg rainfall over past 1mo (12w lag);
ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po_rainfall-1mo-lag12_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=po X=rainfall_avg_lag12_past_1mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
RUN;
ODS GRAPHICS OFF;

*po-avg rainfall over past 3mo;
ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po_rainfall-3mo-lag0_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=po X=rainfall_avg_lag0_past_3mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
RUN;
ODS GRAPHICS OFF;

*po-avg rainfall over past 3mo (6w lag);
ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po_rainfall-3mo-lag6_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=po X=rainfall_avg_lag6_past_3mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
RUN;
ODS GRAPHICS OFF;

*po-avg rainfall over past 3mo (9w lag);
ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po_rainfall-3mo-lag9_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=po X=rainfall_avg_lag9_past_3mo/ MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
RUN;
ODS GRAPHICS OFF;

*po-avg rainfall over past 3mo (12w lag);
ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po_rainfall-3mo-lag12_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=po X=rainfall_avg_lag12_past_3mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
RUN;
ODS GRAPHICS OFF;

*po-avg rainfall over past 6mo;
ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po_rainfall-6mo-lag0_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=po X=rainfall_avg_lag0_past_6mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
RUN;
ODS GRAPHICS OFF;

*po-avg rainfall over past 6mo (6w lag);
ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po_rainfall-6mo-lag6_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=po X=rainfall_avg_lag6_past_6mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
RUN;
ODS GRAPHICS OFF;

*po-avg rainfall over past 6mo (9w lag);
ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po_rainfall-6mo-lag9_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=po X=rainfall_avg_lag9_past_6mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
RUN;
ODS GRAPHICS OFF;

*po-avg rainfall over past 6mo (12w lag);
ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po_rainfall-6mo-lag12_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=po X=rainfall_avg_lag12_past_6mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
RUN;
ODS GRAPHICS OFF;


*Check whether difference in rainfall between various timeframes correlates with po positivity;
ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po_rainanm1v3-lag0_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=po X=rainanm_lag0_1mo_vs_3mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po_rainanm1v3-lag6_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=po X=rainanm_lag6_1mo_vs_3mo / MARKERATTRS=(size=0) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.2 by 0.10) LABELATTRS=(size=14);
	XAXIS LABELATTRS=(size=14);
RUN;
ODS GRAPHICS OFF;

*po- difference in average rainfall between last 1mo and last 3mo (6w lag) during wet season;
ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po_rainanm1v3-lag6_wet-only_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=po X=rainanm_lag6_1mo_vs_3mo / MARKERATTRS=(size=0) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.2 by 0.10) LABELATTRS=(size=14);
	XAXIS LABELATTRS=(size=14);
	WHERE seasonlag6 = 0;
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po_rainanm1v3-lag6_dry-only_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=po X=rainanm_lag6_1mo_vs_3mo / MARKERATTRS=(size=0) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.4 by 0.10) LABELATTRS=(size=14);
	XAXIS LABELATTRS=(size=14);
	WHERE seasonlag6 = 1;
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po_rainanm1v3-lag9_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=po X=rainanm_lag9_1mo_vs_3mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
	XAXIS VALUES = (-7 to 7 by 1);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po_rainanm1v3-lag12_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=po X=rainanm_lag12_1mo_vs_3mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
	XAXIS VALUES = (-7 to 7 by 1);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po_rainanm1v6-lag0_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=po X=rainanm_lag0_1mo_vs_6mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po_rainanm1v6-lag6_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=po X=rainanm_lag6_1mo_vs_6mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
	XAXIS VALUES = (-7 to 7 by 1);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po_rainanm1v6-lag9_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=po X=rainanm_lag9_1mo_vs_6mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
	XAXIS VALUES = (-7 to 7 by 1);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po_rainanm1v6-lag12_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=po X=rainanm_lag12_1mo_vs_6mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
	XAXIS VALUES = (-7 to 7 by 1);
RUN;
ODS GRAPHICS OFF;


*Investigate associations between Pf positivity and other variables;
ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf_date_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=pf X=sasdate / MARKERATTRS=(size=0) SMOOTH=0.25;
	YAXIS LABELATTRS=(Size=14) VALUEATTRS=(size=12);
	XAXIS LABELATTRS=(Size=14) VALUEATTRS=(size=12);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf_age_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=pf X=enrollment_age / MARKERATTRS=(size=0) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10) LABELATTRS=(size=14) VALUEATTRS=(size=12) LABEL="P. falciparum prevalence";
	XAXIS VALUES = (0 to 100 by 10) LABELATTRS=(size=14) VALUEATTRS=(size=12) LABEL="Participant age";
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf_rainfall-1mo-lag0_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=pf X=rainfall_avg_lag0_past_1mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf_rainfall-1mo-lag6_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=pf X=rainfall_avg_lag6_past_1mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf_rainfall-1mo-lag9_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=pf X=rainfall_avg_lag9_past_1mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf_rainfall-1mo-lag12_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=pf X=rainfall_avg_lag12_past_1mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf_rainfall-3mo-lag0_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=pf X=rainfall_avg_lag0_past_3mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf_rainfall-3mo-lag6_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=pf X=rainfall_avg_lag6_past_3mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf_rainfall-3mo-lag9_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=pf X=rainfall_avg_lag9_past_3mo/ MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf_rainfall-3mo-lag12_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=pf X=rainfall_avg_lag12_past_3mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf_rainfall-6mo-lag0_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=pf X=rainfall_avg_lag0_past_6mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
	*XAXIS VALUES = (0 to 300 by 10);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf_rainfall-6mo-lag6_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=pf X=rainfall_avg_lag6_past_6mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
	*XAXIS VALUES = (0 to 300 by 10);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf_rainfall-6mo-lag9_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=pf X=rainfall_avg_lag9_past_6mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
	*XAXIS VALUES = (0 to 300 by 10);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf_rainfall-6mo-lag12_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=pf X=rainfall_avg_lag12_past_6mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
	*XAXIS VALUES = (0 to 300 by 10);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf_rainanm1v3-lag0_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=pf X=rainanm_lag0_1mo_vs_3mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
	*XAXIS VALUES = (0 to 300 by 10);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf_rainanm1v3-lag6_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=pf X=rainanm_lag6_1mo_vs_3mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10) LABELATTRS=(size=14) VALUEATTRS=(size=12);
	XAXIS LABELATTRS=(size=14) VALUEATTRS=(size=12);
RUN;
ODS graphics OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf_rainanm1v3-lag6_wet-only_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=pf X=rainanm_lag6_1mo_vs_3mo / MARKERATTRS=(size=0) SMOOTH=0.25;
	YAXIS LABELATTRS=(size=14);
	XAXIS LABELATTRS=(size=14);
	WHERE seasonlag6 = 0;
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf_rainanm1v3-lag6_dry-only_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=pf X=rainanm_lag6_1mo_vs_3mo / MARKERATTRS=(size=0) SMOOTH=0.25;
	YAXIS LABELATTRS=(size=14);
	XAXIS LABELATTRS=(size=14);
	WHERE seasonlag6 = 1;
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf_rainanm1v3-lag9_loess" imagefmt=jpeg height=6in width=8in noborder;;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=pf X=rainanm_lag9_1mo_vs_3mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
	XAXIS VALUES = (-7 to 7 by 1);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf_rainanm1v3-lag12_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=pf X=rainanm_lag12_1mo_vs_3mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
	XAXIS VALUES = (-7 to 7 by 1);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf_rainanm1v6-lag0_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=pf X=rainanm_lag0_1mo_vs_6mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
	*XAXIS VALUES = (0 to 300 by 10);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf_rainanm1v6-lag6_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=pf X=rainanm_lag6_1mo_vs_6mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
	XAXIS VALUES = (-7 to 7 by 1);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf_rainanm1v6-lag9_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=pf X=rainanm_lag9_1mo_vs_6mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
	XAXIS VALUES = (-7 to 7 by 1);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf_rainanm1v6-lag12_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=pf X=rainanm_lag12_1mo_vs_6mo / MARKERATTRS=(size=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10);
	XAXIS VALUES = (-7 to 7 by 1);
RUN;
ODS GRAPHICS OFF;

*Model Po and Pf positivity by enrollment sex;
PROC GENMOD DATA= transmit_cleaned;
	MODEL po = enrollment_gender / LINK=LOG;
	ESTIMATE "RR of F vs. M" enrollment_gender 1 / EXP;
	OUTPUT OUT=sexfig PRED = sexpred L = sexl U = sexu;
RUN;

PROC FREQ DATA=transmit_cleaned;
	TABLES po*enrollment_gender;
RUN;

PROC GENMOD DATA= transmit_cleaned;
	MODEL pf = enrollment_gender / LINK=LOG;
	ESTIMATE "RR of F vs. M" enrollment_gender 1 / EXP;
	OUTPUT OUT=sexfig PRED = sexpred L = sexl U = sexu;
RUN;

*Model Po parasite density by sex and season;
PROC GENMOD DATA= transmit_cleaned;
	MODEL po18s_tert0 = enrollment_gender / LINK=LOG;
RUN;

PROC GENMOD DATA= transmit_cleaned;
	MODEL pf18s_tert0 = enrollment_gender / LINK=LOG;
RUN;

PROC GENMOD DATA= transmit_cleaned;
	MODEL po18s_tert2 = enrollment_gender / LINK=LOG;
RUN;

PROC GENMOD DATA= transmit_cleaned;
	MODEL pf18s_tert2 = enrollment_gender / LINK=LOG;
RUN;

PROC GENMOD DATA= transmit_cleaned;
	MODEL po18s_tert0 = seasonlag6 / LINK=LOG;
RUN;

PROC GENMOD DATA= transmit_cleaned;
	MODEL pf18s_tert0 = seasonlag6 / LINK=LOG;
RUN;

PROC GENMOD DATA= transmit_cleaned;
	MODEL po18s_tert2 = seasonlag6 / LINK=LOG;
RUN;

PROC GENMOD DATA= transmit_cleaned;
	MODEL pf18s_tert2 = seasonlag6 / LINK=LOG;
RUN;

*Calculate Pearson's Chi-Square test for parasite density vs other variables;
PROC FREQ DATA=transmit_cleaned;
	TABLES po18s_tert*seasonlag6 / CHISQ;
	TABLES pf18s_tert*seasonlag6 / CHISQ;
	TABLES po18s_tert*enrollment_gender / CHISQ;
	TABLES pf18s_tert*enrollment_gender / CHISQ;
RUN;

*****Model Po positivity by various functional forms of age;
*linear age;
PROC GENMOD DATA= transmit_cleaned;
	MODEL po = enrollment_age / LINK=LOG;
	ESTIMATE "RR of a given age compared to one year younger" enrollment_age 1 / EXP;
	ESTIMATE "RR of 18 vs 12" enrollment_age 6 / EXP;
	ESTIMATE "RR of 30 vs 18" enrollment_age 12 / EXP;
	OUTPUT OUT=linearfig PRED = linearpred L = linearl U = linearu;
RUN;

*squared age;
PROC GENMOD DATA= transmit_cleaned;
	MODEL po = sqage enrollment_age / LINK=LOG;
	ESTIMATE "RR of 18 vs 12" enrollment_age 6 sqage 180 / EXP;
	ESTIMATE "RR of 30 vs 18" enrollment_age 12 sqage 576/ EXP;
	OUTPUT OUT=quadfig PRED = quadpred L = quadl U = quadu;
RUN;

*cubed age;
PROC GENMOD DATA= transmit_cleaned;
	MODEL po = enrollment_age sqage cubage / LINK=LOG;
	ESTIMATE "RR of 18 vs 12" enrollment_age 6 sqage 180 cubage 4104/ EXP;
	ESTIMATE "RR of 30 vs 18" enrollment_age 12 sqage 576 cubage 21168/ EXP;
	OUTPUT OUT=cubfig PRED = cubpred L = cubl U = cubu;
RUN;

*restricted quadratic splines at quartiles;
PROC GENMOD DATA = transmit_cleaned;
MODEL po = enrollment_age rqspline11 rqspline18 / LINK=LOG;
	ESTIMATE "RR of 18 vs 12" enrollment_age 6 rqspline11 48 / EXP;
	ESTIMATE "RR of 30 vs 18" enrollment_age 12 rqspline11 312 rqspline18 144/ EXP;
	OUTPUT OUT=splinefig PRED = splinepred L = splinel U = splineu;
RUN;

*restricted quadratic splines at literature breakpoints (11, 15yo);
PROC GENMOD DATA = transmit_cleaned;
MODEL po = enrollment_age rqspline10/ LINK=LOG;
	OUTPUT OUT=litsplinefig PRED = splinepred L = splinel U = splineu;
RUN;

*age categorical at quartiles;
PROC GENMOD DATA = transmit_cleaned;
MODEL po = agecat18 agecat30 agecatover / LINK=LOG;
	ESTIMATE "RR of 18 vs 12" agecat18 1 / EXP;
	ESTIMATE "RR of 30 vs 18" agecat30 1 agecat18 -1 / EXP;
	OUTPUT OUT=cat4fig PRED = catpred L = catl U = catu;
RUN;

*age categorical at literature breakpoints (11, 15yo);
PROC GENMOD DATA = transmit_cleaned;
MODEL po = agecat15 agecat16up / LINK=LOG;
	ESTIMATE "RR of 14 vs 9" agecat15 1 /EXP;
	ESTIMATE "RR of 30 vs 14" agecat16up 1 agecat15 -1 /EXP;
	OUTPUT OUT=cat3fig PRED = catpred L = catl U = catu;
RUN;

****Model Pf positivity by various functional forms of age;
PROC GENMOD DATA= transmit_cleaned;
	MODEL pf = enrollment_age / LINK=LOG;
	ESTIMATE "RR of a given age compared to one year younger" enrollment_age 1 / EXP;
	ESTIMATE "RR of 18 vs 12" enrollment_age 6 / EXP;
	ESTIMATE "RR of 30 vs 18" enrollment_age 12 / EXP;
	OUTPUT OUT=linearfig PRED = linearpred L = linearl U = linearu;
RUN;

PROC GENMOD DATA= transmit_cleaned;
	MODEL pf = enrollment_age sqage / LINK=LOG;
	ESTIMATE "RR of 18 vs 12" enrollment_age 6 sqage 180 / EXP;
	ESTIMATE "RR of 30 vs 18" enrollment_age 12 sqage 576/ EXP;
	OUTPUT OUT=quadfig PRED = quadpred L = quadl U = quadu;
RUN;

PROC GENMOD DATA= transmit_cleaned;
	MODEL pf = enrollment_age sqage cubage / LINK=LOG;
	ESTIMATE "RR of 18 vs 12" enrollment_age 6 sqage 180 cubage 4104/ EXP;
	ESTIMATE "RR of 30 vs 18" enrollment_age 12 sqage 576 cubage 21168/ EXP;
	OUTPUT OUT=cubfig PRED = cubpred L = cubl U = cubu;
RUN;

PROC GENMOD DATA = transmit_cleaned;
	MODEL pf = enrollment_age rqspline11 rqspline18 / LINK=LOG;
	ESTIMATE "RR of 18 vs 12" enrollment_age 6 rqspline11 48 / EXP;
	ESTIMATE "RR of 30 vs 18" enrollment_age 12 rqspline11 312 rqspline18 144/ EXP;
	OUTPUT OUT=splinefig PRED = splinepred L = splinel U = splineu;
RUN;

PROC GENMOD DATA = transmit_cleaned;
MODEL pf = enrollment_age rqspline10/ LINK=LOG;
	OUTPUT OUT=litsplinefig PRED = splinepred L = splinel U = splineu;
RUN;

PROC GENMOD DATA = transmit_cleaned;
MODEL pf = agecat18 agecat30 agecatover / LINK=LOG;
	ESTIMATE "RR of 18 vs 12" agecat18 1 / EXP;
	ESTIMATE "RR of 30 vs 18" agecat30 1 agecat18 -1 / EXP;
	OUTPUT OUT=cat4fig PRED = catpred L = catl U = catu;
RUN;

PROC GENMOD DATA = transmit_cleaned;
MODEL pf = agecat15 agecat16up / LINK=LOG;
	ESTIMATE "RR of 14 vs 9" agecat15 1 /EXP;
	ESTIMATE "RR of 30 vs 14" agecat16up 1 agecat15 -1 /EXP;
	OUTPUT OUT=cat3fig PRED = catpred L = catl U = catu;
RUN;

*****Model Po and Pf positivity by functional forms of rainfall;
PROC GENMOD DATA = transmit_cleaned;
MODEL po = rainanm_lag6_1mo_vs_3mo / LINK=LOG;
	ESTIMATE "RR by an increase in 1mm/day" rainanm_lag6_1mo_vs_3mo 1 /EXP;
	OUTPUT OUT=rain_po_anm1v3_lag6fig PRED = poanmpred L = poanml U = poanmu;
RUN;

PROC GENMOD DATA = transmit_cleaned;
MODEL pf = rainanm_lag6_1mo_vs_3mo / LINK=LOG;
	ESTIMATE "RR by an increase in 1mm/day" rainanm_lag6_1mo_vs_3mo 1 /EXP;
	OUTPUT OUT=rain_pf_anm1v3_lag6fig PRED = pfanmpred L = pfanml U = pfanmu;
RUN;

PROC GENMOD DATA = transmit_cleaned;
MODEL po = rainanm_lag9_1mo_vs_3mo / LINK=LOG;
	ESTIMATE "RR by an increase in 1mm/day" rainanm_lag9_1mo_vs_3mo 1 /EXP;
	OUTPUT OUT=rain_po_anm1v3_lag9fig PRED = poanmpred L = poanml U = poanmu;
RUN;

PROC GENMOD DATA = transmit_cleaned;
MODEL pf = rainanm_lag9_1mo_vs_3mo / LINK=LOG;
	ESTIMATE "RR by an increase in 1mm/day" rainanm_lag9_1mo_vs_3mo 1 /EXP;
	OUTPUT OUT=rain_pf_anm1v3_lag9fig PRED = poanmpred L = poanml U = poanmu;
RUN;

*Model Malaria prevalence by quadratic and cubic version of rain_anm_lag6_1mo_v_3mo variable;
PROC GENMOD DATA = transmit_cleaned;
MODEL po = rainanm_lag6_1mo_vs_3mo rainanm_lag6_1mo_vs_3mo_sq / LINK=LOG;
	OUTPUT OUT=rain_po_anm1v3_sq_lag6fig PRED = poanmpred L = poanml U = poanmu;
RUN;

PROC GENMOD DATA = transmit_cleaned;
MODEL po = rainanm_lag6_1mo_vs_3mo rainanm_lag6_1mo_vs_3mo_sq rainanm_lag6_1mo_vs_3mo_cubic / LINK=LOG;
RUN;

PROC GENMOD DATA = transmit_cleaned;
MODEL pf = rainanm_lag6_1mo_vs_3mo rainanm_lag6_1mo_vs_3mo_sq / LINK=LOG;
	OUTPUT OUT=rain_pf_anm1v3_sq_lag6fig PRED = pfanmpred L = pfanml U = pfanmu;
RUN;


PROC GENMOD DATA = transmit_cleaned;
MODEL pf = rainanm_lag6_1mo_vs_3mo rainanm_lag6_1mo_vs_3mo_sq rainanm_lag6_1mo_vs_3mo_cubic/ LINK=LOG;
RUN;

*Model restricted quadratic splines of rain anomaly variable defined at quartile breakpoints;
PROC GENMOD DATA = transmit_cleaned;
MODEL po = rainanm_lag6_1mo_vs_3mo rainanm1v3_lag6_spline1 rainanm1v3_lag6_spline2 / LINK=LOG;
	OUTPUT OUT=rain_po_anm1v3_lag6_rqsplinefig PRED = poanmpred L = poanml U = poanmu;
RUN;

PROC GENMOD DATA = transmit_cleaned;
MODEL pf = rainanm_lag6_1mo_vs_3mo rainanm1v3_lag6_spline1 rainanm1v3_lag6_spline2 / LINK=LOG;
	OUTPUT OUT=rain_pf_anm1v3_lag6_rqsplinefig PRED = pfanmpred L = pfanml U = pfanmu;
RUN;

*merge rain-anomaly-6wlag for pf and po models for graphing on same plot;
PROC SORT DATA=rain_po_anm1v3_lag6_rqsplinefig;
	BY screening_id;
RUN;

PROC SORT DATA=rain_pf_anm1v3_lag6_rqsplinefig;
	BY screening_id;
RUN;

DATA rain_anm1v3_lag6_rqsplinefig;
	MERGE rain_po_anm1v3_lag6_rqsplinefig (in = in1) rain_pf_anm1v3_lag6_rqsplinefig (keep = screening_id pfanmpred pfanml pfanmu);
	BY screening_id;
	IF in1 = 1;
RUN;

*sort all prediction datasets by age so they plot cleanly when comparing LOESS to model;
PROC SORT DATA=linearfig;
	BY enrollment_age;
RUN;
PROC SORT DATA=quadfig;
	BY enrollment_age;
RUN;

PROC SORT DATA=cubfig;
	BY enrollment_age;
RUN;

PROC SORT DATA=splinefig;
	BY enrollment_age;
RUN;

PROC SORT DATA=cat4fig;
	BY enrollment_age;
RUN;

PROC SORT DATA=cat3fig;
	BY enrollment_age;
RUN;

PROC SORT DATA=rain_po_anm1v3_lag6fig;
	BY rainanm_lag6_1mo_vs_3mo;
RUN;

PROC SORT DATA=rain_pf_anm1v3_lag6fig;
	BY rainanm_lag6_1mo_vs_3mo;
RUN;

PROC SORT DATA=rain_po_anm1v3_sq_lag6fig;
	BY rainanm_lag6_1mo_vs_3mo;
RUN;

PROC SORT DATA=rain_pf_anm1v3_sq_lag6fig;
	BY rainanm_lag6_1mo_vs_3mo;
RUN;

PROC SORT DATA=rain_po_anm1v3_lag6_rqsplinefig;
	BY rainanm_lag6_1mo_vs_3mo;
RUN;

PROC SORT DATA=rain_anm1v3_lag6_rqsplinefig;
	BY rainanm_lag6_1mo_vs_3mo;
RUN;

PROC SORT DATA=rain_pf_anm1v3_lag6_rqsplinefig;
	BY rainanm_lag6_1mo_vs_3mo;
RUN;

***Plot LOESS and regression-modeled relationships between Po/Pf positivity and continuous variables;
/* Must use this line in order for subsquent graphs to be saved to local files*/
ODS LISTING GPATH="\\Mac\Home\Documents\P. ovale Epi Bagamoyo\Figures";


ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po_rainanm1v3-lag6_linear" imagefmt=jpeg height=6in width=6in noborder;
PROC SGPLOT DATA=rain_po_anm1v3_lag6fig;
	SERIES x=rainanm_lag6_1mo_vs_3mo y=poanmpred / markers markerattrs=(symbol="circlefilled" size=0 color="black") lineattrs=(thickness=1 color="blue") legendlabel="Linear Model";
	SERIES x=rainanm_lag6_1mo_vs_3mo y=poanml / lineattrs=(pattern=2 color="red" thickness=1) legendlabel="Lower 95% CI";
	SERIES x=rainanm_lag6_1mo_vs_3mo y=poanmu / lineattrs=(pattern=2 color="red" thickness=1) legendlabel="Upper 95% CI"; 
	LOESS Y=po X=rainanm_lag6_1mo_vs_3mo / MARKERATTRS=(size=0) lineattrs=(pattern=1 color="black" thickness=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.2 by 0.10) LABEL = "Predicted P. ovale prevalence" LABELATTRS=(size=14);
	XAXIS LABELATTRS=(size=14);
	TITLE "Linear model and Loess curve of P. ovale prevalence by 6w-lagged difference in average daily rainfall between previous month and previous 3 months";
RUN;
TITLE;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf_rainanm1v3-lag6_linear" imagefmt=jpeg height=6in width=6in noborder;
PROC SGPLOT DATA=rain_pf_anm1v3_lag6fig;
	SERIES x=rainanm_lag6_1mo_vs_3mo y=pfanmpred / markers markerattrs=(symbol="circlefilled" size=0 color="black") lineattrs=(thickness=1 color="blue") legendlabel="Linear Model";
	SERIES x=rainanm_lag6_1mo_vs_3mo y=pfanml / lineattrs=(pattern=2 color="red" thickness=1) legendlabel="Lower 95% CI";
	SERIES x=rainanm_lag6_1mo_vs_3mo y=pfanmu / lineattrs=(pattern=2 color="red" thickness=1) legendlabel="Upper 95% CI"; 
	LOESS Y=pf X=rainanm_lag6_1mo_vs_3mo / MARKERATTRS=(size=0) lineattrs=(pattern=1 color="black" thickness=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10) LABEL = "Predicted P. falciparum prevalence" LABELATTRS=(size=14);
	XAXIS LABELATTRS=(size=14);
	TITLE "Linear model and Loess curve of P. falciparum prevalence by 6w-lagged difference in average daily rainfall between previous month and previous 3 months";
RUN;
TITLE;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po_rainanm1v3-lag6_sq" imagefmt=jpeg height=6in width=6in noborder;
PROC SGPLOT DATA=rain_po_anm1v3_sq_lag6fig;
	SERIES x=rainanm_lag6_1mo_vs_3mo y=poanmpred / markers markerattrs=(symbol="circlefilled" size=0 color="black") lineattrs=(thickness=1 color="blue") legendlabel="Linear Model";
	SERIES x=rainanm_lag6_1mo_vs_3mo y=poanml / lineattrs=(pattern=2 color="red" thickness=1) legendlabel="Lower 95% CI";
	SERIES x=rainanm_lag6_1mo_vs_3mo y=poanmu / lineattrs=(pattern=2 color="red" thickness=1) legendlabel="Upper 95% CI"; 
	LOESS Y=po X=rainanm_lag6_1mo_vs_3mo / MARKERATTRS=(size=0) lineattrs=(pattern=1 color="black" thickness=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.2 by 0.10) LABEL = "Predicted P. ovale prevalence" LABELATTRS=(size=14);
	XAXIS LABELATTRS=(size=14);
	TITLE "Quadratic model and Loess curve of P. ovale prevalence by 6w-lagged difference in average daily rainfall between previous month and previous 3 months";
RUN;
TITLE;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf_rainanm1v3-lag6_sq" imagefmt=jpeg height=6in width=6in noborder;
PROC SGPLOT DATA=rain_pf_anm1v3_sq_lag6fig;
	SERIES x=rainanm_lag6_1mo_vs_3mo y=pfanmpred / markers markerattrs=(symbol="circlefilled" size=0 color="black") lineattrs=(thickness=1 color="blue") legendlabel="Linear Model";
	SERIES x=rainanm_lag6_1mo_vs_3mo y=pfanml / lineattrs=(pattern=2 color="red" thickness=1) legendlabel="Lower 95% CI";
	SERIES x=rainanm_lag6_1mo_vs_3mo y=pfanmu / lineattrs=(pattern=2 color="red" thickness=1) legendlabel="Upper 95% CI"; 
	LOESS Y=pf X=rainanm_lag6_1mo_vs_3mo / MARKERATTRS=(size=0) lineattrs=(pattern=1 color="black" thickness=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10) LABEL = "Predicted P. falciparum prevalence" LABELATTRS=(size=14);
	XAXIS LABELATTRS=(size=14);
	TITLE "Quadratic model and Loess curve of P. falciparum prevalence by 6w-lagged difference in average daily rainfall between previous month and previous 3 months";
RUN;
TITLE;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po_rainanm1v3-lag6_rqspline" imagefmt=jpeg height=6in width=6in noborder;
PROC SGPLOT DATA=rain_po_anm1v3_lag6_rqsplinefig;
	SERIES x=rainanm_lag6_1mo_vs_3mo y=poanmpred / markers markerattrs=(symbol="circlefilled" size=0 color="black") lineattrs=(thickness=1 color="blue") legendlabel="Linear Model";
	SERIES x=rainanm_lag6_1mo_vs_3mo y=poanml / lineattrs=(pattern=2 color="red" thickness=1) legendlabel="Lower 95% CI";
	SERIES x=rainanm_lag6_1mo_vs_3mo y=poanmu / lineattrs=(pattern=2 color="red" thickness=1) legendlabel="Upper 95% CI"; 
	LOESS Y=po X=rainanm_lag6_1mo_vs_3mo / MARKERATTRS=(size=0) lineattrs=(pattern=1 color="black" thickness=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.2 by 0.10) LABEL = "Predicted P. ovale prevalence" LABELATTRS=(size=14);
	XAXIS LABELATTRS=(size=14);
	TITLE "Restricted quadratic spline model and Loess curve of P. ovale prevalence by 6w-lagged difference in average daily rainfall between previous month and previous 3 months";
RUN;
TITLE;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf_rainanm1v3-lag6_rqspline" imagefmt=jpeg height=6in width=6in noborder;
PROC SGPLOT DATA=rain_pf_anm1v3_lag6_rqsplinefig;
	SERIES x=rainanm_lag6_1mo_vs_3mo y=pfanmpred / markers markerattrs=(symbol="circlefilled" size=0 color="black") lineattrs=(thickness=1 color="blue") legendlabel="Linear Model";
	SERIES x=rainanm_lag6_1mo_vs_3mo y=pfanml / lineattrs=(pattern=2 color="red" thickness=1) legendlabel="Lower 95% CI";
	SERIES x=rainanm_lag6_1mo_vs_3mo y=pfanmu / lineattrs=(pattern=2 color="red" thickness=1) legendlabel="Upper 95% CI"; 
	LOESS Y=pf X=rainanm_lag6_1mo_vs_3mo / MARKERATTRS=(size=0) lineattrs=(pattern=1 color="black" thickness=1) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.5 by 0.10) LABEL = "Predicted P. falciparum prevalence" LABELATTRS=(size=14);
	XAXIS LABELATTRS=(size=14);
	TITLE "Restricted quadratic spline model and Loess curve of P. falciparum prevalence by 6w-lagged difference in average daily rainfall between previous month and previous 3 months";
RUN;
TITLE;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pfpo_rainanm1v3-lag6_rqspline" imagefmt=jpeg height=6in width=6in noborder;
PROC SGPLOT DATA=rain_anm1v3_lag6_rqsplinefig;
	SERIES x=rainanm_lag6_1mo_vs_3mo y=pfanmpred / markers markerattrs=(symbol="circlefilled" size=0 color="blue") lineattrs=(pattern=2 thickness=1 color="blue") legendlabel="Pf spline model";
	SERIES x=rainanm_lag6_1mo_vs_3mo y=pfanml / lineattrs=(pattern=3 color="blue" thickness=1) legendlabel="Lower 95% CI";
	SERIES x=rainanm_lag6_1mo_vs_3mo y=pfanmu / lineattrs=(pattern=3 color="blue" thickness=1) legendlabel="Upper 95% CI"; 
	LOESS Y=pf X=rainanm_lag6_1mo_vs_3mo / MARKERATTRS=(size=0) lineattrs=(pattern=1 color="blue" thickness=1) SMOOTH=0.25 LEGENDLABEL="Pf LOESS curve";
	SERIES x=rainanm_lag6_1mo_vs_3mo y=poanmpred / markers markerattrs=(symbol="circlefilled" size=0 color="orange") lineattrs=(pattern=2 thickness=1 color="orange") legendlabel="Po spline model";
	SERIES x=rainanm_lag6_1mo_vs_3mo y=poanml / lineattrs=(pattern=3 color="orange" thickness=1) legendlabel="Lower 95% CI";
	SERIES x=rainanm_lag6_1mo_vs_3mo y=poanmu / lineattrs=(pattern=3 color="orange" thickness=1) legendlabel="Upper 95% CI"; 
	LOESS Y=po X=rainanm_lag6_1mo_vs_3mo / MARKERATTRS=(size=0) lineattrs=(pattern=1 color="orange" thickness=1) SMOOTH=0.25 LEGENDLABEL="Po LOESS curve";
	YAXIS VALUES = (0 to 0.5 by 0.10) LABEL = "Prevalence" LABELATTRS=(size=16) VALUEATTRS=(size=14);
	XAXIS LABELATTRS=(size=16) LABEL="Change in rainfall" VALUEATTRS=(size=14);
	*TITLE "Restricted quadratic spline model and Loess curve of P. falciparum prevalence by 6w-lagged difference in average daily rainfall between previous month and previous 3 months";
	KEYLEGEND / VALUEATTRS=(SIZE=13) LOCATION=INSIDE;
RUN;
TITLE;
ODS GRAPHICS OFF;

*******Perform univariate analysis of Po positivity by covariates;

*Bivariate distributions of Po positivity, Pf positivity, mono-infections;
PROC FREQ DATA=transmit_cleaned;
	TABLES enrollment_gender*po / NOPERCENT NOCOL;
	TABLES region*po / NOPERCENT NOCOL;
	TABLES catvillage*po / NOPERCENT NOCOL;
	TABLES enrollment_gender*pf/ NOPERCENT NOCOL;
	TABLES region*pf / NOPERCENT NOCOL;
	TABLES catvillage*pf / NOPERCENT NOCOL OUT=village_pfprev OUTPCT;
	TABLES seasoncat_lag6;
	TABLES seasoncat_lag6*po / NOPERCENT NOCOL CHISQ;
	TABLES seasoncat_lag6*pf / NOPERCENT NOCOL CHISQ;
	TABLES seasoncat_lag6*po*pf / NOPERCENT;
	TABLES region*po_mono / NOPERCENT NOCOL;
	TABLES region*pf_mono / NOPERCENT NOCOL;
	TABLES seasoncat_lag6*po_mono / NOPERCENT NOCOL CHISQ;
	TABLES seasoncat_lag6*pf_mono / NOPERCENT NOCOL CHISQ;
	TABLES enrollment_gender*po_mono / NOPERCENT NOCOL CHISQ;
	TABLES enrollment_gender*pf_mono / NOPERCENT NOCOL CHISQ;
	TABLES agecat*po_mono / NOPERCENT NOCOL CHISQ;
	TABLES agecat*pf_mono / NOPERCENT NOCOL CHISQ;
RUN;


*evaluate bivariate associations within seasons;
*Dry season associations;
PROC FREQ DATA=transmit_cleaned;
	WHERE seasoncat_lag6 = 0;
	TABLES region*po / NOPERCENT NOCOL CHISQ;
	TABLES region*pf / NOPERCENT NOCOL CHISQ;
	*slight Po increase in children (both groups) and decrease in adults, not sig;
	TABLES agecat*po / NOPERCENT NOCOL CHISQ;
	TABLES agecat*pf / NOPERCENT NOCOL CHISQ;
	TABLES enrollment_gender*po / NOPERCENT NOCOL CHISQ;
	TABLES enrollment_gender*pf / NOPERCENT NOCOL CHISQ;
RUN;

*Long wet season associations;
PROC FREQ DATA=transmit_cleaned;
	WHERE seasoncat_lag6 = 1;
	TABLES region*po / NOPERCENT NOCOL CHISQ;
	TABLES region*pf / NOPERCENT NOCOL CHISQ;
	*strong Po enrichment in adults;
	TABLES agecat*po / NOPERCENT NOCOL CHISQ;
	TABLES agecat*pf / NOPERCENT NOCOL CHISQ;
	*increased detection in females;
	TABLES enrollment_gender*po / NOPERCENT NOCOL CHISQ;
	TABLES enrollment_gender*pf / NOPERCENT NOCOL CHISQ;
RUN;

*Short wet season associations;
PROC FREQ DATA=transmit_cleaned;
	WHERE seasoncat_lag6 = 2;
	TABLES region*po / NOPERCENT NOCOL CHISQ;
	TABLES region*pf / NOPERCENT NOCOL CHISQ;
	TABLES agecat*po / NOPERCENT NOCOL CHISQ;
	TABLES agecat*pf / NOPERCENT NOCOL CHISQ;
	TABLES enrollment_gender*po / NOPERCENT NOCOL CHISQ;
	TABLES enrollment_gender*pf / NOPERCENT NOCOL CHISQ;
RUN;

*merge village-wide Pf prevalence into broader dataset;
*village will be missing from published dataset;
*alter output of pf prevalence table dataset so it only contains positive prevalence for Pf;
DATA village_pfprev;
	SET village_pfprev;
	WHERE pf = 1;
	KEEP catvillage PCT_ROW;
	RENAME PCT_ROW=village_pf_prev;
RUN;

PROC SORT DATA=village_pfprev;
	BY catvillage;
RUN;

PROC SORT DATA=transmit_cleaned;
	BY catvillage;
RUN;

DATA transmit_cleaned;
	MERGE transmit_cleaned (in = in1) village_pfprev;
	BY catvillage;
	IF in1;
	LABEL village_pf_prev="Village-level Pf Prevalence";
RUN;

*Determine whether Po density was predicted by Pf village prevalence;
ODS LISTING GPATH="\\Mac\Home\Documents\P. ovale Epi Bagamoyo\Figures";
ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="po-density_pf-village-prev_loess" imagefmt=jpeg height=6in width=6in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=po18s_tert X=village_pf_prev / MARKERATTRS=(size=1) SMOOTH=0.25;
RUN;
ODS GRAPHICS OFF;

*Estimate RR between males and females;
PROC GENMOD DATA= transmit_cleaned;
	MODEL po = enrollment_gender / LINK=LOG;
	ESTIMATE "RR of Females vs. Males" enrollment_gender 1 / EXP;
RUN;

*Estimate region Po+ RR after adjusting for categorical age and sex;
PROC GENMOD DATA= transmit_cleaned;
	MODEL po =  agecat15 agecat16up enrollment_gender region / LINK=LOG;
	ESTIMATE "RR of South vs North" region 1 / EXP;
RUN;

*Estimate region Pf+ RR after adjusting for categorical age and sex;
PROC GENMOD DATA= transmit_cleaned;
	MODEL pf = region agecat15 agecat16up enrollment_gender / LINK=LOG;
	ESTIMATE "RR of South vs North" region 1 / EXP;
RUN;	

*Estimate whether po18S Ct was predictive of Pf co-infection among Po+ samples;
PROC GENMOD DATA=transmit_cleaned;
	MODEL pf = po18s_tert1 po18s_tert2 / LINK= LOG;
RUN;

*determine whether Pf positivity is predictive of Po positivity after adjusting for other relevant variables;
*adjusted covariates: rainfall (restricted quadratic splines of rainfall anomaly w 6w lag), age categories at 11 and 15yo breakpoints, year of study, season with 6w lag);
PROC GENMOD DATA= transmit_cleaned;
	MODEL po = pf agecat15 agecat16up wetshort_lag6 wetlong_lag6 rainanm_lag6_1mo_vs_3mo rainanm1v3_lag6_spline1 rainanm1v3_lag6_spline2 year18 year19 year20 year21/ LINK=LOG;
	ESTIMATE "RR of Pf+ vs. Pf-" pf 1 / EXP;
	ESTIMATE "PR of 11-15yo vs. 5-10yo" agecat15 1 / EXP;
	ESTIMATE "PR of 16+yo vs. 5-10yo" agecat16up 1 / EXP;
	ESTIMATE "PR of short wet season vs. dry season" wetshort_lag6 1 / EXP;
	ESTIMATE "PR of long wet season vs. dry season" wetlong_lag6 1 / EXP;
RUN;

*determine whether Pf positivity is predictive of Po positivity when stratifying by dry, long wet, and short wet season, adjusting for other covariates, by including interaction terms;
PROC GENMOD DATA= transmit_cleaned;
	MODEL po = pf wetshort_lag6 pf_wetshort_lag6 wetlong_lag6 pf_wetlong_lag6 rainanm_lag6_1mo_vs_3mo rainanm1v3_lag6_spline1 rainanm1v3_lag6_spline2 agecat15 agecat16up year18 year19 year20 year21/ LINK=LOG;
	ESTIMATE "PR of Pf+ vs. Pf- in dry season" pf 1 / EXP;
	ESTIMATE "PR of Pf+ vs. Pf- in short wet season" pf 1 pf_wetshort_lag6 1 / EXP;
	ESTIMATE "PR of Pf+ vs. Pf- in long wet season" pf 1 pf_wetlong_lag6 1 / EXP;
RUN;

*determine whether Pf positivity is predictive of Po positivity when stratifying by dry, long wet, and short wet season, adjusting for other covariates, by including interaction terms;
*add sex to determine necessity for adjustment;
PROC GENMOD DATA= transmit_cleaned;
	MODEL po = pf enrollment_gender wetshort_lag6 pf_wetshort_lag6 wetlong_lag6 pf_wetlong_lag6 rainanm_lag6_1mo_vs_3mo rainanm1v3_lag6_spline1 rainanm1v3_lag6_spline2 agecat15 agecat16up year18 year19 year20 year21/ LINK=LOG;
	ESTIMATE "PR of Pf+ vs. Pf- in dry season" pf 1 / EXP;
	ESTIMATE "PR of Pf+ vs. Pf- in short wet season" pf 1 pf_wetshort_lag6 1 / EXP;
	ESTIMATE "PR of Pf+ vs. Pf- in long wet season" pf 1 pf_wetlong_lag6 1 / EXP;
RUN;

*calculate interaction between Pf and Po among age-strata, adjusting for other covariates, by including interaction terms;
PROC GENMOD DATA= transmit_cleaned;
	MODEL po = pf agecat15 agecat15pf agecat16up agecat16uppf rainanm_lag6_1mo_vs_3mo rainanm1v3_lag6_spline1 rainanm1v3_lag6_spline2 year18 year19 year20 year21/ LINK=LOG;
	ESTIMATE "RR of Pf+ vs. Pf- among kids under 11" pf 1 / EXP;
	ESTIMATE "RR of Pf+ vs. Pf- among kids 11-15" pf 1 agecat15pf 1/ EXP;
	ESTIMATE "RR of Pf+ vs. Pf- among people over 15yo" pf 1 agecat16uppf 1 / EXP;
RUN;

*determine whether Pf positivity is predictive of Po positivity using conservative Po Ct cutoff of 40 cycles (rather than 45;
*adjusted covariates: rainfall (restricted quadratic splines of rainfall anomaly w 6w lag), age categories at 11 and 15yo breakpoints, year of study);
PROC GENMOD DATA= transmit_cleaned;
	MODEL po_40 = pf agecat15 agecat16up wetshort_lag6 wetlong_lag6 rainanm_lag6_1mo_vs_3mo rainanm1v3_lag6_spline1 rainanm1v3_lag6_spline2 year18 year19 year20 year21/ LINK=LOG;
	ESTIMATE "RR of Pf+ vs. Pf-" pf 1 / EXP;
RUN;

*investigate association of parasite density and age;
ODS GRAPHICS / LOESSMAXOBS=70000;
ODS LISTING GPATH="\\Mac\Home\Documents\P. ovale Epi Bagamoyo\Figures";
goptions GSFMODE =REPLACE;
ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pfdens_age_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=qpcr_pfdens_screening X=enrollment_age / MARKERATTRS=(size=0) SMOOTH=0.25;
	YAXIS LABEL = 'Pf parasite density (parasites/uL)' VALUES=(0 to 10000 by 3000) VALUEATTRS=(size=12) LABELATTRS=(Size=14);
	XAXIS LABELATTRS=(Size=14) VALUEATTRS=(size=12) VALUES= (0 to 80 by 10);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="podens_age_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=po18s_ct X=enrollment_age / MARKERATTRS=(size=0) SMOOTH=0.25;
	YAXIS VALUES = (36 to 42 by 1) LABEL = 'po18s Ct' REVERSE VALUEATTRS=(size=12) LABELATTRS=(Size=14);
	XAXIS LABELATTRS=(Size=14) VALUEATTRS=(size=12) VALUES= (0 to 80 by 10);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pfct_age_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=pf18s_ct X=enrollment_age / MARKERATTRS=(size=0) SMOOTH=0.25;
	YAXIS VALUES = (30 to 42 by 1)LABEL = 'pf18s Ct' REVERSE VALUEATTRS=(size=12) LABELATTRS=(Size=14);
	XAXIS LABELATTRS=(Size=14) VALUEATTRS=(size=12) VALUES= (0 to 80 by 10);
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf-po-ct_age_loess" imagefmt=jpeg height=6in width=6in noborder;
PROC SGPLOT DATA=transmit_cleaned;
	LOESS Y=pf18s_ct X=enrollment_age / MARKERATTRS=(size=0) SMOOTH=0.25 legendlabel="pf18S";
	LOESS Y=po18s_ct X=enrollment_age / MARKERATTRS=(size=0) SMOOTH=0.25 legendlabel="po18S";
	YAXIS VALUES = (30 to 42 by 1) LABEL = '18S rRNA cyclic threshold' REVERSE VALUEATTRS=(size=14) LABELATTRS=(Size=16);
	XAXIS LABELATTRS=(Size=16) LABEL = 'Participant age' VALUEATTRS=(size=14) VALUES= (0 to 85 by 10);
	KEYLEGEND / VALUEATTRS=(SIZE=14) LOCATION=INSIDE;
RUN;
ODS GRAPHICS OFF;



***Perform indirect standardization of transmit age-specific rates to overall bagamoyo population;

*determine stratum-specific counts and number of individuals per stratum;
PROC FREQ DATA=transmit_cleaned;
	TABLES age5bin*po;
RUN;

*Load in data of stratum specific event counts (# positive) and number of individuals per stratum;
data transmit_age_noyoung;
	input age $10. total_pos total_persons;
	*data from TranSMIT;
	datalines;
	05-09 106 946
	10-14 258 2317
	15-19 53 556
	20-24 100 742
	25-29 79 688
	30-34 74 607
	35-39 50 412
	40-44 40 324
	45-49 21 172
	50-54 18 134
	55-59 6 69
	60-64 16 109
	65-69 3 43
	70-74 7 32
	75-79 0 11
	80plus 0 10
	;

*Load in data of bagamoyo age strata. Events coded as zero but will not be used in standardization;
*derived from 2022 census data;
data bagamoyo_age_noyoung;
	input age $10. total_pos total_persons;
	*data is from census, 0 events are arbitrarily placed here;
	datalines;
	05-09 0 18379
	10-14 0 16065
	15-19 0 13579
	20-24 0 12632
	25-29 0 11558
	30-34 0 10135
	35-39 0 8571
	40-44 0 7445
	45-49 0 6036
	50-54 0 4770
	55-59 0 3683
	60-64 0 3387
	65-69 0 2572
	70-74 0 1904
	75-79 0 1209
	80plus 0 1047
	;

*Calculate indirectly standardized estimate for population over the age of 4;
ods graphics on;
PROC STDRATE DATA=bagamoyo_age_noyoung REFDATA=transmit_age_noyoung
METHOD=indirect
STAT=risk
PLOTS=all
;
POPULATION EVENT=total_pos TOTAL=total_persons;
REFERENCE EVENT=total_pos TOTAL=total_persons;
STRATA age / stats smr;
RUN;
ods graphics off;

*Calclate upper bound of age prevalence estimate;
data transmit_age;
	input age $10. total_pos total_persons;
	*count in 0-4 is highest prevalence (raw count made randomly) among other TranSMIT age strata, which is 60-64;
	*remainder of counts are data from TranSMIT;
	datalines;
	00-04 16 113
	05-09 106 946
	10-14 258 2317
	15-19 53 556
	20-24 100 742
	25-29 79 688
	30-34 74 607
	35-39 50 412
	40-44 40 324
	45-49 21 172
	50-54 18 134
	55-59 6 69
	60-64 16 109
	65-69 3 43
	70-74 7 32
	75-79 0 11
	80plus 0 10
	;

*Load in data of bagamoyo age strata. Events coded as zero but will not be used in standardization;
*derived from census data from 2022;
data bagamoyo_age;
	input age $10. total_pos total_persons;
	*data is from census, 0 events are arbitrarily placed here;
	datalines;
	00-04 0 19831
	05-09 0 18379
	10-14 0 16065
	15-19 0 13579
	20-24 0 12632
	25-29 0 11558
	30-34 0 10135
	35-39 0 8571
	40-44 0 7445
	45-49 0 6036
	50-54 0 4770
	55-59 0 3683
	60-64 0 3387
	65-69 0 2572
	70-74 0 1904
	75-79 0 1209
	80plus 0 1047
	;

ods graphics on;
PROC STDRATE DATA=bagamoyo_age REFDATA=transmit_age
METHOD=indirect
STAT=risk
PLOTS=all
;
POPULATION EVENT=total_pos TOTAL=total_persons;
REFERENCE EVENT=total_pos TOTAL=total_persons;
STRATA age / stats smr;
RUN;
ods graphics off;

*Generate lower bound of age prevalence estimates assuming no positivity among childen under 5;
data transmit_age_lowerbound;
	input age $10. total_pos total_persons;
	*count in 0-4 assumes no positivity among young kids;
	*remainder of counts are data from TranSMIT;
	datalines;
	00-04 0 248
	05-09 106 946
	10-14 258 2317
	15-19 53 556
	20-24 100 742
	25-29 79 688
	30-34 74 607
	35-39 50 412
	40-44 40 324
	45-49 21 172
	50-54 18 134
	55-59 6 69
	60-64 16 109
	65-69 3 43
	70-74 7 32
	75-79 0 11
	80plus 0 10
	;


ods graphics on;
PROC STDRATE DATA=bagamoyo_age REFDATA=transmit_age_lowerbound
METHOD=indirect
STAT=risk
PLOTS=all
;
POPULATION EVENT=total_pos TOTAL=total_persons;
REFERENCE EVENT=total_pos TOTAL=total_persons;
STRATA age / stats smr;
RUN;
ods graphics off;

*Generate prevalence estimate by age AND SEX strata;
PROC FREQ DATA=transmit_cleaned;
	TABLES age5bin*enrollment_gender;
	TABLES po*age5bin*enrollment_gender;
RUN;

*Create dataset for prevalence estimate of males and females over 4 years old;
DATA transmit_age_sex_noyoung;
	input age $10. sex $2. total_pos total_persons;
	*data from TranSMIT;
	datalines;
	05-09 F 54 555
	05-09 M 48 391
	10-14 F 134 1180
	10-14 M 124 1137
	15-19 F 27 330
	15-19 M 26 226
	20-24 F 90 644
	20-24 M 10 98
	25-29 F 74 615
	25-29 M 5 73
	30-34 F 63 529
	30-34 M 11 78
	35-39 F 40 347
	35-39 M 10 65
	40-44 F 31 263
	40-44 M 9 61
	45-49 F 15 128
	45-49 M 6 44
	50-54 F 11 82
	50-54 M 7 52
	55-59 F 2 44
	55-59 M 4 25
	60-64 F 12 76
	60-64 M 4 33
	65-69 F 3 28
	65-69 M 0 15
	70-74 F 4 20
	70-74 M 3 12
	75-79 F 0 6
	75-79 M 0 5
	80plus F 0 5
	80plus M 0 5
	;

data bagamoyo_age_sex_noyoung;
	input age $10. sex $2. total_pos total_persons;
	*data is from census, 0 events are arbitrarily placed here;
	datalines;
	05-09 F 0 8982
	05-09 M 0 9397
	10-14 F 0 7776
	10-14 M 0 8289
	15-19 F 0 6542
	15-19 M 0 7037
	20-24 F 0 6260
	20-24 M 0 6372
	25-29 F 0 5687
	25-29 M 0 5871
	30-34 F 0 4893
	30-34 M 0 5242
	35-39 F 0 4016
	35-39 M 0 4555
	40-44 F 0 3540
	40-44 M 0 3905
	45-49 F 0 2733
	45-49 M 0 3303
	50-54 F 0 2292
	50-54 M 0 2478
	55-59 F 0 1785
	55-59 M 0 1898
	60-64 F 0 1839
	60-64 M 0 1548
	65-69 F 0 1500
	65-69 M 0 1072
	70-74 F 0 1088
	70-74 M 0 816
	75-79 F 0 673
	75-79 M 0 536
	80plus F 0 607
	80plus M 0 440
	;

*Calculate standardized prevalence among population over 4;
ods graphics on;
PROC STDRATE DATA=bagamoyo_age_sex_noyoung REFDATA=transmit_age_sex_noyoung
METHOD=indirect
STAT=risk
PLOTS=all
;
POPULATION EVENT=total_pos TOTAL=total_persons;
REFERENCE EVENT=total_pos TOTAL=total_persons;
STRATA age sex/ stats smr;
RUN;
ods graphics off;

DATA transmit_age_sex_upperbound;
	input age $10. sex $2. total_pos total_persons;
	*count in 0-4 is highest TranSMIT stratum-specific prevalence w/ >75 participants, being F 60-64;
	*this count is kept identical between male and female children (original study does not indicate;
	*remainder of counts are data from TranSMIT;
	datalines;
	00-04 F 15789 100000
	00-04 M 15789 100000
	05-09 F 54 555
	05-09 M 48 391
	10-14 F 134 1180
	10-14 M 124 1137
	15-19 F 27 330
	15-19 M 26 226
	20-24 F 90 644
	20-24 M 10 98
	25-29 F 74 615
	25-29 M 5 73
	30-34 F 63 529
	30-34 M 11 78
	35-39 F 40 347
	35-39 M 10 65
	40-44 F 31 263
	40-44 M 9 61
	45-49 F 15 128
	45-49 M 6 44
	50-54 F 11 82
	50-54 M 7 52
	55-59 F 2 44
	55-59 M 4 25
	60-64 F 12 76
	60-64 M 4 33
	65-69 F 3 28
	65-69 M 0 15
	70-74 F 4 20
	70-74 M 3 12
	75-79 F 0 6
	75-79 M 0 5
	80plus F 0 5
	80plus M 0 5
	;

DATA bagamoyo_age_sex;
	INPUT age $10. sex $2. total_pos total_persons;
	*data is from census, 0 events are arbitrarily placed here;
	DATALINES;
	00-04 F 0 9657
	00-04 M 0 10174
	05-09 F 0 8982
	05-09 M 0 9397
	10-14 F 0 7776
	10-14 M 0 8289
	15-19 F 0 6542
	15-19 M 0 7037
	20-24 F 0 6260
	20-24 M 0 6372
	25-29 F 0 5687
	25-29 M 0 5871
	30-34 F 0 4893
	30-34 M 0 5242
	35-39 F 0 4016
	35-39 M 0 4555
	40-44 F 0 3540
	40-44 M 0 3905
	45-49 F 0 2733
	45-49 M 0 3303
	50-54 F 0 2292
	50-54 M 0 2478
	55-59 F 0 1785
	55-59 M 0 1898
	60-64 F 0 1839
	60-64 M 0 1548
	65-69 F 0 1500
	65-69 M 0 1072
	70-74 F 0 1088
	70-74 M 0 816
	75-79 F 0 673
	75-79 M 0 536
	80plus F 0 607
	80plus M 0 440
	;

*Calculate upper bound of standardized prevalence among whole population;
ods graphics on;
PROC STDRATE DATA=bagamoyo_age_sex REFDATA=transmit_age_sex_upperbound
METHOD=indirect
STAT=risk
PLOTS=all
;
POPULATION EVENT=total_pos TOTAL=total_persons;
REFERENCE EVENT=total_pos TOTAL=total_persons;
STRATA age sex/ stats smr;
RUN;
ods graphics off;

*Generate lower bound of prevalence estimate in whole population;
DATA transmit_age_sex_lowerbound;
	input age $10. sex $2. total_pos total_persons;
	*Count in ages 0-4 assumes no positivity among this age group;
	*remainder of counts are data from TranSMIT;
	datalines;
	00-04 F 0 248
	00-04 M 0 248
	05-09 F 54 555
	05-09 M 48 391
	10-14 F 134 1180
	10-14 M 124 1137
	15-19 F 27 330
	15-19 M 26 226
	20-24 F 90 644
	20-24 M 10 98
	25-29 F 74 615
	25-29 M 5 73
	30-34 F 63 529
	30-34 M 11 78
	35-39 F 40 347
	35-39 M 10 65
	40-44 F 31 263
	40-44 M 9 61
	45-49 F 15 128
	45-49 M 6 44
	50-54 F 11 82
	50-54 M 7 52
	55-59 F 2 44
	55-59 M 4 25
	60-64 F 12 76
	60-64 M 4 33
	65-69 F 3 28
	65-69 M 0 15
	70-74 F 4 20
	70-74 M 3 12
	75-79 F 0 6
	75-79 M 0 5
	80plus F 0 5
	80plus M 0 5
	;

*Calculate lower bound of standardized prevalence among whole population;
ods graphics on;
PROC STDRATE DATA=bagamoyo_age_sex REFDATA=transmit_age_sex_lowerbound
METHOD=indirect
STAT=risk
PLOTS=all
;
POPULATION EVENT=total_pos TOTAL=total_persons;
REFERENCE EVENT=total_pos TOTAL=total_persons;
STRATA age sex/ stats smr;
RUN;
ods graphics off;

****calculate transmit wet season prevalence of Pf and Po;

*First, determine distribution of Po positivity across both wet seasons;
PROC FREQ DATA=transmit_cleaned;
	WHERE wet_either_lag6 = 1;
	TABLES age5bin*enrollment_gender;
	TABLES po*age5bin*enrollment_gender;
RUN;

*Enter stratum totals and po cases;
DATA transmit_agesexnoyoung_wet_po;
	input age $10. sex $2. total_pos total_persons;
	*data from TranSMIT;
	datalines;
	05-09 F 30 336
	05-09 M 26 227
	10-14 F 53 576
	10-14 M 63 589
	15-19 F 16 195
	15-19 M 13 104
	20-24 F 58 364
	20-24 M 6 59
	25-29 F 46 346
	25-29 M 3 41
	30-34 F 38 292
	30-34 M 8 50
	35-39 F 27 210
	35-39 M 8 47
	40-44 F 20 149
	40-44 M 5 34
	45-49 F 14 76
	45-49 M 3 31
	50-54 F 6 46
	50-54 M 5 31
	55-59 F 1 27
	55-59 M 3 15
	60-64 F 4 44
	60-64 M 2 17
	65-69 F 2 18
	65-69 M 0 12
	70-74 F 3 13
	70-74 M 2 10
	75-79 F 0 2
	75-79 M 0 5
	80plus F 0 4
	80plus M 0 5
	;

*determine distribution of Pf positivity across both wet seasons;
PROC FREQ DATA=transmit_cleaned;
	WHERE wet_either_lag6 = 1;
	TABLES age5bin*enrollment_gender;
	TABLES pf*age5bin*enrollment_gender;
RUN;

*Enter stratum totals and pf cases;
DATA transmit_agesexnoyoung_wet_pf;
	input age $10. sex $2. total_pos total_persons;
	*data from TranSMIT;
	datalines;
	05-09 F 84 336
	05-09 M 59 227
	10-14 F 159 576
	10-14 M 186 589
	15-19 F 63 195
	15-19 M 49 104
	20-24 F 85 364
	20-24 M 22 59
	25-29 F 85 346
	25-29 M 14 41
	30-34 F 73 292
	30-34 M 18 50
	35-39 F 48 210
	35-39 M 27 47
	40-44 F 40 149
	40-44 M 11 34
	45-49 F 24 76
	45-49 M 9 31
	50-54 F 12 46
	50-54 M 11 31
	55-59 F 4 27
	55-59 M 6 15
	60-64 F 10 44
	60-64 M 4 17
	65-69 F 3 18
	65-69 M 2 12
	70-74 F 3 13
	70-74 M 3 10
	75-79 F 0 2
	75-79 M 1 5
	80plus F 2 4
	80plus M 2 5
	;

*Calculate expected Pf prevalence to bagamoyo population (over age of 5);
ods graphics on;
PROC STDRATE DATA=bagamoyo_age_sex_noyoung REFDATA=transmit_agesexnoyoung_wet_pf
METHOD=indirect
STAT=risk
PLOTS=all
;
POPULATION EVENT=total_pos TOTAL=total_persons;
REFERENCE EVENT=total_pos TOTAL=total_persons;
STRATA age sex/ stats smr;
RUN;
ods graphics off;


*Calculate expected Po prevalence to bagamoyo population (over age of 5);
ods graphics on;
PROC STDRATE DATA=bagamoyo_age_sex_noyoung REFDATA=transmit_agesexnoyoung_wet_po
METHOD=indirect
STAT=risk
PLOTS=all
;
POPULATION EVENT=total_pos TOTAL=total_persons;
REFERENCE EVENT=total_pos TOTAL=total_persons;
STRATA age sex/ stats smr;
RUN;
ods graphics off;

***Analyze the species-specific distributions of P. ovale infections in Bagamoyo;

*rename variables reflecting screening for individual Po species;
DATA transmit_cleaned;
	SET transmit_cleaned;
	RENAME qpcr_po_speciation___1 = species_unknown qpcr_po_speciation___2 = poc qpcr_po_speciation___3 = pow;
RUN;

*Encode variables showing mixed infections and samples which were not tested;
DATA transmit_cleaned;
	SET transmit_cleaned;
	IF (poc = 1 AND pow = 1) THEN mixed = 1;
		ELSE IF (poc = 0 OR pow = 0) THEN mixed = 0;
	IF (poc = 0 and pow = 0 and species_unknown = 0 and po = 1) THEN unspeciated = 1;
		ELSE unspeciated = 0;
RUN;

*create species-status (not for regressions) showing species-identity of all samples;
PROC FORMAT;
	VALUE species_status_ 0='Po-negative' .='no assay run' 1='species unknown' 2='mixed' 3='Poc' 4='Pow';
RUN;

*Encode additional analysis variables for species composition;
DATA transmit_cleaned;
	SET transmit_cleaned;
	*Encode species status compositie variable for tabulations (NOT REGRESSION);
	IF po = 0 THEN species_status = 0;
		ELSE IF unspeciated = 1 THEN species_status = .;
		ELSE IF species_unknown = 1 THEN species_status = 1;
		ELSE IF mixed = 1 THEN species_status = 2;
		ELSE IF poc = 1 THEN species_status = 3;
		ELSE IF pow = 1 THEN species_status = 4;
	LABEL species_status = "Po species identity";
	FORMAT species_status species_status_.;
	*determine conservative species positivity;
	*first, encode conservative positivity as initial positivity;
	pow_conservative = pow;
	poc_conservative = poc;
	mixed_conservative = mixed;
	*then, among mixed infections, if the difference in Ct values exceeds 3 cycles, encode as only positive for the more abundant species;
	IF qpcr_po_speciation_curtisi > . AND qpcr_po_speciation_wallikeri > . THEN DO;
		IF mixed = 1 AND ABS(qpcr_po_speciation_curtisi - qpcr_po_speciation_wallikeri) < 3 THEN mixed_conservative = 1;
			ELSE IF qpcr_po_speciation_curtisi < qpcr_po_speciation_wallikeri THEN DO;
				mixed_conservative = 0;
				pow_conservative = 0;
				poc_conservative = 1;
			END;
			ELSE IF qpcr_po_speciation_wallikeri < qpcr_po_speciation_curtisi THEN DO;
				mixed_conservative = 0;
				pow_conservative = 1;
				poc_conservative = 0;
			END;
	END;
	*encode interaction terms for each species and short/long wet seasons;
	IF pow = 1 & wetlong_lag6 = 1 THEN pow_wetlong_lag6 = 1;
		ELSE pow_wetlong_lag6 = 0;
	IF pow = 1 & wetshort_lag6 = 1 THEN pow_wetshort_lag6 = 1;
		ELSE pow_wetshort_lag6 = 0;
RUN;

*create dataset which doesn't include unspeciated samples so they aren't erroneously treated as negative for Poc and Pow;
DATA transmit_no_missing_po_species;
	SET transmit_cleaned;
	IF unspeciated = 1 THEN DELETE;
RUN;

*Merge in parasite densities (see supplemental materials);
PROC IMPORT OUT=poc_pow_parasite_density
		DATAFILE="\\Mac\Home\Documents\P. ovale Epi Bagamoyo\TranSMIT_Poc-Pow_parasite-density.csv"
		DBMS=csv REPLACE;
	GETNAMES=YES;
RUN;

PROC SORT DATA=poc_pow_parasite_density;
	BY screening_id;
RUN;

DATA poc_pow_parasite_density;
	SET poc_pow_parasite_density;
	KEEP screening_id poc_parasite_density pow_parasite_density;
RUN;

PROC SORT DATA=transmit_cleaned;
	BY screening_id;
RUN;

*merge parasite density;
*recode species-specific parasite density to be 0 if the individual is labelled negative for that species;
*This is necessary for test results that are possible false-positives (see cross-reactivity algorithm and sensitivity analysis);
DATA transmit_cleaned;
	MERGE transmit_cleaned (in = in1) poc_pow_parasite_density;
	BY screening_id;
	IF in1;
	IF poc_parasite_density > 0 & poc = 0 THEN poc_parasite_density = 0;
	IF pow_parasite_density > 0 & pow = 0 THEN pow_parasite_density = 0;
RUN;

*create dataset for exporting species-specific parasite densities;
DATA transmit_parasite_density_final;
	SET transmit_cleaned;
	WHERE poc = 1 | pow = 1;
	KEEP screening_id poc pow mixed poc_parasite_density pow_parasite_density;
RUN;

PROC EXPORT DATA=transmit_parasite_density_final dbms=csv outfile="\\Mac\Home\Documents\P. ovale Epi Bagamoyo\transmit_parasite_density_final.csv" REPLACE;
RUN;


*examine sample distributions by species status;
PROC FREQ DATA=transmit_cleaned;
	TABLES po*poc;
	TABLES po*pow;
	TABLES poc*pow;
	TABLES po*mixed;
	TABLES po*species_unknown;
	TABLES po*species_status / MISSING;
	TABLES pow_conservative*poc_conservative;
RUN;

*assess po18S distributions by Po species composition;
PROC FREQ DATA=transmit_no_missing_po_species;
	WHERE po = 1;
	TABLES species_status*po18s_tert / CHISQ;
RUN;

*Evaluatw Pf co-infection distributions by Poc/Pow and conservative Poc/Pow designation;
PROC FREQ DATA=transmit_no_missing_po_species;
	TABLES pf*poc*pow;
	TABLES pf*poc_conservative*pow_conservative;
RUN;

*Create dataset of only mono-infections to investigate species-specific assocations with demographics;
DATA transmit_pocpow_monoinfections;
	SET transmit_cleaned;
		WHERE poc = 1 | pow = 1;
		IF poc = 1 & pow = 1 THEN DELETE;
RUN;	

*evaluate association with covariates among Poc and Pow mono-infections;
PROC FREQ DATA=transmit_pocpow_monoinfections;
	TABLES species_status*enrollment_gender / CHISQ;
	TABLES species_status*seasoncat_lag6 / CHISQ;
	TABLES species_status*region / CHISQ;
	TABLES species_status*agecat / CHISQ;
RUN;

*evaluate distributions of covariates among Po-positives;
PROC FREQ DATA=transmit_no_missing_po_species;
	WHERE po=1;
	TABLES species_status*enrollment_gender / NOROW NOPERCENT MISSING CHISQ;
	TABLES species_status*season / NOROW NOPERCENT MISSING;
	TABLES species_status*seasonlag6 / NOROW NOPERCENT MISSING;
	TABLES species_status*region / NOROW NOPERCENT CHISQ;
	TABLES species_status*year_nolag / NOPERCENT NOROW MISSING;
RUN;

*Test association of Poc/Pow positivity with demographic variables;
PROC FREQ DATA=transmit_no_missing_po_species;
	TABLES poc*enrollment_gender / NOPERCENT NOROW CHISQ;
	TABLES pow*enrollment_gender / NOPERCENT NOROW CHISQ;
	TABLES poc*seasoncat_lag6 / NOPERCENT NOROW CHISQ;
	TABLES pow*seasoncat_lag6 / NOPERCENT NOROW CHISQ;
	TABLES poc*agecat / NOPERCENT NOROW CHISQ;
	TABLES pow*agecat / NOPERCENT NOROW CHISQ;
	TABLES poc*region / NOPERCENT NOROW CHISQ;
	TABLES pow*region / NOPERCENT NOROW CHISQ;
RUN;

*compare poc and pow prevalence between children and adolescents;
PROC FREQ DATA=transmit_no_missing_po_species;
	WHERE agecat ^= 2;
	TABLES poc*agecat / NOPERCENT NOROW CHISQ;
	TABLES pow*agecat / NOPERCENT NOROW CHISQ;
RUN;

*examine distribution of species positivity by year and by concomitant parasitemia;
*cohort does not include 196 individuals for whom the species-identification assays were not performed;
PROC FREQ DATA=transmit_no_missing_po_species;
	TABLES year_nolag*poc / MISSING NOCOL NOPERCENT;
	TABLES year_nolag*pow / MISSING NOCOL NOPERCENT;
	TABLES poc*pow /MISSING NOROW NOCOL;
	TABLES pf*pow*poc /MISSING;
RUN;

*determine Pf and Po prevalences in wet season by study year in full dataset;
PROC FREQ DATA=transmit_cleaned;
	WHERE wet_either_lag6 = 1;
	TABLES year_nolag*pf / MISSING NOCOL NOPERCENT;
	TABLES year_nolag*po / MISSING NOCOL NOPERCENT;
RUN;

*Determine Poc and Pow prevalences in wet season by study year among samples which could be tested;
PROC FREQ DATA=transmit_no_missing_po_species;
	WHERE wet_either_lag6 = 1;
	TABLES year_nolag*poc / MISSING NOCOL NOPERCENT;
	TABLES year_nolag*pow / MISSING NOCOL NOPERCENT;
RUN;


*model pow prevalence by pf;
*adjusting for categorical age, year of study, season with 6w lag, and restricted quadratic splines of difference in rainfall between prior 1mo and 3mo with 6w lag;
PROC GENMOD DATA=transmit_no_missing_po_species;
	MODEL pow = pf agecat15 agecat16up year18 year19 year20 year21 wetlong_lag6 wetshort_lag6 rainanm_lag6_1mo_vs_3mo rainanm1v3_lag6_spline1 rainanm1v3_lag6_spline2  / LINK=LOG;
	ESTIMATE "RR of Pf+ vs. Pf-" pf 1 / EXP;
RUN;

*estimate pow-prevalence between short, long, and dry seasons;
*sample size too small to decompose by season;
PROC GENMOD DATA=transmit_no_missing_po_species;
	MODEL pow = pf agecat15 agecat16up year18 year19 year20 year21 wetlong_lag6 pf_wetlong_lag6 wetshort_lag6 pf_wetshort_lag6 rainanm_lag6_1mo_vs_3mo rainanm1v3_lag6_spline1 rainanm1v3_lag6_spline2 / LINK=LOG;
	ESTIMATE "RR of Pf+ vs. Pf- in dry season" pf 1 / EXP;
	ESTIMATE "PR of Pf+ vs. Pf- in short wet season" pf 1 pf_wetshort_lag6 1 / EXP;
	ESTIMATE "PR of Pf+ vs. Pf- in long wet season" pf 1 pf_wetlong_lag6 1 / EXP;
RUN;

*model poc prevalence using pf;
*adjusting for categorical age, year of study, season with 6w lag, and restricted quadratic splines of difference in rainfall between prior 1mo and 3mo with 6w lag;
PROC GENMOD DATA=transmit_no_missing_po_species;
	MODEL poc = pf agecat15 agecat16up year18 year19 year20 year21 wetlong_lag6 wetshort_lag6 rainanm_lag6_1mo_vs_3mo rainanm1v3_lag6_spline1 rainanm1v3_lag6_spline2 / LINK=LOG;
	ESTIMATE "RR of Pf+ vs. Pf-" pf 1 / EXP;
RUN;

*estimate poc-prevalence between short, long, and dry seasons;
*sample size too small to decompose by season;
PROC GENMOD DATA=transmit_no_missing_po_species;
	MODEL poc = pf agecat15 agecat16up year18 year19 year20 year21 wetlong_lag6 pf_wetlong_lag6 wetshort_lag6 pf_wetshort_lag6 rainanm_lag6_1mo_vs_3mo rainanm1v3_lag6_spline1 rainanm1v3_lag6_spline2 / LINK=LOG;
	ESTIMATE "RR of Pf+ vs. Pf- in dry season" pf 1 / EXP;
	ESTIMATE "PR of Pf+ vs. Pf- in short wet season" pf 1 pf_wetshort_lag6 1 / EXP;
	ESTIMATE "PR of Pf+ vs. Pf- in long wet season" pf 1 pf_wetlong_lag6 1 / EXP;
RUN;

*model covariate-adjusted PR between Pow- and Pow+;
*adjusting for Pf, categorical age, year of study, season with 6w lag;
PROC GENMOD DATA=transmit_no_missing_po_species;
	MODEL poc = pow pf agecat15 agecat16up year18 year19 year20 year21 wetlong_lag6 wetshort_lag6 / LINK=LOG;
	ESTIMATE "PR of Pow+ vs. Pow-" pow 1 / EXP;
RUN;

*model within seasons using interaction terms;
*adjusting for Pf, categorical age, year of study, season with 6w lag;
PROC GENMOD DATA=transmit_no_missing_po_species;
	MODEL poc = pow pf agecat15 agecat16up year18 year19 year20 year21 wetlong_lag6 pow_wetlong_lag6 wetshort_lag6 pow_wetshort_lag6/ LINK=LOG;
	ESTIMATE "PR of Pow+ vs. Pow- in dry season" pow 1 / EXP;
	ESTIMATE "PR of Pow+ vs. Pow- in long wet season" pow 1 pow_wetlong_lag6 1 / EXP;
	ESTIMATE "PR of Pow+ vs. Pow- in short wet season" pow 1 pow_wetshort_lag6 1 / EXP;
RUN;

*model poc-pow interaction using more conservative cut-offs for possible cross-reactivity;
PROC GENMOD DATA=transmit_no_missing_po_species;
	MODEL poc_conservative = pow_conservative pf agecat15 agecat16up year18 year19 year20 year21 wetlong_lag6 wetshort_lag6 / LINK=LOG;
	ESTIMATE "PR of Pow+ vs. Pow-" pow_conservative 1 / EXP;
RUN;

*model pow prevalence using pf with more conservative cut-offs for cross-reactivity;
PROC GENMOD DATA=transmit_no_missing_po_species;
	MODEL pow_conservative = pf agecat15 agecat16up year18 year19 year20 year21 wetlong_lag6 wetshort_lag6 rainanm_lag6_1mo_vs_3mo rainanm1v3_lag6_spline1 rainanm1v3_lag6_spline2 / LINK=LOG;
	ESTIMATE "PR of Pf+ vs. Pf-" pf 1 / EXP;
RUN;

*model pow prevalence using pf with more conservative cut-offs for cross-reactivity;
PROC GENMOD DATA=transmit_no_missing_po_species;
	MODEL poc_conservative = pf agecat15 agecat16up year18 year19 year20 year21 wetlong_lag6 wetshort_lag6 rainanm_lag6_1mo_vs_3mo rainanm1v3_lag6_spline1 rainanm1v3_lag6_spline2/ LINK=LOG;
	ESTIMATE "PR of Pf+ vs. Pf-" pf 1 / EXP;
RUN;

*Possible biasing from missing data in samples which failed to speciate;
*Mixed infections were enriched in observed data, but this might not be unusual if the samples which failed are all mono-infections;
*Synergism would be pulled to null if all the missing samples were mono-infections, making the observed excess of mixed infections expected;
*perhaps the mixed infections had higher parasite density (which we observed) and were more likely to be successfully speciated;
*Greatest attenuation would occur if our 255 failed samples were actually all mono-infections;
*Specifically, if the 255 missing samples were all mono-infections of each species so that the overall prevalence of each is the same (which maximizes expected mixed infections);
*In this dataset, this would occur with an 83:172 ratio of Poc:Pow mono-infections;
*The following code will randomly assign missing individuals to either status in this ratio;
*If we iterate multiple times, this will show is a general lower bound for the corrected version of this bias if the underlying data is extreme;
*In bivariate analysis, this roughly doubles expected mixed ifnections (7.4 to 17.85);

*First, sort, all observations by merge variable screening id;
PROC SORT DATA=transmit_no_missing_po_species;
	BY screening_id;
RUN;

*Assign missing samples to groups in 83:172 ratio;
*random seed should be changed for each iteration;
PROC SURVEYSELECT DATA=transmit_no_missing_po_species seed=76378 OUT=transmit_missingspecies groups = (83 172);
	WHERE species_unknown = 1;
RUN;

*prep dataset for merge;
DATA transmit_missingspecies;
	SET transmit_missingspecies;
	KEEP screening_id GroupID;
RUN;

*merge missing assignment variable into broader dataset;
DATA transmit_missingspecies;
	MERGE transmit_no_missing_po_species (in = in1) transmit_missingspecies;
	BY screening_id;
	IF in1;
RUN;

*assign new bounding poc and pow variables including the assignment of missing samples;
DATA transmit_missingspecies;
	SET transmit_missingspecies;
	poc_bound = poc;
	pow_bound = pow;
	poc_bound_conservative = poc_conservative;
	pow_bound_conservative = pow_conservative;
	IF GroupID > . THEN DO;
		IF GroupID = 1 THEN poc_bound = 1;
		ELSE IF GroupID = 2 THEN pow_bound = 1;
		END;
	IF GroupID > . THEN DO;
		IF GroupID = 1 THEN poc_bound_conservative = 1;
		ELSE IF GroupID = 2 THEN pow_bound_conservative = 1;
		END;
RUN;

*Calculate synergism in extreme dataset;
*adjusting for Pf, categorical age, year of study, season with 6w lag;
PROC GENMOD DATA=transmit_missingspecies;
	MODEL poc_bound = pow_bound pf agecat15 agecat16up year18 year19 year20 year21 rainanm_lag6_1mo_vs_3mo rainanm1v3_lag6_spline1 rainanm1v3_lag6_spline2 / LINK=LOG;
	ESTIMATE "RR of Pow+ vs. Pow-" pow_bound 1 / EXP;
RUN;

*calculate synergism in extreme dataset with conservative mixed infection assignment;
*adjusting for Pf, categorical age, year of study, season with 6w lag;
PROC GENMOD DATA=transmit_missingspecies;
	MODEL poc_bound_conservative = pow_bound_conservative pf agecat15 agecat16up year18 year19 year20 year21 rainanm_lag6_1mo_vs_3mo rainanm1v3_lag6_spline1 rainanm1v3_lag6_spline2 / LINK=LOG;
	ESTIMATE "RR of Pow+ vs. Pow-" pow_bound_conservative 1 / EXP;
RUN;

*Show randomized frequencies;
PROC FREQ DATA=transmit_missingspecies;
	TABLES pf*poc_bound_conservative*pow_bound_conservative;
RUN;

*Loess plots of Poc and Pow status by various covariates;
ODS LISTING GPATH="\\Mac\Home\Documents\P. ovale Epi Bagamoyo\Figures";
goptions GSFMODE =REPLACE;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="poc_age_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_no_missing_po_species;
	LOESS Y=poc X=enrollment_age / MARKERATTRS=(size=0) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.1 by 0.05) LABELATTRS=(size=14) VALUEATTRS=(size=12) LABEL="P. ovale curtisi prevalence";
	XAXIS VALUES = (0 to 100 by 10) LABELATTRS=(size=14) VALUEATTRS=(size=12) LABEL="Participant age";
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="poc_rainfall-1mo-lag0_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_no_missing_po_species;
	LOESS Y=poc X=rainfall_avg_lag0_past_1mo / MARKERATTRS=(size=0) SMOOTH=0.25;
	YAXIS VALUES=(0 to 0.15 by 0.03) LABELATTRS=(size=14) VALUEATTRS=(size=12);
	XAXIS LABELATTRS=(size=14) VALUEATTRS=(size=12);
RUN;
ODS GRAPHICS OFF;


ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="poc_rainfallanm1v3-lag0_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_no_missing_po_species;
	LOESS Y=poc X=rainanm_lag6_1mo_vs_3mo / MARKERATTRS=(size=0) SMOOTH=0.25;
	YAXIS VALUES=(0 to 0.15 by 0.03) LABELATTRS=(size=14) VALUEATTRS=(size=12);
	XAXIS LABELATTRS=(size=14) VALUEATTRS=(size=12);
RUN;
ODS GRAPHICS OFF;


ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pow_age_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_no_missing_po_species;
	LOESS Y=pow X=enrollment_age / MARKERATTRS=(size=0) SMOOTH=0.25;
	YAXIS VALUES = (0 to 0.1 by 0.05) LABELATTRS=(size=14) VALUEATTRS=(size=12) LABEL="P. ovale wallikeri prevalence";
	XAXIS VALUES = (0 to 100 by 10) LABELATTRS=(size=14) VALUEATTRS=(size=12) LABEL="Participant age";
RUN;
ODS GRAPHICS OFF;

ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pow_rainfall-1mo-lag0_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_no_missing_po_species;
	LOESS Y=pow X=rainfall_avg_lag0_past_1mo / MARKERATTRS=(size=0) SMOOTH=0.25;
	YAXIS VALUES=(0 to 0.15 by 0.03) LABELATTRS=(size=14) VALUEATTRS=(size=12);
	XAXIS LABELATTRS=(size=14) VALUEATTRS=(size=12);
RUN;
ODS GRAPHICS OFF;


ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pow_rainfallanm1v3-lag0_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_no_missing_po_species;
	LOESS Y=pow X=rainanm_lag6_1mo_vs_3mo / MARKERATTRS=(size=0) SMOOTH=0.25;
	YAXIS VALUES=(0 to 0.15 by 0.03) LABELATTRS=(size=14) VALUEATTRS=(size=12);
	XAXIS LABELATTRS=(size=14) VALUEATTRS=(size=12);
RUN;
ODS GRAPHICS OFF;

*composite graph of LOESS association of Pf, Po, Poc, and Pow by participant age;
ODS GRAPHICS /  reset=all LOESSMAXOBS=70000 antialiasmax=68000 imagename="pf-po-poc-pow_age_loess" imagefmt=jpeg height=6in width=8in noborder;
PROC SGPLOT DATA=transmit_no_missing_po_species;
	LOESS Y=pf X=enrollment_age / MARKERATTRS=(size=0) SMOOTH=0.25 LINEATTRS=(COLOR= Blue PATTERN=SOLID) LEGENDLABEL="P. falciparum";
	LOESS Y=po X=enrollment_age / MARKERATTRS=(size=0) SMOOTH=0.25 LINEATTRS=(COLOR= Orange PATTERN=SOLID) LEGENDLABEL="P. ovale";
	LOESS Y=poc X=enrollment_age / MARKERATTRS=(size=0) SMOOTH=0.25 LINEATTRS=(COLOR= Green PATTERN=SOLID) LEGENDLABEL="P. ovale curtisi";
	LOESS Y=pow X=enrollment_age / MARKERATTRS=(size=0) SMOOTH=0.25 LINEATTRS=(COLOR= Purple PATTERN=SOLID) LEGENDLABEL="P. ovale wallikeri";
	YAXIS VALUES = (0 to 0.4 by 0.05) LABELATTRS=(size=16) VALUEATTRS=(size=14) LABEL="Prevalence";
	XAXIS VALUES = (0 to 100 by 10) LABELATTRS=(size=16) VALUEATTRS=(size=14) LABEL="Participant age";
	KEYLEGEND / VALUEATTRS=(SIZE=14) LOCATION=INSIDE;
RUN;
ODS GRAPHICS OFF;

***Prep data for outputting and graphing;

*First, using time of screening for overall prevalences across the study (Pf and Po);
DATA transmit_cleaned;
	SET transmit_cleaned;
	month_nolag_char = PUT(month_nolag, $8.);
	IF LENGTH(STRIP(month_nolag_char)) = 1 THEN month_nolag_char = CAT('0',month_nolag);
	year_nolag_char = PUT(year_nolag,$8.);
	FORMAT study_month $8.;
	study_month = CAT(STRIP(year_nolag_char),"-",STRIP(month_nolag_char));
	pfpo = 0;
	IF (pf = 1 AND po = 1) THEN pfpo = 1;
RUN;

*prep speciated dataset for poc, pow, and mixed prevalences (not including individuals who were not speciated);
DATA transmit_no_missing_po_species;
	SET transmit_no_missing_po_species;
	month_nolag_char = PUT(month_nolag, $8.);
	IF LENGTH(STRIP(month_nolag_char)) = 1 THEN month_nolag_char = CAT('0',month_nolag);
	year_nolag_char = PUT(year_nolag,$8.);
	FORMAT study_month $8.;
	study_month = CAT(STRIP(year_nolag_char),"-",STRIP(month_nolag_char));
	pfpo = 0;
	IF (pf = 1 AND po = 1) THEN pfpo = 1;
RUN;

*calculate study month totals and prevalences of Pf, Po, Pf-Po co-infections;
PROC FREQ DATA=transmit_cleaned;
	TABLES study_month / OUT= transmit_month (rename = (COUNT = N) DROP = PERCENT);
	TABLES study_month*pf / OUTPCT OUT= transmit_pf_month (rename = (COUNT = pf_count PCT_ROW = pf_prev) DROP = PERCENT PCT_COL);
	TABLES study_month*po  / OUTPCT OUT= transmit_po_month (rename = (COUNT = po_count PCT_ROW = po_prev) DROP = PERCENT PCT_COL);
	TABLES study_month*pfpo / OUTPCT OUT= transmit_pfpo_month (rename = (COUNT = pfpo_count PCT_ROW = pfpo_prev) DROP = PERCENT PCT_COL);
RUN;

*calculate study month totals and prevalences of Poc, Pow, Poc-Pow mixed infections;
*Poc, Pow, and mixed infections should only be calculated among samples in which speciation was attemped;
PROC FREQ DATA=transmit_no_missing_po_species;
	TABLES study_month*poc  / OUTPCT OUT= transmit_poc_month (rename = (COUNT = poc_count PCT_ROW = poc_prev) DROP = PERCENT PCT_COL);
	TABLES study_month*pow / OUTPCT OUT= transmit_pow_month (rename = (COUNT = pow_count PCT_ROW = pow_prev) DROP = PERCENT PCT_COL);
	TABLES study_month*mixed / OUTPCT OUT= transmit_mixed_month (rename = (COUNT = mixed_count PCT_ROW = mixed_prev) DROP = PERCENT PCT_COL);
RUN;

*merge prevalences by study month;
DATA transmit_prevalences;
	MERGE transmit_month transmit_pf_month transmit_po_month transmit_pfpo_month transmit_poc_month transmit_pow_month transmit_mixed_month;
	BY study_month;
RUN;

*fix SAS behavior where months without screening are coded as having 100% prevalence. Recode as no observations or positives;
DATA transmit_prevalences;
	SET transmit_prevalences;
	IF po_prev = 100 THEN DO;
		po = 1;
		po_count = 0;
		po_prev = 0;
		END;
	IF pfpo_prev = 100 THEN DO;
		pfpo = 1;
		pfpo_count = 0;
		pfpo_prev = 0;
		END;
	IF poc_prev = 100 THEN DO;
		poc = 1;
		poc_count = 0;
		poc_prev = 0;
		END;
	IF pow_prev = 100 THEN DO;
		pow = 1;
		pow_count = 0;
		pow_prev = 0;
		END;
	IF mixed_prev = 100 THEN DO;
		mixed = 1;
		mixed_count = 0;
		mixed_prev = 0;
		END;
	IF pf = 0 THEN DELETE;
RUN;

*Export dataset for plotting;
proc export 
  data=transmit_prevalences
  dbms=csv 
  outfile="\\Mac\Home\Documents\P. ovale Epi Bagamoyo\transmit_rainfall_malaria-prev_studied-months.csv" 
  replace;
run;


*Prep 6w lag data for outputting and graphing;
***WITH 6w lag for comparison of wet seasons and dry season;
*Participant month categorization reflects the month/season 6w prior to their screening date;
DATA transmit_cleaned;
	SET transmit_cleaned;
	month_lag6_char = PUT(month_lag6, $8.);
	IF LENGTH(STRIP(month_lag6_char)) = 1 THEN month_lag6_char = CAT('0',month_lag6);
	year_lag6_char = PUT(year_lag6,$8.);
	FORMAT study_month_lag6 $8.;
	study_month_lag6 = CAT(STRIP(year_lag6_char),"-",STRIP(month_lag6_char));
	pfpo = 0;
	IF (pf = 1 AND po = 1) THEN pfpo = 1;
RUN;

*prep speciated dataset for poc, pow, and mixed prevalences (not including individuals who were not speciated);
DATA transmit_no_missing_po_species;
	SET transmit_no_missing_po_species;
	month_lag6_char = PUT(month_lag6, $8.);
	IF LENGTH(STRIP(month_lag6_char)) = 1 THEN month_lag6_char = CAT('0',month_lag6);
	year_lag6_char = PUT(year_lag6,$8.);
	FORMAT study_month_lag6 $8.;
	study_month_lag6 = CAT(STRIP(year_lag6_char),"-",STRIP(month_lag6_char));
	pfpo = 0;
	IF (pf = 1 AND po = 1) THEN pfpo = 1;
RUN;

PROC FREQ DATA=transmit_cleaned;
	TABLES study_month_lag6 / OUT= transmit_month_lag6 (rename = (COUNT = N) DROP = PERCENT);
	TABLES study_month_lag6*pf / OUTPCT OUT= transmit_pf_month_lag6 (rename = (COUNT = pf_count PCT_ROW = pf_prev) DROP = PERCENT PCT_COL);
	TABLES study_month_lag6*po  / OUTPCT OUT= transmit_po_month_lag6 (rename = (COUNT = po_count PCT_ROW = po_prev) DROP = PERCENT PCT_COL);
	TABLES study_month_lag6*pfpo / OUTPCT OUT= transmit_pfpo_month_lag6 (rename = (COUNT = pfpo_count PCT_ROW = pfpo_prev) DROP = PERCENT PCT_COL);
RUN;

*Poc, Pow, and mixed infections should only be calculated among samples in which speciation was attemped;
PROC FREQ DATA=transmit_no_missing_po_species;
	TABLES study_month_lag6*poc  / OUTPCT OUT= transmit_poc_month_lag6 (rename = (COUNT = poc_count PCT_ROW = poc_prev) DROP = PERCENT PCT_COL);
	TABLES study_month_lag6*pow / OUTPCT OUT= transmit_pow_month_lag6 (rename = (COUNT = pow_count PCT_ROW = pow_prev) DROP = PERCENT PCT_COL);
	TABLES study_month_lag6*mixed / OUTPCT OUT= transmit_mixed_month_lag6 (rename = (COUNT = mixed_count PCT_ROW = mixed_prev) DROP = PERCENT PCT_COL);
RUN;

DATA transmit_prevalences_lag6;
	MERGE transmit_month_lag6 transmit_pf_month_lag6 transmit_po_month_lag6 transmit_pfpo_month_lag6 transmit_poc_month_lag6 transmit_pow_month_lag6 transmit_mixed_month_lag6;
	BY study_month_lag6;
RUN;

DATA transmit_prevalences_lag6;
	SET transmit_prevalences_lag6;
	IF po_prev = 100 THEN DO;
		po = 1;
		po_count = 0;
		po_prev = 0;
		END;
	IF pfpo_prev = 100 THEN DO;
		pfpo = 1;
		pfpo_count = 0;
		pfpo_prev = 0;
		END;
	IF poc_prev = 100 THEN DO;
		poc = 1;
		poc_count = 0;
		poc_prev = 0;
		END;
	IF pow_prev = 100 THEN DO;
		pow = 1;
		pow_count = 0;
		pow_prev = 0;
		END;
	IF mixed_prev = 100 THEN DO;
		mixed = 1;
		mixed_count = 0;
		mixed_prev = 0;
		END;
	IF pf = 0 THEN DELETE;
RUN;

PROC PRINT DATA=transmit_prevalences_lag6;
RUN;

proc export 
  data=transmit_prevalences_lag6
  dbms=csv 
  outfile="\\Mac\Home\Documents\P. ovale Epi Bagamoyo\transmit_rainfall_malaria-prev_studied-months_lag6.csv" 
  replace;
run;
