DM "log; clear; odsresults; clear;";
proc datasets library=work kill nolist;
quit;
/*****************************************************************************
	Name: Kelly Carey-Ewend
	Program: P._ovale_epi_bagaomoyo_datacleaning.sas
   	Date: 06-05-2024
	Description: SAS program file for cleaning dataset for TranSMIT
*****************************************************************************/
OPTIONS MERGENOBY=warn NODATE NONUMBER FORMCHAR="|----|+|---+=|-/\<>*";
FOOTNOTE "P._ovale_epi_bagaomoyo_datacleaning.sas run at %SYSFUNC(DATETIME(), DATETIME.) by Kelly Carey-Ewend";
/******************************* begin program ******************************/
TITLE;
LIBNAME dir '\\Mac\Home\Documents\P. ovale Epi Bagamoyo';

/* Edit the following line to reflect the full path to your CSV file*/
%let csv_file = '\\Mac\Home\Documents\P. ovale Epi Bagamoyo\TRANSMIT-Screening_KellyCE_01-08-25.csv';
OPTIONS nofmterr;

*format input variables;
proc format;
	value enrollment_gender_ 0='Male' 1='Female';
	value enrollment_antimalarial_yn_ 0='No' 1='Yes' 
		98='Don''t know' 99='Refused';
	value season_ 0='Wet' 1='Dry';
	value region_ 0='North' 1='South';
	run;

*read in input dataset;
*NOTE: if drawing from publication, the dataset will be missing the enrollment_village variable for de-identification;
*adjust the lines of this procedure as needed, or contact the authors for a full dataset;
data work.redcap; %let _EFIERR_ = 0;
infile &csv_file  delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;

	informat screening_id $500. ;
	informat enrollment_id $500. ;
	informat enrollment_age best32. ;
	informat enrollment_gender best32. ;
	informat enrollment_screen_date $500. ;
	informat enrollment_village $500. ;
	informat enrollment_time_clinic time5. ;
	informat enrollment_antimalarial_yn best32. ;
	informat qpcr_reps_pos_screening best32. ;
	informat qpcr_reps_pos_screening_po best32. ;

	format screening_id $500. ;
	format enrollment_id $500. ;
	format enrollment_age best12. ;
	format enrollment_gender best12. ;
	format enrollment_screen_date $500. ;
	format enrollment_village $500. ;
	format enrollment_time_clinic time5. ;
	format enrollment_antimalarial_yn best12. ;
	format qpcr_reps_pos_screening best12. ;
	format qpcr_reps_pos_screening_po best12. ;
	format qpcr_pfdens_screening best12.;
	format enrollment_symptoms_yn best12.;

input
	screening_id $
	enrollment_id $
	enrollment_age
	enrollment_gender
	enrollment_screen_date $
	enrollment_village $
	enrollment_time_clinic
	enrollment_antimalarial_yn
	qpcr_reps_pos_screening
	qpcr_ct_screening
	qpcr_pfdens_screening
	qpcr_reps_pos_screening_po
	qpcr_ct_screening_po
	qpcr_po_speciation___1
	qpcr_po_speciation___2
	qpcr_po_speciation___3
	qpcr_po_speciation_curtisi
	qpcr_po_speciation_wallikeri
	enrollment_symptoms_yn
;
if _ERROR_ then call symput('_EFIERR_',"1");
run;

*view contents of dataset;
PROC CONTENTS DATA=redcap;
RUN;

*Apply labels to variable;
data redcap;
	set redcap;
	label screening_id='Screening ID';
	label enrollment_id='Subject ID';
	label enrollment_age='Participant age:';
	label enrollment_gender='Gender';
	label enrollment_screen_date='Screening date:';
	label enrollment_village='Village name:';
	label enrollment_time_clinic='Time to clinic:';
	label enrollment_antimalarial_yn='Has the participant used an antimalarial in the last 28 days?';
	label qpcr_reps_pos_screening='Number of P. falciparum positive replicates';
	label qpcr_reps_pos_screening_po='Number of P. ovale positive replicates';
	label qpcr_po_speciation___1="Failed to speciate";
	label qpcr_po_speciation___2="P. ovale curtisi";
	label qpcr_po_speciation___3="P. ovale wallikeri";
	label qpcr_po_speciation_curtisi="Poc speciation Ct value";
	label qpcr_po_speciation_wallikeri="Pow speciation Ct value";
	format enrollment_gender enrollment_gender_.;
	format enrollment_antimalarial_yn enrollment_antimalarial_yn_.;
run;

*Create new dataset for processing;
DATA screening;
	SET redcap;
RUN;

*Load in daily CHIRPS rainfall for bagamoyo (with lagged comparisons of average rainfall in past _ months to average rainfall in past _ months;
PROC IMPORT OUT=chirps_daily_anomaly
		DATAFILE="\\Mac\Home\Documents\P. ovale Epi Bagamoyo\CHIRPS_Bagamoyo_Daily_Rainfall.xls"
		DBMS=XLS REPLACE;
	GETNAMES=YES;
RUN;

*Format date in daily rainfall set as sasdate for merging;
DATA chirps_daily_anomaly;
	SET chirps_daily_anomaly;
	FORMAT Date date9.;
	*sasdate = INPUT(Date,MMDDYY10.);
	RENAME Date = sasdate;
	rainanm_lag0_1mo_vs_3mo =  INPUT(rainanm_lag0_1mo_vs_3mo_char,8.);
	rainanm_lag6_1mo_vs_3mo =  INPUT(rainanm_lag6_1mo_vs_3mo_char,8.);
	rainanm_lag9_1mo_vs_3mo =  INPUT(rainanm_lag9_1mo_vs_3mo_char,8.);
	rainanm_lag12_1mo_vs_3mo = INPUT(rainanm_lag12_1mo_vs_3mo_char,8.);
	rainanm_lag0_1mo_vs_6mo =  INPUT(rainanm_lag0_1mo_vs_6mo_char,8.);
	rainanm_lag6_1mo_vs_6mo =  INPUT(rainanm_lag6_1mo_vs_6mo_char,8.);
	rainanm_lag9_1mo_vs_6mo =  INPUT(rainanm_lag9_1mo_vs_6mo_char,8.);
	rainanm_lag12_1mo_vs_6mo = INPUT(rainanm_lag12_1mo_vs_6mo_char,8.);
	rainfall_avg_lag0_past_1mo = INPUT(rainfall_avg_lag0_past_1mo_char,8.);
	rainfall_avg_lag0_past_3mo = INPUT(rainfall_avg_lag0_past_3mo_char,8.);
	rainfall_avg_lag0_past_6mo = INPUT(rainfall_avg_lag0_past_6mo_char,8.);
	rainfall_avg_lag6_past_1mo = INPUT(rainfall_avg_lag6_past_1mo_char,8.);
	rainfall_avg_lag6_past_3mo = INPUT(rainfall_avg_lag6_past_3mo_char,8.);
	rainfall_avg_lag6_past_6mo = INPUT(rainfall_avg_lag6_past_6mo_char,8.);
	rainfall_avg_lag9_past_1mo = INPUT(rainfall_avg_lag9_past_1mo_char,8.);
	rainfall_avg_lag9_past_3mo = INPUT(rainfall_avg_lag9_past_3mo_char,8.);
	rainfall_avg_lag9_past_6mo = INPUT(rainfall_avg_lag9_past_6mo_char,8.);
	rainfall_avg_lag12_past_1mo = INPUT(rainfall_avg_lag12_past_1mo_char,8.);
	rainfall_avg_lag12_past_3mo = INPUT(rainfall_avg_lag12_past_3mo_char,8.);
	rainfall_avg_lag12_past_6mo = INPUT(rainfall_avg_lag12_past_6mo_char,8.);
	DROP rainanm_lag0_1mo_vs_3mo_char rainanm_lag6_1mo_vs_3mo_char rainanm_lag9_1mo_vs_3mo_char rainanm_lag12_1mo_vs_3mo_char rainanm_lag0_1mo_vs_6mo_char rainanm_lag6_1mo_vs_6mo_char rainanm_lag9_1mo_vs_6mo_char rainanm_lag12_1mo_vs_6mo_char;
RUN;

PROC SORT DATA=chirps_daily_anomaly;;
	BY sasdate;
RUN;

*Create binary exposure and outcomes for Pf and Po, as well as date, season, and rainfall variables;
DATA screening;
	SET screening;
	*if positive Pf18S qPCR reps is greater than 0, is Pf-positive;
	IF qpcr_reps_pos_screening > 0 THEN pf = 1;
		ELSE IF qpcr_reps_pos_screening = 0 THEN pf = 0;
	*if positive Po18S qPCR reps is greater than 0, is Po-positive;
	IF qpcr_reps_pos_screening_po > 0 THEN po = 1;
		ELSE IF qpcr_reps_pos_screening_po = 0 THEN po = 0;
	*Must remove observations over 7567 (PfTz, PoTz, and unused TranSMIT ID numbers 7566-8000);
	IF SUBSTR(screening_id,1,1) in ("P") THEN DELETE;
	IF INPUT(SUBSTR(screening_id,6,4),5.) > 7556 THEN DELETE;
	*Remove erroneous observations with all missing data, such as missing enrollment date;
	IF enrollment_screen_date = " " THEN DELETE;
	*Label plasmodium-positivity variables;
	LABEL pf = "P. falciparum positivity";
	LABEL po = "P. ovale positivity";
	***Village processing. This data will be absent from published dataset for anonymity;
	***Please contact authors if you would like access to a full dataset;
	*initialize village variables;
	LENGTH village $12.;
	LENGTH catvillage $12.;
	LENGTH region 4.;;
	village = enrollment_village;
	*recode villages to only feature the 15 most common villages, consider additional cleaning;
	IF village = " " THEN catvillage = " ";
		*Fix typos or known sub-villages;
		ELSE IF village = "Fukayosi-mng" THEN catvillage = "Fukayosi";
		ELSE IF village = "Fukayosi-Mng'Ongo" THEN catvillage = "Fukayosi";
		ELSE IF village = "Kwa Mkolea" THEN catvillage = "Kwa Mkorea";
		ELSE IF village = "Kwamkolea" THEN catvillage = "Kwa Mkorea";
		ELSE IF village = "Mtakuja A" THEN catvillage = "Mtakuja";
		ELSE IF village = "Mtakuja B" THEN catvillage = "Mtakuja";
		ELSE IF village = "Mwavi" THEN catvillage = "Mwavi Proper";
		ELSE IF village = "Mwanamvulu" THEN catvillage = "Mwanamvuli";
		*keep catvillage = village for subjects in villages with >=75 subjects;
		ELSE IF village = "Chasimba" THEN catvillage = village;
		ELSE IF village = "Fukayosi" THEN catvillage = village;
		ELSE IF village = "Kiegea" THEN catvillage = village;
		ELSE IF village = "Kijiweni" THEN catvillage = village;
		ELSE IF village = "Kisumbi" THEN catvillage = village;
		ELSE IF village = "Kiwangwa" THEN catvillage = village;
		ELSE IF village = "Kudibaha" THEN catvillage = village;
		ELSE IF village = "Kwa Mkorea" THEN catvillage = village;
		ELSE IF village = "Lwazi" THEN catvillage = village;
		ELSE IF village = "Magoza" THEN catvillage = village;
		ELSE IF village = "Mtakuja" THEN catvillage = village;
		ELSE IF village = "Mwanamvuli" THEN catvillage = village;
		ELSE IF village = "Mwavi Proper" THEN catvillage = village;
		ELSE IF village = "Mwavi A" THEN catvillage = village;
		ELSE IF village = "Mwavi B" THEN catvillage = village;
		ELSE IF village = "Vigwaza" THEN catvillage = village;
		ELSE IF village = "Yombo" THEN catvillage = village;
		ELSE IF village = "Zemba" THEN catvillage = village;
		ELSE catvillage = " ";
	*Create region variable from cleaned village data. 0 = North, 1 = South;
	IF catvillage = " " THEN region = .;
		ELSE IF catvillage = "Chasimba" THEN region = 1;
		ELSE IF catvillage = "Fukayosi" THEN region = 0;
		ELSE IF catvillage = "Kiegea" THEN region = 1;
		ELSE IF catvillage = "Kijiweni" THEN region = 1;
		*Kisumbi seems to be a sub-village of Yombo;
		ELSE IF catvillage = "Kisumbi" THEN region = 1;
		ELSE IF catvillage = "Kiwangwa" THEN region = 0;
		ELSE IF catvillage = "Kudibaha" THEN region = 0;
		*Kwa Mkorea not yet identified in local or regional maps;
		ELSE IF catvillage = "Kwa Mkorea" THEN region = 0;
		ELSE IF catvillage = "Lwazi" THEN region = 1;
		ELSE IF catvillage = "Magoza" THEN region = 1;
		ELSE IF catvillage = "Mtakuja" THEN region = 0;
		*Mwananamvuli may also be spelled as Mwanamvulu, sub-village of Fukayosi;
		ELSE IF catvillage = "Mwanamvuli" THEN region = 0;
		ELSE IF catvillage = "Mwavi Proper" THEN region = 0;
		ELSE IF catvillage = "Mwavi A" THEN region = 0;
		ELSE IF catvillage = "Mwavi B" THEN region = 0;
		ELSE IF catvillage = "Vigwaza" THEN region = 1;
		ELSE IF catvillage = "Yombo" THEN region = 1;
		ELSE IF catvillage = "Zemba" THEN region = 0;
	FORMAT region region_.;
	*Create date variables compatible with sas;
	LENGTH date $8.;
	date = SUBSTR(enrollment_screen_date,3,8);
	sasdate = INPUT(date,YYMMDD8.);
	*add lagged sasdate for use in rainfall;
	sasdatelag6 = sasdate - 42;
	sasdatelag9 = sasdate - 63;
	sasdatelag12 = sasdate - 84;
	*Create unlagged and lagged month and year variables for merging with rainfall datasets;
	month_nolag = MONTH(sasdate);
	year_nolag = YEAR(sasdate);
	month_lag6 = MONTH(sasdatelag6);
	year_lag6 = YEAR(sasdatelag6);
	month_lag9 = MONTH(sasdatelag9);
	year_lag9 = YEAR(sasdatelag9);
	month_lag12 = MONTH(sasdatelag12);
	year_lag12 = YEAR(sasdatelag12);
	*recode screening date to not include specific time stamps;
	screen_date = SUBSTR(enrollment_screen_date,1,10);
	*code ct values so missing for negative samples;
	IF pf ^= 1 THEN pf18s_ct = .;
		ELSE pf18s_ct = qpcr_ct_screening;
	po18s_ct = qpcr_ct_screening_po;
RUN;

*Sort screening data by date for merge;
PROC SORT DATA=screening;
	BY sasdate;
RUN;

*Merge daily rainfall anomaly;
DATA screening;
	MERGE screening (in = in1) chirps_daily_anomaly;
	IF in1;
	BY sasdate;
RUN;

*Add formatting for rainfall variables;
DATA screening;
	SET screening;
	FORMAT season seasonlag6 seasonlag9 seasonlag12 season_.;
	LABEL village = "Local Village";
	LABEL rainfall_avg_lag0_past_1mo = "Average daily rainfall over previous mo (mm/day)";
	LABEL rainfall_avg_lag0_past_3mo = "Average daily rainfall over previous 3 months (mm/day)";
	LABEL rainanm_lag0_1mo_vs_3mo = "Difference in average daily rainfall between previous 1 and 3mo (mm/day)";
	LABEL rainfall_avg_lag6_past_1mo = "Average daily rainfall over previous mo, starting 6w earlier (mm/day)";
	LABEL rainfall_avg_lag6_past_3mo = "Average daily rainfall over previous 3 months, starting 6w earlier  (mm/day)";
	LABEL rainanm_lag6_1mo_vs_3mo = "Difference in average daily rainfall between previous 1 and 3mo, starting 6w earlier  (mm/day)";
	LABEL rainfall_avg_lag9_past_1mo = "Average daily rainfall over previous mo, starting 9w earlier (mm/day)";
	LABEL rainfall_avg_lag9_past_3mo = "Average daily rainfall over previous 3 months, starting 9w earlier  (mm/day)";
	LABEL rainanm_lag9_1mo_vs_3mo = "Difference in average daily rainfall between previous 1 and 3mo, starting 9w earlier  (mm/day)";
	LABEL rainfall_avg_lag12_past_1mo = "Average daily rainfall over previous mo, starting 12w earlier (mm/day)";
	LABEL rainfall_avg_lag12_past_3mo = "Average daily rainfall over previous 3 months, starting 12w earlier  (mm/day)";
	LABEL rainanm_lag12_1mo_vs_3mo = "Difference in average daily rainfall between previous 1 and 3mo, starting 12w earlier  (mm/day)";
	LABEL pf18s_ct = "Pf18S qPCR Ct value";
	LABEL po18s_ct = "Po18S qPCR Ct value";
RUN;

PROC CONTENTS DATA=screening;
RUN;

*determine missingness of variables including Pf and Po;
PROC FREQ DATA=screening;
	TABLES pf / MISSPRINT;
	TABLES po / MISSPRINT;
	TABLES enrollment_gender / MISSPRINT;
	TABLES enrollment_age / MISSPRINT;
	TABLES date / MISSPRINT;
	TABLES catvillage / MISSPRINT;
RUN;

PROC PRINT DATA=screening;
	WHERE (po=. OR pf=.);
RUN;
	

*generate dataset that drops those missing either pf or po;
*will be saved to directory noted as "dir" at top of this file;
DATA dir.transmit_cleaned;
	SET screening;
	IF (po=. OR pf=.) THEN DELETE;
RUN;

*modify dataset to drop individuals with symptoms at screening;
DATA dir.transmit_cleaned;
	SET dir.transmit_cleaned;
	IF (enrollment_symptoms_yn = 1) THEN DELETE;
RUN;

PROC FREQ DATA=dir.transmit_cleaned;
	TABLES pf / MISSPRINT;
	TABLES po / MISSPRINT;
	TABLES enrollment_gender / MISSPRINT;
	TABLES enrollment_age / MISSPRINT;
	TABLES catvillage / MISSPRINT;
RUN;

PROC MEANS DATA=dir.transmit_cleaned;
	VAR enrollment_age;
RUN;

PROC MEANS DATA=dir.transmit_cleaned N P25 P50 P75;
	VAR enrollment_age;
	VAR rainanm_lag0_1mo_vs_3mo;
RUN;



PROC CONTENTS DATA=dir.transmit_cleaned;
RUN;









