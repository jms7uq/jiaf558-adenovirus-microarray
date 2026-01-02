/* adenovirus_microarray_analysis.sas
   JID doi:10.1093/infdis/jiaf558

   This program:
   - Imports GEO-curated sample annotation + processed matrix + (optional) raw matrix
   - Generates Table 1 descriptive statistics (by infection group) with KW/chi-square tests
   - Performs rank-based PCA (METHOD=PRIN) retaining 8 PCs
   - Runs logistic regressions for year-2 infection using PC scores
   - Runs severity analyses using Ruuska scoring derived from TAC episode-level clinical components (available in GEO)
   - Writes outputs to ./outputs/tables
*/

options nodate nonumber mprint mlogic symbolgen;

%let ROOT=.;  /* set to repo root when running in SAS */
%let DATADIR=&ROOT./data;
%let OUTDIR=&ROOT./outputs/tables;

libname out "&OUTDIR.";

/* ---------- Import sample annotation ---------- */
proc import datafile="&DATADIR./GEO_sample_annotation_AdV_year1_serum_Table1_aligned_with_1yr.tsv"
    out=work.sample_annot dbms=dlm replace;
    delimiter='09'x;
    getnames=yes;
    guessingrows=max;
run;

/* Standardize key variables */
data work.sample_annot;
    set work.sample_annot;
    /* ensure numeric infection flags */
    year1_infection_num = input(year1_infection, best12.);
    year2_infection_num = input(year2_infection, best12.);
    /* create group variable (as in manuscript) */
    length group $6;
    group = infection_group;
run;

/* ---------- Table 1: Continuous variables (median/IQR) ---------- */
ods listing close;
ods output summary=out.table1_cont_summary;

proc means data=work.sample_annot median q1 q3 n nmiss;
    class infection_group;
    var child_age_upon_enrollment_days
        enrollment_haz
        enrollment_waz
        child_haz_at_1_year
        child_waz_at_1_year
        maternal_age_upon_enrollment_years
        maternal_bmi_upon_enrollment
        household_total_yearly_income_taka
        household_number_of_rooms
        household_number_of_children
        household_number_of_siblings_under_5;
run;

/* Kruskal-Wallis tests across 4 infection groups */
ods output KruskalWallisTest=out.table1_kw_pvalues;
proc npar1way data=work.sample_annot wilcoxon;
    class infection_group;
    var child_age_upon_enrollment_days
        enrollment_haz
        enrollment_waz
        child_haz_at_1_year
        child_waz_at_1_year
        maternal_age_upon_enrollment_years
        maternal_bmi_upon_enrollment
        household_total_yearly_income_taka
        household_number_of_rooms
        household_number_of_children
        household_number_of_siblings_under_5;
run;

/* ---------- Table 1: Categorical variables (n/%) + chi-square ---------- */
ods output CrossTabFreqs=out.table1_cat_crosstabs
           ChiSq=out.table1_cat_chisq;

proc freq data=work.sample_annot;
    tables infection_group*(child_sex breastfed_at_enrollment
           maternal_education_grouped maternal_occupation_grouped
           open_drain_beside_house toilet_shared_with_other_households boiled_or_treated_household_water) / chisq;
run;

ods listing;

/* ---------- Import processed matrix and transpose to sample x antigen ---------- */
/* Processed matrix is features (ID_REF) x samples (S_####) */
proc import datafile="&DATADIR./GEO_processed_matrix_AdV_year1_serum.tsv"
    out=work.proc_matrix_long dbms=dlm replace;
    delimiter='09'x;
    getnames=yes;
    guessingrows=max;
run;

/* Transpose: one row per sample, one column per antigen.
   Create SAS-safe variable names from ID_REF (<=32 chars) and ensure uniqueness.
*/
data work.proc_idmap;
    set work.proc_matrix_long(keep=ID_REF);
    length id_sas $32;
    id_sas = compress(translate(ID_REF, '_', '-'),,'kw'); /* replace hyphen with underscore; keep word chars */
    if length(id_sas) > 32 then id_sas = substr(id_sas,1,32);
run;

proc sort data=work.proc_idmap nodupkey; by ID_REF; run;

/* Assign unique suffixes if collisions occur */
proc sort data=work.proc_idmap; by id_sas ID_REF; run;
data work.proc_idmap;
    set work.proc_idmap;
    by id_sas;
    retain suffix;
    if first.id_sas then suffix=0;
    else suffix+1;
    if suffix>0 then id_sas = cats(substr(id_sas,1,28), "_", put(suffix, z3.));
run;

/* Merge map back and reshape */
proc sql;
    create table work.proc_mapped as
    select a.*, b.id_sas
    from work.proc_matrix_long as a
    left join work.proc_idmap as b
    on a.ID_REF = b.ID_REF;
quit;

/* Convert to long then transpose: 
   First, turn wide (one row per antigen) into long (sample,value)
*/
proc transpose data=work.proc_mapped out=work.proc_long name=sample_id;
    by ID_REF id_sas;
    var S_:;
run;

data work.proc_long;
    set work.proc_long;
    /* _NAME_ contains sample column name; sample_id already contains that */
    value = col1;
    keep sample_id id_sas value;
run;

/* Transpose to sample x antigen */
proc transpose data=work.proc_long out=work.proc_wide(drop=_NAME_) prefix=Ag_;
    by sample_id;
    id id_sas;
    var value;
run;

/* Merge outcomes */
proc sql;
    create table work.analysis as
    select a.*, b.year1_infection_num as year1_infection,
              b.year2_infection_num as year2_infection,
              b.infection_group,
              b.child_sex,
              b.enrollment_haz,
              b.maternal_bmi_upon_enrollment,
              b.toilet_shared_with_other_households,
              b.child_waz_at_1_year
    from work.proc_wide as a
    left join work.sample_annot as b
    on a.sample_id = b.sample_id;
quit;

/* ---------- Rank-transform antibody variables then PCA (retain 8 PCs) ---------- */
proc contents data=work.analysis out=work._vars noprint; run;

/* Create macro variable list of antigen columns (all that start with Ag_) */
proc sql noprint;
    select name into :AGLIST separated by ' '
    from work._vars
    where upcase(name) like 'AG_%';
quit;

proc rank data=work.analysis out=work.ranked_data;
    var &AGLIST.;
run;

proc princomp data=work.ranked_data out=work.pca_scores n=8;
    var &AGLIST.;
run;

/* Merge PCA scores back with covariates */
proc sql;
    create table work.pca_out as
    select a.*, b.year1_infection_num as year1_infection,
              b.year2_infection_num as year2_infection,
              b.maternal_bmi_upon_enrollment, b.child_sex,
              b.enrollment_haz, b.toilet_shared_with_other_households,
              b.child_waz_at_1_year
    from work.pca_scores as a
    left join work.sample_annot as b
    on a.sample_id = b.sample_id;
quit;

/* ---------- Logistic regression: year2 infection ~ PCs + year1 infection ---------- */
ods output ParameterEstimates=out.logit_pcs_unadjusted;
proc logistic data=work.pca_out descending;
    model year2_infection = year1_infection Prin1-Prin8;
run;

ods output ParameterEstimates=out.logit_pcs_adjusted;
proc logistic data=work.pca_out descending;
    model year2_infection = year1_infection Prin1-Prin8
        maternal_bmi_upon_enrollment
        child_sex
        enrollment_haz
        toilet_shared_with_other_households
        child_waz_at_1_year;
run;

/* =======================================================================
   Severity analyses (Ruuska score) – derived from TAC episode components
   ======================================================================= */

proc import datafile="&DATADIR./GEO_TAC_episode_level_clean.tsv"
    out=work.episodes dbms=dlm replace;
    delimiter='09'x;
    getnames=yes;
    guessingrows=max;
run;

data work.episodes_ruu;
    set work.episodes;
    length sample_id $20;
    sample_id = cats("S_", SID);

    if Diarr_Specage <= 365 then year=1;
    else if Diarr_Specage <= 730 then year=2;
    else year=.;

    adv_ct  = input(Adenovirus_40_41, best12.);
    adv_afe = input(adenovirus_40_41_afe, best12.);
    pcr_pos = (adv_ct ne . and adv_ct < 35);

    ruu_dur_diarr = input(DURADIAR, best12.);
    ruu_max_diarr = input(MAXDIAR, best12.);
    ruu_dur_vomit = input(DURAVOMI, best12.);
    ruu_max_vomit = input(MAXVOMI, best12.);
    ruu_fever_raw = input(TEMP, best12.);
    ruu_dehyd_raw = input(DEHYD, best12.);
    ruu_treat_raw = input(REHYD, best12.);

    if ruu_dur_vomit = . or ruu_dur_vomit = 9 then ruu_dur_vomit = 0;
    if ruu_max_vomit = . or ruu_max_vomit = 9 then ruu_max_vomit = 0;

    if ruu_fever_raw = 4 then ruu_fever = .;
    else ruu_fever = ruu_fever_raw;

    if ruu_dehyd_raw = 1 then ruu_dehyd = 0;
    else if ruu_dehyd_raw = 2 then ruu_dehyd = 2;
    else ruu_dehyd = .;

    if ruu_treat_raw = 1 then ruu_treat = 0;
    else if ruu_treat_raw = 2 then ruu_treat = 1;
    else if ruu_treat_raw = 3 then ruu_treat = 2;
    else ruu_treat = .;

    ruuska_score = sum(ruu_dur_diarr, ruu_max_diarr, ruu_dur_vomit, ruu_max_vomit, ruu_fever, ruu_dehyd, ruu_treat);

    length ruuska_cat $20;
    if ruuska_score = . then ruuska_cat = "";
    else if ruuska_score <= 6 then ruuska_cat = "Mild (0-6)";
    else if ruuska_score <= 10 then ruuska_cat = "Moderate (7-10)";
    else ruuska_cat = "Severe (>=11)";
run;

/* Select lowest-Ct AdV40/41 episode per child-year */
proc sort data=work.episodes_ruu;
    where pcr_pos=1 and year in (1,2);
    by SID year adv_ct;
run;

data work.adv_epi_selected;
    set work.episodes_ruu;
    where pcr_pos=1 and year in (1,2);
    by SID year adv_ct;
    if first.year then output;
run;

proc sql;
    create table work.adv_epi_selected2 as
    select a.*, b.infection_group, b.year1_infection_num as year1_infection,
              b.year2_infection_num as year2_infection
    from work.adv_epi_selected as a
    left join work.sample_annot as b
    on a.sample_id = b.sample_id;
quit;

data work.year2_severity;
    set work.adv_epi_selected2;
    where year=2 and year2_infection=1;
    keep sample_id SID Diarr_Specage adv_ct adv_afe ruuska_score ruuska_cat infection_group;
run;

proc export data=work.year2_severity outfile="&OUTDIR./year2_adv_selected_episode_severity.csv" dbms=csv replace; run;

/* Import platform annotation so we can identify penton base / fiber targets */
proc import datafile="&DATADIR./GEO_platform_annotation_AdV_array.tsv"
    out=work.platform dbms=dlm replace;
    delimiter='09'x;
    getnames=yes;
    guessingrows=max;
run;

/* Long-form processed matrix with ID_REF preserved */
proc transpose data=work.proc_mapped out=work.proc_long2 name=sample_id;
    by ID_REF;
    var S_:;
run;

data work.proc_long2;
    set work.proc_long2;
    value = col1;
    keep sample_id ID_REF value;
run;

proc sql;
    create table work.proc_long_annot as
    select a.sample_id, a.ID_REF, a.value, b.Description
    from work.proc_long2 as a
    left join work.platform as b
    on a.ID_REF = b.ID_REF;
quit;

data work.top_targets;
    set work.proc_long_annot;
    length target_key $32;
    desc_up = upcase(Description);

    if index(desc_up, "HADV-40") and index(desc_up, "PENTON BASE") then target_key="AdV40_PB";
    else if index(desc_up, "HADV-41") and index(desc_up, "PENTON BASE") then target_key="AdV41_PB";
    else if index(desc_up, "HADV-40") and (index(desc_up, "FIBER-2") or index(desc_up,"SHORT FIBER")) then target_key="AdV40_short_fiber";
    else if index(desc_up, "HADV-41") and index(desc_up, "SHORT FIBER") then target_key="AdV41_short_fiber";
    else if index(desc_up, "HADV-41") and index(desc_up, "LONG FIBER") then target_key="AdV41_long_fiber";
    else target_key="";

    if target_key ne "";
run;

proc sort data=work.top_targets; by sample_id target_key; run;

proc means data=work.top_targets noprint median;
    by sample_id target_key;
    var value;
    output out=work.top_targets_summ(drop=_TYPE_ _FREQ_) median=si_median;
run;

proc sql;
    create table work.top_targets_year2 as
    select a.*, b.ruuska_cat, b.ruuska_score, b.infection_group
    from work.top_targets_summ as a
    inner join work.year2_severity as b
    on a.sample_id = b.sample_id;
quit;

ods output KruskalWallisTest=out.year2_severity_kw;
proc npar1way data=work.top_targets_year2 wilcoxon;
    class ruuska_cat;
    var si_median;
    by target_key;
run;

proc export data=work.top_targets_year2 outfile="&OUTDIR./top_targets_year2_severity_long.csv" dbms=csv replace; run;
proc export data=out.year2_severity_kw outfile="&OUTDIR./top_targets_year2_severity_kw.csv" dbms=csv replace; run;

/* ---------- Save key outputs ---------- */
proc export data=out.table1_cont_summary outfile="&OUTDIR./table1_cont_summary.csv" dbms=csv replace; run;
proc export data=out.table1_kw_pvalues outfile="&OUTDIR./table1_kw_pvalues.csv" dbms=csv replace; run;
proc export data=out.table1_cat_crosstabs outfile="&OUTDIR./table1_cat_crosstabs.csv" dbms=csv replace; run;
proc export data=out.table1_cat_chisq outfile="&OUTDIR./table1_cat_chisq.csv" dbms=csv replace; run;

proc export data=out.logit_pcs_unadjusted outfile="&OUTDIR./logit_pcs_unadjusted.csv" dbms=csv replace; run;
proc export data=out.logit_pcs_adjusted outfile="&OUTDIR./logit_pcs_adjusted.csv" dbms=csv replace; run;
