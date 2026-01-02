# Adenovirus 40/41 microarray analysis (JID doi:10.1093/infdis/jiaf558)

This repository contains reproducible analysis code for the JID manuscript **“Adenovirus 40 and 41 Antibodies Associated With Protection From Infection in a Bangladeshi Birth Cohort”** (doi:10.1093/infdis/jiaf558).

## What’s here
- **R**: all manuscript figures + random forest analysis (caret + randomForest) and supporting analyses (PCA, antibody tertiles, correlations, coinfections).
- **SAS**: Table 1 descriptive statistics, PCA, and logistic regressions (unadjusted/adjusted) as described in the manuscript.

## Data inputs
This code expects the GEO-curated tabular files (TSV/CSV) that you prepared for submission:
- `GEO_processed_matrix_AdV_year1_serum.tsv`
- `GEO_raw_matrix_AdV_year1_serum.tsv`
- `GEO_platform_annotation_AdV_array.tsv`
- `GEO_sample_annotation_AdV_year1_serum_Table1_aligned_with_1yr.tsv`
- `GEO_TAC_episode_level_clean.tsv`
- `GEO_TAC_column_name_map.csv` (optional)
- `GEO_Clinical_Data_Dictionary_AllSheets.tsv` (optional)

Place these files in `data/` (or change paths in the scripts).



## How to run (R)
1. Install packages (first run may take a few minutes):
   ```r
   install.packages(c(
     "tidyverse","data.table","readr","janitor","ggplot2","patchwork","broom",
     "caret","randomForest","pheatmap","forcats","scales"
   ))
   ```
2. Run:
   ```r
   source("R/analysis_adenovirus_microarray.R")
   ```
Figures and tables will be written to `outputs/figures/` and `outputs/tables/`.

## How to run (SAS)
Open and run:
- `sas/adenovirus_microarray_analysis.sas`

This script imports the same GEO-curated files and writes results tables to `outputs/tables/`.

## Repositories (for RPPR / DMSP)
- Data: **NIH GEO** (accession pending / to be inserted when assigned)
- Code: this GitHub repo (recommended to archive a release to Zenodo for a DOI)
