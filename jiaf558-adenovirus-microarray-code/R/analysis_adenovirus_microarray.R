#Rscript

# Adenovirus 40/41 microarray analysis (JID doi:10.1093/infdis/jiaf558)
# - Generates manuscript figures and runs random forest analysis
# - Uses GEO-curated tabular inputs (processed + raw matrices, platform annotation, sample annotation, TAC episode-level table)
#
# Expected inputs (default: ./data):
#   GEO_processed_matrix_AdV_year1_serum.tsv
#   GEO_raw_matrix_AdV_year1_serum.tsv
#   GEO_platform_annotation_AdV_array.tsv
#   GEO_sample_annotation_AdV_year1_serum_Table1_aligned_with_1yr.tsv
#   GEO_TAC_episode_level_clean.tsv
#
# Outputs (default: ./outputs):
#   outputs/figures/*.png
#   outputs/tables/*.csv
#
# Manuscript-guided analysis:
# - PCA on ranked antibody reactivities; retain leading PCs (n=8) and evaluate association with year-2 AdV40/41 infection.
# - Logistic regression: year-2 infection ~ PC scores + year-1 infection; adjusted models include covariates (sex, enrollment HAZ, shared toilet).
# - Forest plot of PC odds ratios (unadjusted/adjusted) with significance highlighting.
# - Top PC loadings for PC2 include AdV40/41 penton base and fiber targets; compare across 4 infection groups (KW tests).
# - Boxplots of top target reactivities by infection group; volcano-style plots for tertile associations.
# - Timeline plot of TAC episodes and serum collection.
# - Coinfection summaries, including AdV episode counts and pathogen AFe coinfection rates by year/group.
# - Antibody comparisons by coinfection status (year-2 AdV+ episodes) and supplemental histograms for top antigens.
# - Random forest with 10-fold CV, grid search for ntree/mtry; evaluate variable importance.
#
# References: Hendrick et al. J Infect Dis. doi:10.1093/infdis/jiaf558

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(janitor)
  library(broom)
  library(caret)
  library(randomForest)
  library(pheatmap)
  library(forcats)
  library(scales)
  library(patchwork)
  library(ggrepel)
})

# ------------------------
# User-configurable paths
# ------------------------
DATA_DIR   <- "data"
OUT_DIR    <- "outputs"
FIG_DIR    <- file.path(OUT_DIR, "figures")
TAB_DIR    <- file.path(OUT_DIR, "tables")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

FILES <- list(
  processed_matrix = file.path(DATA_DIR, "GEO_processed_matrix_AdV_year1_serum.tsv"),
  raw_matrix       = file.path(DATA_DIR, "GEO_raw_matrix_AdV_year1_serum.tsv"),
  platform         = file.path(DATA_DIR, "GEO_platform_annotation_AdV_array.tsv"),
  sample_annot      = file.path(DATA_DIR, "GEO_sample_annotation_AdV_year1_serum_Table1_aligned_with_1yr.tsv"),
  tac_episode      = file.path(DATA_DIR, "GEO_TAC_episode_level_clean.tsv")
)

# ------------------------
# Helper functions
# ------------------------
stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path, call. = FALSE)
}

read_tsv_dt <- function(path) {
  data.table::fread(path, sep = "\t", data.table = FALSE, showProgress = FALSE) %>% as_tibble()
}

# Convert reactivity matrix (ID_REF + S_#### columns) to long format joined to annotation
matrix_to_long <- function(mat_df) {
  mat_df %>%
    pivot_longer(cols = -ID_REF, names_to = "sample_id", values_to = "value")
}

# Make tertiles (0/1/2 -> label)
make_tertiles <- function(x) {
  # x numeric, return ordered factor with levels low/mid/high
  q <- quantile(x, probs = c(1/3, 2/3), na.rm = TRUE)
  cut(x, breaks = c(-Inf, q[1], q[2], Inf), labels = c("T1_low","T2_mid","T3_high"), include.lowest = TRUE, ordered_result = TRUE)
}

# Spearman correlation heatmap convenience
save_heatmap <- function(mat, file, main = NULL) {
  png(file, width = 1400, height = 1100, res = 150)
  pheatmap::pheatmap(
    mat,
    main = main,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation"
  )
  dev.off()
}

# ------------------------
# Load inputs
# ------------------------
purrr::walk(FILES, stop_if_missing)

processed <- read_tsv_dt(FILES$processed_matrix) %>% clean_names()
raw_mat   <- read_tsv_dt(FILES$raw_matrix) %>% clean_names()
platform  <- read_tsv_dt(FILES$platform) %>% clean_names()
samples   <- read_tsv_dt(FILES$sample_annot) %>% clean_names()
episodes  <- read_tsv_dt(FILES$tac_episode) %>% clean_names()

# Standardize key column names
# Matrices: id_ref + sample columns "s_####"
processed <- processed %>% rename(id_ref = id_ref)
raw_mat   <- raw_mat %>% rename(id_ref = id_ref)

# Sample annotation must include: sample_id, infection_group, year1_infection, year2_infection
samples <- samples %>%
  mutate(
    infection_group = factor(infection_group, levels = c("Y1+Y2-","Y1-Y2+","Y1+Y2+","Y1-Y2-")),
    year1_infection = as.integer(year1_infection),
    year2_infection = as.integer(year2_infection)
  )

# ------------------------
# Antigen set definitions
# ------------------------
# Targets of interest (Table/Fig references): penton base & fiber targets
platform_targets <- platform %>%
  mutate(desc = tolower(description)) %>%
  mutate(
    is_penton = str_detect(desc, "penton base"),
    is_fiber  = str_detect(desc, "fiber"),
    is_short_fiber = str_detect(desc, "short fiber|fiber-2"),
    is_long_fiber  = str_detect(desc, "long fiber"),
    is_capsid_piiia = str_detect(desc, "p\\s*iiia|piiia|iii[a]")
  )

# Canonical “top” AdV 40/41 targets (matches the started code list and manuscript narrative)
# (AdV40 short fiber corresponds to HAdV-40 L5A fiber-2 in this platform)
top_targets_regex <- list(
  adv40_penton = "HAdV-40.*penton base",
  adv41_penton = "HAdV-41.*penton base",
  adv41_short_fiber = "HAdV-41.*short fiber",
  adv41_long_fiber  = "HAdV-41.*long fiber",
  adv40_short_fiber = "HAdV-40.*fiber-2",
  adv41_iiia        = "HAdV-41.*IIIa|HAdV-41.*pIIIa|HAdV-41.*III A",
  adv40_capsid_piiia = "HAdV-40.*pIIIa|HAdV-40.*IIIa|HAdV-40.*capsid.*IIIa"
)

select_targets <- function(regex) {
  platform %>%
    filter(str_detect(description, regex) | str_detect(id_ref, regex)) %>%
    distinct(id_ref, description, spot_type, virus_hint)
}

top_targets <- bind_rows(
  purrr::imap(top_targets_regex, ~ select_targets(.x) %>% mutate(target_key = .y))
) %>% distinct()

write_csv(top_targets, file.path(TAB_DIR, "top_targets_selected.csv"))

# ------------------------
# Prepare matrices joined to sample annotation
# ------------------------
proc_long <- processed %>%
  rename_with(~ str_replace_all(.x, "^s_", "S_")) %>% # support if cleaned to lowercase
  { . } %>%
  { 
    # identify sample columns robustly
    samp_cols <- names(.)[names(.) != "id_ref"]
    colnames(.) <- c("id_ref", samp_cols)
    .
  } %>%
  matrix_to_long() %>%
  rename(sample_id = sample_id) %>%
  left_join(samples, by = "sample_id") %>%
  left_join(platform %>% select(id_ref, description, spot_type, virus_hint), by = "id_ref")

# Wide format for machine learning: samples x features (processed)
proc_wide <- processed %>%
  rename_with(~ str_replace_all(.x, "^s_", "S_"))

# ------------------------
# PCA on ranked antibody reactivities
# ------------------------
# Rank-transform each feature across samples (ties averaged), then PCA
proc_mat <- processed %>%
  rename_with(~ str_replace_all(.x, "^s_", "S_")) %>%
  column_to_rownames("id_ref") %>%
  as.matrix()

# rows = features; cols = samples. Rank within feature across samples:
ranked_mat <- apply(proc_mat, 1, function(v) rank(v, ties.method = "average", na.last = "keep"))
ranked_mat <- t(ranked_mat) # now: samples x features

# PCA
pca <- prcomp(ranked_mat, center = TRUE, scale. = TRUE)
pca_var <- (pca$sdev^2) / sum(pca$sdev^2)
pca_scree <- tibble(PC = paste0("PC", seq_along(pca_var)), var_explained = pca_var)

write_csv(pca_scree, file.path(TAB_DIR, "pca_variance_explained.csv"))

# Keep first 8 PCs (manuscript)
pcs_keep <- 1:8
pc_scores <- as_tibble(pca$x[, pcs_keep, drop = FALSE]) %>%
  mutate(sample_id = rownames(pca$x)) %>%
  left_join(samples, by = "sample_id")

# Logistic regression: year2 infection ~ PC1..PC8 + year1 infection (unadjusted and adjusted)
fit_logit <- function(df, rhs, label) {
  f <- as.formula(paste("year2_infection ~", rhs))
  m <- glm(f, data = df, family = binomial())
  broom::tidy(m, conf.int = TRUE, exponentiate = TRUE) %>%
    mutate(model = label)
}

rhs_unadj <- paste(c(paste0("PC", 1:8), "year1_infection"), collapse = " + ")
rhs_adj   <- paste(c(paste0("PC", 1:8), "year1_infection",
                     "maternal_bmi_upon_enrollment",
                     "child_sex",
                     "enrollment_haz",
                     "toilet_shared_with_other_households",
                     "child_waz_at_1_year"), collapse = " + ")

# Some covariates may be missing in public metadata; handle gracefully
pc_scores2 <- pc_scores %>%
  mutate(
    child_sex = factor(child_sex),
    toilet_shared_with_other_households = factor(toilet_shared_with_other_households)
  )

# Drop variables not present
safe_rhs <- function(df, rhs) {
  vars <- str_split(rhs, "\\s*\\+\\s*")[[1]] %>% str_trim()
  vars_present <- vars[vars %in% names(df)]
  paste(vars_present, collapse = " + ")
}

rhs_unadj2 <- safe_rhs(pc_scores2, rhs_unadj)
rhs_adj2   <- safe_rhs(pc_scores2, rhs_adj)

logit_unadj <- fit_logit(pc_scores2, rhs_unadj2, "Unadjusted")
logit_adj   <- fit_logit(pc_scores2, rhs_adj2, "Adjusted")

pc_logit_out <- bind_rows(logit_unadj, logit_adj) %>%
  filter(term != "(Intercept)") %>%
  mutate(term = str_replace(term, "year1_infection", "Year 1 infection"),
         term = factor(term, levels = c(paste0("PC", 1:8), "Year 1 infection")))

write_csv(pc_logit_out, file.path(TAB_DIR, "logistic_pcs_year2infection.csv"))

# Forest plot of PC ORs (Supplementary)
p_forest <- pc_logit_out %>%
  mutate(
    significant = if_else(p.value < 0.05, "Significant", "Not Significant"),
    term_label = if_else(term == "Year 1 infection", "Year 1 Infection Status", as.character(term)),
    term_label = factor(term_label, levels = c(paste0("PC", 1:8), "Year 1 Infection Status"))
  ) %>%
  ggplot(aes(x = term_label, y = estimate, ymin = conf.low, ymax = conf.high, shape = model, color = significant)) +
  geom_pointrange(position = position_dodge(width = 0.6)) +
  geom_hline(yintercept = 1, linetype = 2) +
  coord_flip() +
  scale_y_log10() +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  labs(x = NULL, y = "Odds ratio (95% CI)", color = NULL, shape = NULL) +
  theme_minimal(base_size = 12)

ggsave(file.path(FIG_DIR, "supp_forestplot_pcs_year2infection.png"), p_forest, width = 8, height = 6, dpi = 300)

# Identify top loadings for PC2 (manuscript highlights PC2)
loadings <- as_tibble(pca$rotation[, pcs_keep, drop = FALSE], rownames = "id_ref") %>%
  left_join(platform %>% select(id_ref, description), by = "id_ref")

pc2_top <- loadings %>%
  arrange(desc(abs(PC2))) %>%
  slice_head(n = 25)

write_csv(pc2_top, file.path(TAB_DIR, "pc2_top25_loadings.csv"))

# ------------------------
# Figure: distribution of top target reactivities by infection group
# ------------------------
# Extract per-sample values for the selected top targets (by regex matches above)
top_ids <- top_targets$id_ref %>% unique()

top_vals <- proc_long %>%
  filter(id_ref %in% top_ids) %>%
  mutate(target_key = top_targets$target_key[match(id_ref, top_targets$id_ref)])

# Summarize + KW p-values
kw <- top_vals %>%
  group_by(target_key) %>%
  summarize(p_kw = kruskal.test(value ~ infection_group)$p.value, .groups = "drop") %>%
  mutate(p_kw = signif(p_kw, 3))

write_csv(kw, file.path(TAB_DIR, "kruskalwallis_top_targets_by_group.csv"))

p_top_targets <- top_vals %>%
  left_join(kw, by = "target_key") %>%
  mutate(target_key = fct_reorder(target_key, value, .fun = median, na.rm = TRUE)) %>%
  ggplot(aes(x = infection_group, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  facet_wrap(~ target_key, scales = "free_y") +
  labs(x = NULL, y = "Normalized signal intensity (SI)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(FIG_DIR, "fig_top_targets_boxplots_by_infection_group.png"),
       p_top_targets, width = 11, height = 7, dpi = 300)

# ------------------------
# Figure 1-like: timeline plot of episodes + serum
# ------------------------
# Episode file structure is TAC-centric; key columns used in plot:
# - sid (child)
# - diarr_specage (age at stool / diarrhea episode in days)
# - adenovirus_40_41 (Ct value; 40 often indicates negative)
# - adenovirus_40_41_afe (AFe)
episodes2 <- episodes %>%
  rename_with(~ str_replace_all(.x, "\\.", "_")) %>%
  mutate(sample_id = paste0("S_", sid)) %>%
  left_join(samples %>% select(sample_id, infection_group, serum_age_days), by = "sample_id") %>%
  mutate(
    year = case_when(
      diarr_specage <= 365 ~ 1L,
      diarr_specage <= 730 ~ 2L,
      TRUE ~ NA_integer_
    ),
    adv_ct = suppressWarnings(as.numeric(adenovirus_40_41)),
    adv_afe = suppressWarnings(as.numeric(adenovirus_40_41_afe)),
    pcr_pos = !is.na(adv_ct) & adv_ct < 35,
    afe_pos = !is.na(adv_afe) & adv_afe >= 0.5,
    episode_type = case_when(
      pcr_pos & afe_pos ~ "PCR+ & AFe+",
      pcr_pos & !afe_pos ~ "PCR+ only",
      TRUE ~ "PCR-"
    )
  )

# basic child range for segments (0 to max observed episode age)
child_ranges <- episodes2 %>%
  group_by(sample_id, infection_group) %>%
  summarize(min_age = min(diarr_specage, na.rm = TRUE),
            max_age = max(diarr_specage, na.rm = TRUE),
            .groups = "drop")

# Plot
p_timeline <- ggplot() +
  geom_segment(data = child_ranges,
               aes(x = 0, xend = max_age, y = fct_reorder(sample_id, max_age), yend = fct_reorder(sample_id, max_age)),
               linewidth = 0.25, alpha = 0.4) +
  geom_point(data = episodes2,
             aes(x = diarr_specage, y = fct_reorder(sample_id, diarr_specage), shape = episode_type),
             alpha = 0.7, size = 1.3) +
  geom_point(data = samples,
             aes(x = serum_age_days, y = fct_reorder(sample_id, serum_age_days)),
             color = "red", shape = 17, size = 1.8, inherit.aes = FALSE) +
  facet_wrap(~ infection_group, scales = "free_y") +
  labs(x = "Age (days)", y = "Child (study ID)", shape = "Episode") +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom")

ggsave(file.path(FIG_DIR, "fig_timeline_episodes_serum.png"), p_timeline, width = 11, height = 7, dpi = 300)

# ------------------------
# Random forest: predict year-2 infection from antibody reactivities
# ------------------------
# Prepare wide feature set: samples x features
# - transpose processed matrix so each feature is a column
proc_wide2 <- processed %>%
  rename_with(~ str_replace_all(.x, "^s_", "S_"))

# transpose with base R
feat_mat <- proc_wide2 %>%
  as.data.frame()
rownames(feat_mat) <- feat_mat$id_ref
feat_mat$id_ref <- NULL
feat_mat_t <- as.data.frame(t(as.matrix(feat_mat)))
feat_mat_t$sample_id <- rownames(feat_mat_t)

rf_df <- feat_mat_t %>%
  as_tibble() %>%
  left_join(samples %>% select(sample_id, year2_infection, year1_infection), by = "sample_id") %>%
  mutate(year2_infection = factor(year2_infection, levels = c(0,1), labels = c("No","Yes"))) %>%
  select(year2_infection, everything(), -sample_id)

# caret setup
set.seed(123)
ctrl <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = twoClassSummary, savePredictions = "final")

# Grid based on manuscript started code idea (mtry, ntree)
mtry_grid <- c(4, 6, 8, 10, 12)
# caret's rf uses ntree via ...; we loop over ntree values and select best AUC
ntrees <- c(500, 1000, 1500, 2000, 2500, 3000)

rf_results <- list()
best_auc <- -Inf
best_model <- NULL
best_params <- NULL

for (nt in ntrees) {
  tune <- expand.grid(mtry = mtry_grid)
  model <- train(
    year2_infection ~ .,
    data = rf_df,
    method = "rf",
    metric = "ROC",
    trControl = ctrl,
    tuneGrid = tune,
    ntree = nt,
    importance = TRUE
  )
  # pick best ROC for this ntree
  res <- model$results %>% mutate(ntree = nt)
  rf_results[[as.character(nt)]] <- res
  best_row <- res %>% arrange(desc(ROC)) %>% slice(1)
  if (best_row$ROC[1] > best_auc) {
    best_auc <- best_row$ROC[1]
    best_model <- model
    best_params <- best_row
  }
}

rf_grid_results <- bind_rows(rf_results)
write_csv(rf_grid_results, file.path(TAB_DIR, "random_forest_cv_grid.csv"))

# Fit final RF with best parameters (caret already fitted best_model with best mtry for best ntree)
# Extract variable importance
imp <- varImp(best_model)$importance %>%
  rownames_to_column("feature") %>%
  arrange(desc(Overall)) %>%
  slice_head(n = 20)

write_csv(imp, file.path(TAB_DIR, "random_forest_top20_importance.csv"))

p_imp <- ggplot(imp, aes(x = fct_reorder(feature, Overall), y = Overall)) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = "Variable importance (scaled)") +
  theme_minimal(base_size = 12)

ggsave(file.path(FIG_DIR, "fig_random_forest_variable_importance_top20.png"),
       p_imp, width = 9, height = 7, dpi = 300)

# ------------------------
# Volcano-style plot: association of top targets with year-2 infection
# ------------------------
# Use tertiles for each target
top_wide <- top_vals %>%
  select(sample_id, id_ref, value) %>%
  pivot_wider(names_from = id_ref, values_from = value) %>%
  left_join(samples, by = "sample_id")

# Choose up to 7 targets (as in started code)
top7_ids <- top_ids %>% unique() %>% head(7)

# Create tertiles
for (id in top7_ids) {
  top_wide[[paste0(id, "_tert")]] <- make_tertiles(top_wide[[id]])
}

# Fit unadjusted and adjusted logistic models for each antibody tertile
fit_antibody <- function(df, id) {
  term <- paste0(id, "_tert")
  if (!term %in% names(df)) return(NULL)
  df2 <- df %>% filter(!is.na(.data[[term]]))
  # unadjusted: year2 ~ year1 + tertile
  f1 <- as.formula(paste("year2_infection ~ year1_infection +", term))
  m1 <- glm(f1, data = df2, family = binomial())
  # adjusted
  covars <- c("year1_infection","maternal_bmi_upon_enrollment","child_sex","enrollment_haz","toilet_shared_with_other_households","child_waz_at_1_year")
  covars <- covars[covars %in% names(df2)]
  f2 <- as.formula(paste("year2_infection ~", paste(c(covars, term), collapse = " + ")))
  m2 <- glm(f2, data = df2, family = binomial())
  list(unadj = m1, adj = m2)
}

ant_res <- purrr::map(top7_ids, ~ fit_antibody(top_wide, .x))
names(ant_res) <- top7_ids

tidy_one <- function(mod, id, model_label) {
  broom::tidy(mod, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(str_detect(term, "_tert")) %>%
    mutate(id_ref = id, model = model_label)
}

ant_tidy <- bind_rows(
  purrr::imap(ant_res, ~ tidy_one(.x$unadj, .y, "Unadjusted")),
  purrr::imap(ant_res, ~ tidy_one(.x$adj,   .y, "Adjusted"))
) %>%
  left_join(platform %>% select(id_ref, description), by = "id_ref") %>%
  mutate(log_or = log(estimate),
         neglog10_p = -log10(p.value),
         significant = if_else(p.value < 0.05, "Significant", "Not Significant"))

write_csv(ant_tidy, file.path(TAB_DIR, "logistic_top7_targets_tertiles.csv"))

p_volcano <- ggplot(ant_tidy, aes(x = log_or, y = neglog10_p, shape = model, color = significant)) +
  geom_point(alpha = 0.8, size = 2.6) +
  ggrepel::geom_text_repel(aes(label = description), size = 3, max.overlaps = 30) +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "black")) +
  labs(x = "log(OR) (tertile effect)", y = "-log10(p)", color = NULL, shape = "Model p value") +
  theme_minimal(base_size = 12)

ggsave(file.path(FIG_DIR, "fig_volcano_top7_targets_tertiles.png"), p_volcano, width = 7, height = 6, dpi = 300)

# ------------------------
# Coinfection definition (AdV40/41 + >=1 additional pathogen with AFe >= 0.5)
# ------------------------
afe_cols <- names(episodes2) %>% keep(~ str_detect(.x, "_afe$"))
afe_other <- setdiff(afe_cols, "adenovirus_40_41_afe")

episodes2 <- episodes2 %>%
  mutate(
    coinfection = if_else(
      afe_pos & rowSums(across(all_of(afe_other), ~ as.numeric(.x) >= 0.5), na.rm = TRUE) > 0,
      TRUE, FALSE
    )
  )

coinf_by_group <- episodes2 %>%
  filter(afe_pos) %>%
  group_by(infection_group, year) %>%
  summarize(coinfection_rate = mean(coinfection, na.rm = TRUE), n = n(), .groups = "drop")

write_csv(coinf_by_group, file.path(TAB_DIR, "coinfection_rate_by_infection_group.csv"))

# AdV episode counts per child-year (AdV+ episodes only)
adv_episodes <- episodes2 %>%
  filter(pcr_pos, !is.na(year)) %>%
  group_by(sid, year, infection_group) %>%
  summarize(adv_pos_episodes = n(), .groups = "drop") %>%
  mutate(more_than_one = adv_pos_episodes > 1)

adv_episode_summary <- adv_episodes %>%
  group_by(infection_group, year) %>%
  summarize(
    total_children = n(),
    children_more_than_one = sum(more_than_one),
    percent_more_than_one = round(100 * children_more_than_one / total_children, 1),
    median_episodes = median(adv_pos_episodes),
    iqr_lower = quantile(adv_pos_episodes, 0.25),
    iqr_upper = quantile(adv_pos_episodes, 0.75),
    .groups = "drop"
  )

write_csv(adv_episode_summary, file.path(TAB_DIR, "adv_episode_counts_by_group_year.csv"))

# Coinfection counts by number of AFe-positive pathogens among AdV+ episodes
coinfection_counts <- episodes2 %>%
  filter(pcr_pos) %>%
  mutate(across(all_of(afe_other), ~ if_else(as.numeric(.x) >= 0.5, 1, 0))) %>%
  rowwise() %>%
  mutate(coinfection_count = sum(c_across(all_of(afe_other)), na.rm = TRUE)) %>%
  ungroup()

coinfection_count_summary <- coinfection_counts %>%
  count(coinfection_count) %>%
  arrange(desc(coinfection_count))

write_csv(coinfection_count_summary, file.path(TAB_DIR, "coinfection_count_summary.csv"))

# Antibody comparison by coinfection status (year-2 AdV+ episodes)
antibody_vars <- c(
  "HAdV_40_E2A_L_single_stranded_DNA_binding_pro",
  "HAdV_40_IVa2_encap_pro_IVa2",
  "HAdV_40_L1_capsid_pro_precur_pIIIa",
  "HAdV_40_L1_encap_pro_52K",
  "HAdV_40_L2_penton_base",
  "HAdV_40_L5A_fiber_2",
  "HAdV_40_L5_fiber",
  "HAdV_41_E2A_DBP",
  "HAdV_41_IIIa",
  "HAdV_41_long_fiber_pro",
  "HAdV_41_penton_base",
  "HAdV_41_short_fiber_pro",
  "HAdV_41_L2_pMu"
)

antibody_vars_present <- intersect(antibody_vars, processed$id_ref)
df_ab <- proc_long %>%
  filter(id_ref %in% antibody_vars_present) %>%
  select(sample_id, id_ref, value) %>%
  pivot_wider(names_from = id_ref, values_from = value)

coinfection_status_year2 <- episodes2 %>%
  filter(year == 2L, pcr_pos) %>%
  mutate(across(all_of(afe_other), ~ if_else(as.numeric(.x) >= 0.5, 1, 0))) %>%
  rowwise() %>%
  mutate(coinfection_count = sum(c_across(all_of(afe_other)), na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(sample_id) %>%
  summarize(any_coinfection = any(coinfection_count >= 1), .groups = "drop")

antibody_compare_df <- coinfection_status_year2 %>%
  left_join(df_ab, by = "sample_id") %>%
  mutate(coinfection_group = if_else(any_coinfection, "Coinfection", "No Coinfection"))

antibody_long <- antibody_compare_df %>%
  pivot_longer(cols = all_of(antibody_vars_present), names_to = "antigen", values_to = "antibody_level")

antibody_tests <- antibody_long %>%
  group_by(antigen) %>%
  summarize(
    p_value = {
      grp_n <- n_distinct(coinfection_group)
      if (grp_n < 2) NA_real_ else suppressWarnings(wilcox.test(antibody_level ~ coinfection_group)$p.value)
    },
    .groups = "drop"
  ) %>%
  mutate(
    p_signif = case_when(
      p_value <= 0.001 ~ "***",
      p_value <= 0.01 ~ "**",
      p_value <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

write_csv(antibody_tests, file.path(TAB_DIR, "coinfection_antibody_wilcox.csv"))

p_coinfection_antibodies <- ggplot(antibody_long, aes(x = coinfection_group, y = antibody_level)) +
  geom_boxplot() +
  facet_wrap(~ antigen, scales = "free_y") +
  labs(x = NULL, y = "Antibody level") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(FIG_DIR, "fig_coinfection_antibody_boxplots.png"),
       p_coinfection_antibodies, width = 12, height = 8, dpi = 300)

# ------------------------
# Correlation heatmap among external capsid targets (example heuristic)
# ------------------------
external_capsid_ids <- platform_targets %>%
  filter(is_penton | is_fiber | is_capsid_piiia) %>%
  pull(id_ref) %>% unique()

capsid_long <- proc_long %>% filter(id_ref %in% external_capsid_ids)
capsid_wide <- capsid_long %>%
  select(sample_id, id_ref, value) %>%
  pivot_wider(names_from = id_ref, values_from = value)

# Spearman correlation on features
cap_mat <- capsid_wide %>% select(-sample_id) %>% as.matrix()
cor_mat <- suppressWarnings(cor(cap_mat, use = "pairwise.complete.obs", method = "spearman"))

save_heatmap(cor_mat, file.path(FIG_DIR, "supp_heatmap_external_capsid_spearman.png"),
             main = "External capsid targets (Spearman correlation)")

# ------------------------
# Supplemental histograms of top antigens (PC2 loadings)
# ------------------------
top15_antigens <- loadings %>%
  mutate(abs_pc2 = abs(PC2)) %>%
  arrange(desc(abs_pc2)) %>%
  slice_head(n = 15) %>%
  pull(id_ref)

antigen_labels <- platform %>%
  filter(id_ref %in% top15_antigens) %>%
  select(id_ref, description) %>%
  deframe()

hist_long <- proc_long %>%
  filter(id_ref %in% top15_antigens) %>%
  select(id_ref, value) %>%
  mutate(id_ref = factor(id_ref, levels = top15_antigens)) %>%
  rename(antigen = id_ref)

p_hist <- ggplot(hist_long, aes(x = value)) +
  geom_histogram(binwidth = 0.2, fill = "steelblue", color = "black") +
  facet_wrap(~ antigen, labeller = labeller(antigen = antigen_labels), scales = "free_x", ncol = 4) +
  labs(x = "Antibody level", y = "Count") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(size = 10))

ggsave(file.path(FIG_DIR, "supp_hist_top15_antigens.png"), p_hist, width = 12, height = 9, dpi = 300)

# ------------------------
# Severity analyses (Ruuska score derived from TAC episode components)
# ------------------------
episodes_ruu <- episodes2 %>%
  mutate(
    ruu_dur_diarr = suppressWarnings(as.numeric(duradiar)),
    ruu_max_diarr = suppressWarnings(as.numeric(maxdiar)),
    ruu_dur_vomit = suppressWarnings(as.numeric(duravomi)),
    ruu_max_vomit = suppressWarnings(as.numeric(maxvomi)),
    ruu_fever_raw = suppressWarnings(as.numeric(temp)),
    ruu_dehyd_raw = suppressWarnings(as.numeric(dehyd)),
    ruu_treat_raw = suppressWarnings(as.numeric(rehyd))
  ) %>%
  mutate(
    ruu_dur_vomit = if_else(is.na(ruu_dur_vomit) | ruu_dur_vomit == 9, 0, ruu_dur_vomit),
    ruu_max_vomit = if_else(is.na(ruu_max_vomit) | ruu_max_vomit == 9, 0, ruu_max_vomit),
    ruu_fever = if_else(ruu_fever_raw == 4, NA_real_, ruu_fever_raw),
    ruu_dehyd = case_when(
      ruu_dehyd_raw == 1 ~ 0,
      ruu_dehyd_raw == 2 ~ 2,
      TRUE ~ NA_real_
    ),
    ruu_treat = case_when(
      ruu_treat_raw == 1 ~ 0,
      ruu_treat_raw == 2 ~ 1,
      ruu_treat_raw == 3 ~ 2,
      TRUE ~ NA_real_
    )
  ) %>%
  mutate(
    ruuska_score = {
      components <- c("ruu_dur_diarr", "ruu_max_diarr", "ruu_dur_vomit", "ruu_max_vomit", "ruu_fever", "ruu_dehyd", "ruu_treat")
      score <- rowSums(across(all_of(components)), na.rm = TRUE)
      if_else(rowSums(!is.na(across(all_of(components)))) == 0, NA_real_, score)
    },
    ruuska_cat = case_when(
      is.na(ruuska_score) ~ NA_character_,
      ruuska_score <= 6 ~ "Mild (0-6)",
      ruuska_score <= 10 ~ "Moderate (7-10)",
      TRUE ~ "Severe (>=11)"
    )
  )

adv_epi_selected <- episodes_ruu %>%
  filter(pcr_pos, year %in% c(1L, 2L)) %>%
  group_by(sid, year) %>%
  slice_min(adv_ct, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  left_join(samples %>% select(sample_id, infection_group, year1_infection, year2_infection), by = "sample_id")

year2_severity <- adv_epi_selected %>%
  filter(year == 2L, year2_infection == 1L) %>%
  select(sample_id, sid, diarr_specage, adv_ct, adv_afe, ruuska_score, ruuska_cat, infection_group)

write_csv(year2_severity, file.path(TAB_DIR, "year2_adv_selected_episode_severity.csv"))

target_map <- platform %>%
  mutate(desc_up = toupper(description)) %>%
  mutate(
    target_key = case_when(
      str_detect(desc_up, "HADV-40") & str_detect(desc_up, "PENTON BASE") ~ "AdV40_PB",
      str_detect(desc_up, "HADV-41") & str_detect(desc_up, "PENTON BASE") ~ "AdV41_PB",
      str_detect(desc_up, "HADV-40") & (str_detect(desc_up, "FIBER-2") | str_detect(desc_up, "SHORT FIBER")) ~ "AdV40_short_fiber",
      str_detect(desc_up, "HADV-41") & str_detect(desc_up, "SHORT FIBER") ~ "AdV41_short_fiber",
      str_detect(desc_up, "HADV-41") & str_detect(desc_up, "LONG FIBER") ~ "AdV41_long_fiber",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(target_key)) %>%
  select(id_ref, target_key)

top_targets_year2 <- proc_long %>%
  inner_join(target_map, by = "id_ref") %>%
  group_by(sample_id, target_key) %>%
  summarize(si_median = median(value, na.rm = TRUE), .groups = "drop") %>%
  inner_join(year2_severity, by = "sample_id")

safe_kw <- function(df) {
  if (n_distinct(df$ruuska_cat, na.rm = TRUE) < 2) return(NA_real_)
  suppressWarnings(kruskal.test(si_median ~ ruuska_cat, data = df)$p.value)
}

kw_severity <- top_targets_year2 %>%
  group_by(target_key) %>%
  summarize(p_kw = safe_kw(cur_data()), .groups = "drop")

write_csv(top_targets_year2, file.path(TAB_DIR, "top_targets_year2_severity_long.csv"))
write_csv(kw_severity, file.path(TAB_DIR, "top_targets_year2_severity_kw.csv"))

# ------------------------
# Save session info for reproducibility
# ------------------------
sink(file.path(OUT_DIR, "R_sessionInfo.txt"))
print(sessionInfo())
sink()

message("Done. Figures written to: ", FIG_DIR, " | Tables written to: ", TAB_DIR)
