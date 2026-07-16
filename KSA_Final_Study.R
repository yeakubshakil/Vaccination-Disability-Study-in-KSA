# ==============================================================================
# Vaccine-Preventable Childhood Disability & Immunization System Resilience
# Saudi Arabia: WHO/UNICEF WUENIC Time-Series Analysis, 1980-2024
# ==============================================================================
# Author(s): Md. Yeakub Ali
# Date: 2026-07-07
# Manuscript: "Vaccine-Preventable Childhood Disability and Immunization System 
#              Resilience in Saudi Arabia: A 45-Year WHO/UNICEF Time-Series Study"
# ==============================================================================

# ---- 0. Setup ---------------------------------------------------------------

RAW_DATA_FILE <- "KSA Immune Study.csv"
OUTPUT_DIR    <- "outputs"

pkgs <- c("readr","dplyr","stringr","tidyr","janitor","purrr",
          "ggplot2","scales","forcats","tibble","flextable","officer","grid")

invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
  library(p, character.only = TRUE)
}))

for (d in c("data","tables","docx","figures")) dir.create(file.path(OUTPUT_DIR, d), showWarnings = FALSE, recursive = TRUE)

# ---- Helper functions --------------------------------------------------------

format_p      <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
signed_pp     <- function(x, d = 0) ifelse(is.na(x), "NA", ifelse(x > 0, paste0("+", sprintf(paste0("%.", d, "f"), x)), sprintf(paste0("%.", d, "f"), x)))
cap01         <- function(x) pmin(pmax(x, 0), 1)
cap100        <- function(x) pmin(pmax(x, 0), 100)
pop_sd        <- function(x) { x <- x[!is.na(x)]; if (length(x) == 0) return(NA_real_); sqrt(mean((x - mean(x))^2)) }
extract_code  <- function(x) str_extract(x, "^[^:]+")
extract_label <- function(x) str_squish(str_remove(x, "^[^:]+:\s*"))
safe_value    <- function(df, yr) { v <- df %>% filter(year == yr) %>% pull(coverage); if (length(v) == 1) v else NA_real_ }

safe_t_test <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  tryCatch(if (length(x) >= 2 && length(y) >= 2 && length(unique(c(x, y))) > 1) t.test(x, y)$p.value else NA_real_, error = function(e) NA_real_)
}

write_tbl <- function(df, fn) write_csv(df, file.path(OUTPUT_DIR, "tables", fn), na = "NA")

make_ft <- function(df, note = NULL, fs = 8) {
  ft <- flextable(df) %>% fontsize(size = fs, part = "all") %>% bold(part = "header") %>%
    align(align = "center", part = "all") %>% align(j = 1, align = "left", part = "body") %>% autofit()
  if (!is.null(note)) ft <- ft %>% add_footer_lines(note) %>% fontsize(size = fs - 1, part = "footer") %>% align(align = "left", part = "footer")
  ft
}

add_to_doc <- function(doc, title, df, note = NULL, fs = 8) {
  doc <- body_add_par(doc, title, style = "heading 2")
  doc <- body_add_flextable(doc, make_ft(df, note = note, fs = fs))
  body_add_par(doc, "", style = "Normal")
}

theme_pub <- function(bs = 13) {
  theme_minimal(base_size = bs) +
    theme(plot.title = element_text(face = "bold", size = bs + 4, color = "black"),
          plot.subtitle = element_text(size = bs, color = "black"),
          axis.title = element_text(face = "bold", color = "black"),
          axis.text = element_text(color = "#333"),
          legend.title = element_text(face = "bold", color = "black"),
          panel.grid.major = element_line(color = "#D9D9D9", linewidth = 0.45),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          legend.background = element_rect(fill = "white", color = NA),
          strip.background = element_rect(fill = "#F5F5F5", color = NA),
          strip.text = element_text(face = "bold", color = "black"))
}

# Color palette
cols <- list(blue = "#2C7FB8", orange = "#F28E2B", green = "#1B9E77", red = "#D95F02",
             purple = "#7570B3", pink = "#E7298A", olive = "#66A61E", gold = "#E6AB02",
             teal = "#009E73", grey = "#7A7A7A", dark = "#333333", lightblue = "#A6CEE3")

# ==============================================================================
# 1. DATA IMPORT & PREPROCESSING
# ==============================================================================

if (!file.exists(RAW_DATA_FILE)) stop("Raw data file not found: ", RAW_DATA_FILE)

raw_data <- read_csv(RAW_DATA_FILE, show_col_types = FALSE)
raw_clean <- clean_names(raw_data)

req_cols <- c("ref_area_geographic_area","unit_measure_unit_of_measure",
              "obs_status_observation_status","indicator_indicator",
              "vaccine_vaccine","age_current_age","time_period_time_period",
              "obs_value_observation_value")

if (length(missing <- setdiff(req_cols, names(raw_clean))) > 0)
  stop("Missing columns: ", paste(missing, collapse = ", "))

ksa_data <- raw_clean %>%
  filter(str_detect(ref_area_geographic_area, "SAU: Saudi Arabia"),
         unit_measure_unit_of_measure == "PCNT: %",
         obs_status_observation_status == "E: Estimated value") %>%
  transmute(
    country = ref_area_geographic_area,
    indicator_code = extract_code(indicator_indicator),
    indicator_name = extract_label(indicator_indicator),
    vaccine_code = extract_code(vaccine_vaccine),
    vaccine_name = extract_label(vaccine_vaccine),
    age_group = extract_label(age_current_age),
    year = as.integer(time_period_time_period),
    coverage = as.numeric(obs_value_observation_value)
  ) %>% arrange(vaccine_code, year)

write_csv(ksa_data, file.path(OUTPUT_DIR, "data", "ksa_immunization_cleaned.csv"), na = "NA")

# ---- Dataset screening ----

dataset_screening <- tibble(
  Item = c("Raw rows","Raw columns","Geographic areas","Saudi Arabia rows",
           "Vaccine indicators","Year range","Missing coverage values","Duplicate rows"),
  Result = c(nrow(raw_data), ncol(raw_data), n_distinct(raw_clean$ref_area_geographic_area),
             nrow(ksa_data), n_distinct(ksa_data$vaccine_code),
             paste0(min(ksa_data$year, na.rm = TRUE), "-", max(ksa_data$year, na.rm = TRUE)),
             sum(is.na(ksa_data$coverage)),
             sum(duplicated(ksa_data[, c("vaccine_code", "year")])))
)
write_tbl(dataset_screening, "Supplementary_Table_S1_Dataset_Screening.csv")

# ---- Vaccine availability summary ----

availability_summary <- ksa_data %>%
  group_by(vaccine_code, vaccine_name) %>%
  summarise(`Age group` = paste(sort(unique(age_group)), collapse = "; "),
            n = n(), first_year = min(year, na.rm = TRUE), last_year = max(year, na.rm = TRUE),
            first_coverage = coverage[which.min(year)], latest_coverage = coverage[which.max(year)],
            mean_coverage = round(mean(coverage, na.rm = TRUE), 1),
            median_coverage = round(median(coverage, na.rm = TRUE), 1),
            min_coverage = min(coverage, na.rm = TRUE), max_coverage = max(coverage, na.rm = TRUE),
            .groups = "drop") %>% arrange(first_year, vaccine_code)
write_tbl(availability_summary, "Supplementary_Table_S2_Vaccine_Indicator_Availability.csv")

# ==============================================================================
# 2. VACCINE-DISABILITY CLASSIFICATION
# ==============================================================================

classification_lookup <- tribble(
  ~vaccine_code, ~disability_category, ~main_relevance,
  "POL3",  "Direct physical disability prevention",       "Prevention of paralytic poliomyelitis and lifelong motor disability",
  "IPV1",  "Direct physical disability prevention",       "Additional polio/paralysis prevention indicator",
  "IPV2",  "Direct physical disability prevention",       "Additional polio/paralysis prevention indicator; limited observations",
  "MCV1",  "Neurological/sensory disability prevention",  "Prevention of measles-related encephalitis, blindness, hearing-related complications",
  "MCV2",  "Neurological/sensory disability prevention",  "Second-dose measles protection and sustained measles-related disability prevention",
  "RCV1",  "Congenital/developmental disability prevention","Prevention of congenital rubella syndrome-related hearing, visual, cardiac, and developmental impairment",
  "HIB3",  "Meningitis-related disability prevention",      "Prevention of Hib meningitis-related deafness, brain damage, and neurological sequelae",
  "PCV3",  "Meningitis-related disability prevention",      "Prevention of pneumococcal meningitis-related hearing loss and developmental delay",
  "DTP1",  "Severe childhood morbidity prevention",       "Early access indicator for prevention of severe vaccine-preventable childhood disease",
  "DTP3",  "Severe childhood morbidity prevention",       "Completion indicator for prevention of diphtheria, tetanus, and pertussis-related severe morbidity",
  "BCG",   "Severe childhood morbidity prevention",       "Prevention of severe childhood tuberculosis forms and severe infectious morbidity",
  "HEPB3", "Chronic morbidity prevention",                "Prevention of chronic hepatitis B infection and long-term liver disease burden",
  "HEPBB", "Chronic morbidity prevention",                "Early-life prevention of hepatitis B transmission and future chronic morbidity",
  "ROTAC", "Severe childhood morbidity prevention",       "Prevention of severe rotavirus gastroenteritis and hospitalization burden",
  "HPV",   "Emerging adolescent morbidity-prevention gap","Prevention of HPV-related future cancer morbidity; insufficient observations for trend modelling"
)

vaccine_order <- c("POL3","IPV1","IPV2","MCV1","MCV2","RCV1","HIB3","PCV3",
                   "DTP1","DTP3","BCG","HEPB3","HEPBB","ROTAC","HPV")

table1_classification <- availability_summary %>%
  left_join(classification_lookup, by = "vaccine_code") %>%
  transmute(Vaccine = vaccine_code, `Vaccine name` = vaccine_name, `Age group`,
            `Dataset period` = paste0(first_year, "-", last_year), n,
            `Disability-prevention category` = disability_category,
            `Main disability/morbidity relevance` = main_relevance) %>%
  arrange(match(Vaccine, vaccine_order))
write_tbl(table1_classification, "Table_1_Vaccine_Disability_Classification.csv")

# ==============================================================================
# 3. LONG-TERM TREND ANALYSIS
# ==============================================================================

calc_trend <- function(df) {
  df <- df %>% filter(!is.na(year), !is.na(coverage)) %>% arrange(year)
  n <- nrow(df); fy <- min(df$year); ly <- max(df$year)
  fc <- df$coverage[which.min(df$year)]; lc <- df$coverage[which.max(df$year)]
  mc <- mean(df$coverage, na.rm = TRUE); sdc <- if (n <= 1) 0 else sd(df$coverage, na.rm = TRUE)
  
  slope <- ols_p <- mk_p <- NA_real_
  if (n >= 2 && length(unique(df$year)) >= 2) {
    fit <- tryCatch(lm(coverage ~ year, data = df), error = function(e) NULL)
    if (!is.null(fit)) { ct <- summary(fit)$coefficients; if ("year" %in% rownames(ct)) { slope <- unname(ct["year","Estimate"]); ols_p <- unname(ct["year","Pr(>|t|)"]) } }
  }
  if (n >= 4 && length(unique(df$coverage)) > 1) {
    mk_p <- tryCatch(cor.test(df$year, df$coverage, method = "kendall", exact = FALSE)$p.value, error = function(e) NA_real_)
  }
  
  interp <- case_when(
    n == 1 ~ "Scenario analysis only", n < 5 ~ "Descriptive only; stable",
    !is.na(ols_p) && !is.na(mk_p) && ols_p < 0.05 && mk_p < 0.05 && slope > 0 ~ "Significant increase",
    !is.na(ols_p) && !is.na(mk_p) && ols_p < 0.05 && mk_p < 0.05 && slope < 0 ~ "Significant slight decrease",
    lc >= 95 ~ "Stable high coverage", TRUE ~ "No significant trend"
  )
  
  tibble(N = n, first_year = fy, last_year = ly, mean_coverage = mc, sd_coverage = sdc,
         first_coverage = fc, latest_coverage = lc, change_pp = lc - fc,
         slope_pp_year = slope, ols_p = ols_p, mk_p = mk_p, trend_interpretation = interp)
}

trend_raw <- ksa_data %>% group_by(vaccine_code) %>% group_modify(~ calc_trend(.x)) %>% ungroup()

table2_order <- c("BCG","DTP1","DTP3","POL3","MCV1","HEPB3","HEPBB","HIB3","MCV2","PCV3","ROTAC","RCV1","IPV1","IPV2","HPV")

table2_trends <- trend_raw %>%
  transmute(Vaccine = vaccine_code, N, Period = paste0(first_year, "-", last_year),
            `Mean ± SD` = paste0(sprintf("%.1f", mean_coverage), " ± ", sprintf("%.1f", sd_coverage)),
            `First %` = sprintf("%.0f", first_coverage), `Latest %` = sprintf("%.0f", latest_coverage),
            `Change pp` = signed_pp(change_pp, 0),
            `Slope pp/year` = ifelse(is.na(slope_pp_year), "NA", signed_pp(slope_pp_year, 2)),
            `OLS p` = format_p(ols_p), `MK p` = format_p(mk_p), `Trend interpretation` = trend_interpretation) %>%
  arrange(match(Vaccine, table2_order))
write_tbl(table2_trends, "Table_2_Descriptive_Long_Term_Trends.csv")

# ==============================================================================
# 4. ERA-WISE COVERAGE ANALYSIS
# ==============================================================================

era_levels <- c("Foundation (1980-1989)", "Expansion (1990-1999)", "Consolidation (2000-2009)",
                "Optimization (2010-2019)", "Resilience (2020-2024)")

era_data <- ksa_data %>%
  mutate(Era = case_when(
    year >= 1980 & year <= 1989 ~ era_levels[1], year >= 1990 & year <= 1999 ~ era_levels[2],
    year >= 2000 & year <= 2009 ~ era_levels[3], year >= 2010 & year <= 2019 ~ era_levels[4],
    year >= 2020 & year <= 2024 ~ era_levels[5], TRUE ~ NA_character_
  )) %>% filter(!is.na(Era))

era_mean <- era_data %>% mutate(Era = factor(Era, levels = era_levels)) %>%
  group_by(vaccine_code, Era) %>% summarise(mean_coverage = round(mean(coverage, na.rm = TRUE), 1), n = n(), .groups = "drop")

era_p <- era_data %>% group_by(vaccine_code) %>%
  summarise(`Kruskal-Wallis p` = tryCatch(if (n_distinct(Era) >= 2 && n_distinct(coverage) > 1) kruskal.test(coverage ~ Era)$p.value else NA_real_, error = function(e) NA_real_), .groups = "drop") %>%
  mutate(`Kruskal-Wallis p` = format_p(`Kruskal-Wallis p`))

table_s3 <- era_mean %>% select(vaccine_code, Era, mean_coverage) %>%
  pivot_wider(names_from = Era, values_from = mean_coverage) %>%
  left_join(era_p, by = "vaccine_code") %>% rename(Vaccine = vaccine_code)
write_tbl(table_s3, "Supplementary_Table_S3_Era_Wise_Coverage.csv")

table_s4 <- era_mean %>% select(vaccine_code, Era, n) %>%
  pivot_wider(names_from = Era, values_from = n, values_fill = 0) %>% rename(Vaccine = vaccine_code)
write_tbl(table_s4, "Supplementary_Table_S4_Era_Observation_Counts.csv")

# ==============================================================================
# 5. COVID-19 SHOCK, RECOVERY & VPDP-RI
# ==============================================================================

covid_metrics <- function(df) {
  df <- df %>% arrange(year)
  baseline_v <- df %>% filter(year >= 2017, year <= 2019) %>% pull(coverage)
  recovery_v <- df %>% filter(year >= 2022, year <= 2024) %>% pull(coverage)
  baseline <- if (length(baseline_v) > 0) mean(baseline_v, na.rm = TRUE) else NA_real_
  v2019 <- safe_value(df, 2019); v2020 <- safe_value(df, 2020); v2024 <- safe_value(df, 2024)
  
  its_p <- tryCatch({
    its_df <- df %>% filter(year >= 2015, year <= 2024) %>%
      mutate(time = year - min(year, na.rm = TRUE), post_covid = as.integer(year >= 2020),
             time_after_covid = ifelse(year >= 2020, year - 2019, 0))
    if (nrow(its_df) >= 6 && length(unique(its_df$coverage)) > 1) {
      ct <- summary(lm(coverage ~ time + post_covid + time_after_covid, data = its_df))$coefficients
      if ("post_covid" %in% rownames(ct)) ct["post_covid","Pr(>|t|)"] else NA_real_
    } else NA_real_
  }, error = function(e) NA_real_)
  
  stability <- if (length(baseline_v) >= 2) cap01(1 - pop_sd(baseline_v) / 100) else NA_real_
  shock_r <- cap01(1 - max(0, baseline - v2020, na.rm = TRUE) / 100)
  recovery_c <- cap01(1 - max(0, baseline - v2024, na.rm = TRUE) / 100)
  sustain <- case_when(is.na(v2024) ~ NA_real_, v2024 >= 95 ~ 1, v2024 >= 90 ~ 0.5, TRUE ~ 0)
  vpdp_ri <- mean(c(stability, shock_r, recovery_c, sustain), na.rm = TRUE)
  
  tibble(baseline_2017_2019 = baseline, coverage_2019 = v2019, coverage_2020 = v2020, coverage_2024 = v2024,
         shock_pp = v2020 - v2019, recovery_pp = v2024 - v2020, net_change_pp = v2024 - v2019,
         pre_vs_recovery_p = safe_t_test(baseline_v, recovery_v), exploratory_its_p = its_p,
         stability = stability, shock_resistance = shock_r, recovery_component = recovery_c,
         sustainability = sustain, vpdp_ri = vpdp_ri)
}

covid_raw <- ksa_data %>% filter(year >= 2017, year <= 2024) %>% group_by(vaccine_code) %>% group_modify(~ covid_metrics(.x)) %>% ungroup()

table3a_order <- c("DTP1","DTP3","POL3","MCV1","MCV2","HEPB3","HEPBB","HIB3","PCV3","ROTAC","RCV1","IPV1")

table3a_main <- covid_raw %>% filter(!(vaccine_code %in% c("BCG","HPV","IPV2"))) %>%
  transmute(Vaccine = vaccine_code, `Baseline 2017-2019` = sprintf("%.1f", baseline_2017_2019),
            `2020 %` = sprintf("%.0f", coverage_2020), `2024 %` = sprintf("%.0f", coverage_2024),
            `Shock pp` = signed_pp(shock_pp, 0), `Recovery pp` = signed_pp(recovery_pp, 0),
            `Net change pp` = signed_pp(net_change_pp, 0), `Pre vs recovery p` = format_p(pre_vs_recovery_p)) %>%
  arrange(match(Vaccine, table3a_order))
write_tbl(table3a_main, "Table_3A_COVID_Shock_Recovery_Main_BCG_Excluded.csv")

table3b_vpdp <- covid_raw %>% filter(!(vaccine_code %in% c("HPV","IPV2"))) %>% arrange(desc(vpdp_ri)) %>%
  mutate(Rank = row_number()) %>%
  transmute(Rank, Vaccine = vaccine_code, Stability = sprintf("%.3f", stability),
            `Shock resistance` = sprintf("%.3f", shock_resistance), Recovery = sprintf("%.3f", recovery_component),
            Sustainability = sprintf("%.3f", sustainability), `VPDP-RI` = sprintf("%.3f", vpdp_ri))
write_tbl(table3b_vpdp, "Table_3B_VPDP_RI_Resilience_Ranking.csv")

table_s5 <- covid_raw %>% filter(!(vaccine_code %in% c("HPV","IPV2"))) %>%
  transmute(Vaccine = vaccine_code, `Baseline 2017-2019` = sprintf("%.1f", baseline_2017_2019),
            `2019 %` = sprintf("%.0f", coverage_2019), `2020 %` = sprintf("%.0f", coverage_2020),
            `2024 %` = sprintf("%.0f", coverage_2024), `Shock pp` = signed_pp(shock_pp, 0),
            `Recovery pp` = signed_pp(recovery_pp, 0), `Net change pp` = signed_pp(net_change_pp, 0),
            `Pre vs recovery p` = format_p(pre_vs_recovery_p), `Exploratory ITS p` = format_p(exploratory_its_p))
write_tbl(table_s5, "Supplementary_Table_S5_COVID_Full_With_BCG_and_ITS.csv")

# ==============================================================================
# 6. GCC BENCHMARKING
# ==============================================================================

gcc_codes <- c("SAU","BHR","KWT","OMN","QAT","ARE")
core_vaccines <- c("DTP3","POL3","MCV1","MCV2","PCV3","ROTAC")
benchmark_year <- 2023

gcc_data <- raw_clean %>%
  mutate(country_code = extract_code(ref_area_geographic_area), country = extract_label(ref_area_geographic_area),
         vaccine_code = extract_code(vaccine_vaccine), year = as.integer(time_period_time_period),
         coverage = as.numeric(obs_value_observation_value)) %>%
  filter(country_code %in% gcc_codes, unit_measure_unit_of_measure == "PCNT: %",
         obs_status_observation_status == "E: Estimated value") %>%
  select(country_code, country, vaccine_code, vaccine_name, year, coverage)
write_csv(gcc_data, file.path(OUTPUT_DIR, "data", "gcc_immunization_cleaned.csv"), na = "NA")

table4_gcc <- gcc_data %>% filter(year == benchmark_year, vaccine_code %in% core_vaccines) %>%
  select(country, vaccine_code, coverage) %>% pivot_wider(names_from = vaccine_code, values_from = coverage) %>%
  mutate(`Mean VPDP core` = round(rowMeans(across(all_of(core_vaccines)), na.rm = TRUE), 1),
         `Available indicators` = rowSums(!is.na(across(all_of(core_vaccines))))) %>%
  arrange(desc(`Mean VPDP core`)) %>% mutate(Rank = row_number()) %>%
  select(Rank, Country = country, all_of(core_vaccines), `Mean VPDP core`, `Available indicators`)
write_tbl(table4_gcc, "Table_4_GCC_Benchmarking_2023.csv")

table_s6 <- gcc_data %>% filter(year == benchmark_year, vaccine_code %in% core_vaccines) %>% group_by(vaccine_code) %>%
  summarise(`Saudi Arabia %` = coverage[country == "Saudi Arabia"][1],
            `Other GCC mean %` = mean(coverage[country != "Saudi Arabia"], na.rm = TRUE),
            `Difference pp` = `Saudi Arabia %` - `Other GCC mean %`,
            `Exploratory p` = tryCatch({
              sau <- coverage[country == "Saudi Arabia"]; others <- coverage[country != "Saudi Arabia"]
              if (length(sau) == 1 && length(others[!is.na(others)]) >= 2) t.test(others, mu = sau)$p.value else NA_real_
            }, error = function(e) NA_real_), .groups = "drop") %>%
  transmute(Vaccine = vaccine_code, `Saudi Arabia %` = sprintf("%.0f", `Saudi Arabia %`),
            `Other GCC mean %` = sprintf("%.1f", `Other GCC mean %`),
            `Difference pp` = signed_pp(`Difference pp`, 1), `Exploratory p` = format_p(`Exploratory p`))
write_tbl(table_s6, "Supplementary_Table_S6_GCC_Vaccine_Specific_Comparisons.csv")

ksa_core <- gcc_data %>% filter(year == benchmark_year, country == "Saudi Arabia", vaccine_code %in% core_vaccines) %>% pull(coverage)
other_core <- gcc_data %>% filter(year == benchmark_year, country != "Saudi Arabia", vaccine_code %in% core_vaccines) %>% pull(coverage)

table_s7 <- tibble(Comparison = "Saudi Arabia vs other GCC pooled core indicators",
                   `Saudi Arabia mean %` = sprintf("%.1f", mean(ksa_core, na.rm = TRUE)),
                   `Other GCC mean %` = sprintf("%.1f", mean(other_core, na.rm = TRUE)),
                   `Difference pp` = signed_pp(mean(ksa_core, na.rm = TRUE) - mean(other_core, na.rm = TRUE), 1),
                   `Exploratory p` = format_p(safe_t_test(ksa_core, other_core)))
write_tbl(table_s7, "Supplementary_Table_S7_GCC_Pooled_Comparison.csv")

# ==============================================================================
# 7. 2030 READINESS & HPV SCENARIO ANALYSIS
# ==============================================================================

tracer_vaccines <- c("DTP3","MCV2","PCV3","HPV")
readiness_order <- c("DTP3","MCV2","PCV3","HPV","POL3","IPV1","IPV2","MCV1","RCV1","HIB3","BCG","HEPB3","HEPBB","ROTAC")

latest_by_vaccine <- ksa_data %>% group_by(vaccine_code) %>% slice_max(order_by = year, n = 1, with_ties = FALSE) %>% ungroup()

table5a_readiness <- latest_by_vaccine %>% filter(vaccine_code %in% readiness_order) %>%
  transmute(Vaccine = vaccine_code, `Latest year` = as.character(year), `Latest coverage %` = sprintf("%.0f", coverage),
            `IA2030 tracer indicator` = ifelse(vaccine_code %in% tracer_vaccines, "Yes", "No"),
            `IA2030 90% status` = case_when(vaccine_code %in% tracer_vaccines & coverage >= 90 ~ "Met",
                                            vaccine_code %in% tracer_vaccines & coverage < 90 ~ "Not met", TRUE ~ "Not applicable"),
            `95% maintenance status` = ifelse(coverage >= 95, "Met", "Not met"),
            `Readiness interpretation` = case_when(
              vaccine_code == "HPV" ~ "Major 2030 readiness gap",
              vaccine_code == "PCV3" ~ "On track, but monitor slight decline",
              vaccine_code %in% tracer_vaccines ~ "On track; maintain high coverage",
              vaccine_code == "IPV2" ~ "Descriptive only; recent indicator",
              TRUE ~ "Strong disability-prevention readiness")) %>%
  arrange(match(Vaccine, readiness_order))
write_tbl(table5a_readiness, "Table_5A_Latest_Coverage_and_2030_Readiness.csv")

recent_proj <- ksa_data %>% filter(year >= 2015, year <= 2024, vaccine_code != "HPV") %>% group_by(vaccine_code) %>%
  group_modify(~ {
    df <- .x; ly <- max(df$year, na.rm = TRUE); lc <- df$coverage[which.max(df$year)]
    slope <- tryCatch(if (nrow(df) >= 4 && length(unique(df$year)) >= 2) unname(coef(lm(coverage ~ year, data = df))[["year"]]) else NA_real_, error = function(e) NA_real_)
    proj <- if (is.na(slope)) NA_real_ else cap100(lc + slope * (2030 - ly))
    tibble(latest_year = ly, latest_coverage = lc, recent_slope_pp_year = slope, projected_2030 = proj,
           projection_interpretation = case_when(is.na(proj) ~ "Insufficient recent observations",
                                                 proj >= 95 ~ "Likely high-coverage maintenance",
                                                 proj >= 90 ~ "Likely IA2030 threshold maintenance", TRUE ~ "Potential risk; monitor"))
  }) %>% ungroup() %>%
  transmute(Vaccine = vaccine_code, `Latest year` = as.character(latest_year), `Latest coverage %` = sprintf("%.0f", latest_coverage),
            `Recent slope pp/year` = ifelse(is.na(recent_slope_pp_year), "NA", signed_pp(recent_slope_pp_year, 2)),
            `Projected 2030 %` = ifelse(is.na(projected_2030), "NA", sprintf("%.1f", projected_2030)), `Projection interpretation`)
write_tbl(recent_proj, "Supplementary_Table_S8_Recent_Trend_2030_Projection.csv")

# ---- HPV scenarios ----

hpv_latest <- latest_by_vaccine %>% filter(vaccine_code == "HPV") %>% slice(1)
hpv_start_year <- hpv_latest$year; hpv_start_coverage <- hpv_latest$coverage
required_annual <- (90 - hpv_start_coverage) / (2030 - hpv_start_year)

hpv_scenarios <- tribble(
  ~Scenario, ~annual_increase,
  "Status quo slow scale-up", 2.0,
  "Moderate scale-up", 4.0,
  "Required scale-up", required_annual,
  "Accelerated scale-up", 8.0
)

table5b_hpv <- hpv_scenarios %>%
  mutate(`Annual increase pp/year` = round(annual_increase, 1), `2023 baseline %` = hpv_start_coverage,
         `Projected 2030 %` = cap100(hpv_start_coverage + annual_increase * (2030 - hpv_start_year)),
         `Gap to 90% target` = pmax(0, 90 - `Projected 2030 %`),
         `Target status` = ifelse(`Projected 2030 %` >= 90, "Met", "Not met"),
         `Policy meaning` = case_when(Scenario == "Status quo slow scale-up" ~ "Insufficient acceleration",
                                      Scenario == "Moderate scale-up" ~ "Improved but still inadequate",
                                      Scenario == "Required scale-up" ~ "Minimum pathway to the 2030 target",
                                      Scenario == "Accelerated scale-up" ~ "Ambitious school/adolescent strategy")) %>%
  transmute(Scenario, `Annual increase pp/year` = sprintf("%.1f", `Annual increase pp/year`),
            `2023 baseline %` = sprintf("%.0f", `2023 baseline %`), `Projected 2030 %` = sprintf("%.0f", `Projected 2030 %`),
            `Gap to 90% target` = sprintf("%.0f", `Gap to 90% target`), `Target status`, `Policy meaning`)
write_tbl(table5b_hpv, "Table_5B_HPV_Scale_Up_Scenarios.csv")

hpv_yearly <- hpv_scenarios %>% crossing(year = hpv_start_year:2030) %>%
  mutate(coverage = cap100(hpv_start_coverage + annual_increase * (year - hpv_start_year)),
         scenario = factor(Scenario, levels = c("Accelerated scale-up","Moderate scale-up","Required scale-up","Status quo slow scale-up")))
write_tbl(hpv_yearly, "Supplementary_Table_S9_HPV_Yearly_Scenarios.csv")

table_s10 <- table5a_readiness %>% count(`IA2030 90% status`, `95% maintenance status`, name = "n")
write_tbl(table_s10, "Supplementary_Table_S10_Readiness_Status_Summary.csv")

# ==============================================================================
# 8. BCG 2019 OUTLIER SENSITIVITY ANALYSIS
# ==============================================================================

bcg_raw <- ksa_data %>% filter(vaccine_code == "BCG") %>% arrange(year)
bcg_2018 <- safe_value(bcg_raw, 2018); bcg_2019 <- safe_value(bcg_raw, 2019); bcg_2020 <- safe_value(bcg_raw, 2020)
bcg_interp_2019 <- mean(c(bcg_2018, bcg_2020), na.rm = TRUE)

diag_table <- tibble(Indicator = c("BCG coverage 2018","Raw BCG 2019","BCG coverage 2020",
                                   "Linear interpolated 2019","Raw minus interpolated"),
                     Value = c(sprintf("%.1f", bcg_2018), sprintf("%.1f", bcg_2019), sprintf("%.1f", bcg_2020),
                               sprintf("%.1f", bcg_interp_2019), signed_pp(bcg_2019 - bcg_interp_2019, 1)))
write_tbl(diag_table, "Supplementary_Table_S11A_BCG_Outlier_Diagnostic.csv")

bcg_scenarios <- bind_rows(
  bcg_raw %>% mutate(Scenario = "Raw 2019 value retained"),
  bcg_raw %>% filter(year != 2019) %>% mutate(Scenario = "2019 value excluded"),
  bcg_raw %>% mutate(coverage = ifelse(year == 2019, bcg_interp_2019, coverage), Scenario = "2019 value interpolated")
) %>% select(Scenario, everything())

bcg_trend_sens <- bcg_scenarios %>% group_by(Scenario) %>% group_modify(~ calc_trend(.x)) %>% ungroup() %>%
  transmute(Scenario, N, Period = paste0(first_year, "-", last_year), `First %` = sprintf("%.0f", first_coverage),
            `Latest %` = sprintf("%.0f", latest_coverage), `Change pp` = signed_pp(change_pp, 0),
            `Slope pp/year` = signed_pp(slope_pp_year, 2), `OLS p` = format_p(ols_p), `MK p` = format_p(mk_p), `Trend interpretation`)
write_tbl(bcg_trend_sens, "Supplementary_Table_S11B_BCG_Long_Term_Trend_Sensitivity.csv")

bcg_era_sens <- bcg_scenarios %>%
  mutate(Era = case_when(year >= 1980 & year <= 1989 ~ era_levels[1], year >= 1990 & year <= 1999 ~ era_levels[2],
                         year >= 2000 & year <= 2009 ~ era_levels[3], year >= 2010 & year <= 2019 ~ era_levels[4],
                         year >= 2020 & year <= 2024 ~ era_levels[5], TRUE ~ NA_character_)) %>%
  filter(!is.na(Era)) %>% group_by(Scenario, Era) %>%
  summarise(`Mean BCG coverage` = round(mean(coverage, na.rm = TRUE), 1), .groups = "drop") %>%
  pivot_wider(names_from = Era, values_from = `Mean BCG coverage`)
write_tbl(bcg_era_sens, "Supplementary_Table_S11C_BCG_Era_Sensitivity.csv")

bcg_covid_sens <- bcg_scenarios %>% group_by(Scenario) %>% group_modify(~ covid_metrics(.x)) %>% ungroup() %>%
  transmute(Scenario, `Baseline 2017-2019` = sprintf("%.1f", baseline_2017_2019),
            `2019 %` = ifelse(is.na(coverage_2019), "NA", sprintf("%.1f", coverage_2019)),
            `2020 %` = sprintf("%.0f", coverage_2020), `2024 %` = sprintf("%.0f", coverage_2024),
            `Shock pp` = signed_pp(shock_pp, 1), `Recovery pp` = signed_pp(recovery_pp, 1),
            `Net change pp` = signed_pp(net_change_pp, 1), `Pre vs recovery p` = format_p(pre_vs_recovery_p),
            `VPDP-RI` = sprintf("%.3f", vpdp_ri))
write_tbl(bcg_covid_sens, "Supplementary_Table_S11D_BCG_COVID_VPDP_Sensitivity.csv")

# ==============================================================================
# 9. FIGURES
# ==============================================================================

# ---- Figure 1: Conceptual framework ----

nodes <- tribble(
  ~id, ~xmin, ~xmax, ~ymin, ~ymax, ~label, ~fill,
  "i1", 0.05, 0.20, 0.68, 0.78, "Health System\nInputs", cols$blue,
  "i2", 0.05, 0.20, 0.50, 0.60, "Policy &\nGovernance", cols$blue,
  "i3", 0.05, 0.20, 0.32, 0.42, "Community\nEngagement", cols$lightblue,
  "c",  0.33, 0.57, 0.42, 0.66, "CHILDHOOD\nIMMUNIZATION\nCOVERAGE\n(Proxy)", "#4C9ED9",
  "o1", 0.70, 0.84, 0.68, 0.78, "Child Health\nOutcomes", "#4DBD7A",
  "o2", 0.70, 0.84, 0.50, 0.60, "Disease\nElimination", "#4DBD7A",
  "o3", 0.70, 0.84, 0.32, 0.42, "Health System\nResilience", "#3CB7A0",
  "b1", 0.08, 0.23, 0.05, 0.15, "Long-term\nTrends", "#E96A5A",
  "b2", 0.26, 0.41, 0.05, 0.15, "COVID-19\nResilience", "#CC5A4A",
  "b3", 0.44, 0.59, 0.05, 0.15, "GCC\nBenchmarking", "#E99645",
  "b4", 0.62, 0.77, 0.05, 0.15, "2030\nReadiness", "#DE6F23",
  "b5", 0.80, 0.95, 0.05, 0.15, "VPDP-RI\nIndex", "#C45A4E"
)

arrows <- tribble(
  ~x, ~y, ~xend, ~yend, ~lt,
  0.20, 0.73, 0.33, 0.54, "solid", 0.20, 0.55, 0.33, 0.54, "solid",
  0.20, 0.37, 0.33, 0.54, "solid", 0.57, 0.54, 0.70, 0.73, "solid",
  0.57, 0.54, 0.70, 0.55, "solid", 0.57, 0.54, 0.70, 0.37, "solid",
  0.155, 0.15, 0.45, 0.42, "dashed", 0.335, 0.15, 0.45, 0.42, "dashed",
  0.515, 0.15, 0.45, 0.42, "dashed", 0.695, 0.15, 0.45, 0.42, "dashed",
  0.875, 0.15, 0.45, 0.42, "dashed"
)

fig1 <- ggplot() +
  geom_rect(data = nodes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), color = NA) +
  scale_fill_identity() +
  geom_segment(data = arrows, aes(x = x, y = y, xend = xend, yend = yend, linetype = lt),
               arrow = arrow(length = unit(0.18, "cm"), type = "open"), linewidth = 0.8, color = "grey45") +
  scale_linetype_identity() +
  geom_text(data = nodes, aes(x = (xmin+xmax)/2, y = (ymin+ymax)/2, label = label), color = "white", size = 4.5, fontface = "bold", lineheight = 1.05) +
  annotate("text", x = 0.5, y = 0.98, label = "Conceptual Framework", size = 8.5, fontface = "bold") +
  annotate("text", x = 0.5, y = 0.93, label = "Childhood immunization coverage as a proxy indicator of vaccine-preventable disability prevention",
           size = 5.3, fontface = "italic", color = "#4A4A4A") +
  annotate("text", x = 0.5, y = 0.005, label = "Saudi Arabia | 1980-2024 | WHO/UNICEF WUENIC", size = 4.4, color = "#7A7A7A") +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1), clip = "off") + theme_void()

ggsave(file.path(OUTPUT_DIR, "figures", "Figure1_Conceptual_Framework.png"), fig1, width = 14, height = 10, dpi = 600)

# ---- Figure 2: Long-term trends ----

fig2_vax <- c("DTP3","POL3","MCV1","MCV2","RCV1","HIB3","PCV3","ROTAC")
fig2_data <- ksa_data %>% filter(vaccine_code %in% fig2_vax) %>% mutate(vaccine_code = factor(vaccine_code, levels = fig2_vax))

era_shade <- tribble(~xmin, ~xmax, ~fill,
                     1980, 1989.99, "#EAF2FB", 1990, 1999.99, "#DDEBF7",
                     2000, 2009.99, "#E9F1F7", 2010, 2019.99, "#D7ECFA", 2020, 2024.99, "#C9E4F7")

fig2 <- ggplot() +
  geom_rect(data = era_shade, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = era_shade$fill, alpha = 0.8, inherit.aes = FALSE) +
  geom_line(data = fig2_data, aes(x = year, y = coverage, color = vaccine_code), linewidth = 1.15) +
  geom_point(data = fig2_data, aes(x = year, y = coverage, color = vaccine_code), size = 2) +
  geom_hline(yintercept = 95, linetype = "dashed", linewidth = 0.6, color = "#E74C3C") +
  scale_color_manual(values = c("DTP3" = cols$red, "POL3" = cols$blue, "MCV1" = cols$gold, "MCV2" = "#A65628",
                                "RCV1" = cols$pink, "HIB3" = cols$green, "PCV3" = cols$orange, "ROTAC" = cols$purple)) +
  scale_x_continuous(breaks = seq(1980, 2025, 5)) + scale_y_continuous(limits = c(0, 104), breaks = seq(0, 100, 20)) +
  labs(x = "Year", y = "Coverage (%)", color = "Vaccine") +
  theme_pub(14) + theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) +
  annotate("text", x = 2021.2, y = 97, label = "95% threshold", color = "#E74C3C", size = 4)

ggsave(file.path(OUTPUT_DIR, "figures", "Figure2_Long_Term_Trends.png"), fig2, width = 14, height = 9, dpi = 600)

# ---- Figure 3: COVID-19 shock and recovery ----

fig3_data <- covid_raw %>% filter(vaccine_code %in% table3a_main$Vaccine) %>%
  select(vaccine_code, shock_pp, recovery_pp) %>%
  pivot_longer(cols = c(shock_pp, recovery_pp), names_to = "metric", values_to = "pp_change") %>%
  mutate(vaccine_code = factor(vaccine_code, levels = table3a_order),
         metric = recode(metric, shock_pp = "Shock: 2020 minus 2019", recovery_pp = "Recovery: 2024 minus 2020"))

fig3 <- ggplot(fig3_data, aes(x = vaccine_code, y = pp_change, fill = metric)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  geom_hline(yintercept = 0, linewidth = 0.5, color = cols$dark) +
  scale_fill_manual(values = c("Recovery: 2024 minus 2020" = cols$blue, "Shock: 2020 minus 2019" = cols$orange)) +
  scale_y_continuous(breaks = pretty_breaks(8)) +
  labs(x = "Vaccine", y = "Percentage-point change", fill = NULL) +
  theme_pub(14) + theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTPUT_DIR, "figures", "Figure3_COVID_Shock_Recovery_Main.png"), fig3, width = 13, height = 8, dpi = 600)

# ---- Figure 4: GCC benchmarking ----

fig4_data <- table4_gcc %>% mutate(highlight = ifelse(Country == "Saudi Arabia", "Saudi Arabia", "Other GCC"),
                                   Country = factor(Country, levels = Country[order(`Mean VPDP core`)]))

fig4 <- ggplot(fig4_data, aes(x = `Mean VPDP core`, y = Country, fill = highlight)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.1f%%", `Mean VPDP core`)), hjust = -0.15, size = 5.5, color = "black") +
  scale_fill_manual(values = c("Saudi Arabia" = cols$orange, "Other GCC" = cols$blue)) +
  coord_cartesian(xlim = c(89.5, 100.5), clip = "off") +
  labs(x = "Mean coverage (%)", y = NULL, fill = NULL) +
  theme_pub(14) + theme(legend.position = "none")

ggsave(file.path(OUTPUT_DIR, "figures", "Figure4_GCC_Benchmarking.png"), fig4, width = 14, height = 8.5, dpi = 600)

# ---- Figure 5: HPV scenarios ----

fig5 <- ggplot(hpv_yearly, aes(x = year, y = coverage, color = scenario)) +
  geom_line(linewidth = 1.35) + geom_point(size = 3) +
  geom_hline(yintercept = 90, linetype = "dashed", linewidth = 0.8, color = cols$dark) +
  annotate("text", x = 2028.3, y = 91.5, label = "90% target", color = cols$dark, size = 5) +
  scale_color_manual(values = c("Accelerated scale-up" = cols$green, "Moderate scale-up" = cols$blue,
                                "Required scale-up" = cols$orange, "Status quo slow scale-up" = cols$grey)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(40, 102)) +
  labs(x = "Year", y = "HPV coverage", color = NULL) +
  theme_pub(14) + theme(legend.position = "top")

ggsave(file.path(OUTPUT_DIR, "figures", "Figure5_HPV_ScaleUp_Scenarios.png"), fig5, width = 13, height = 8.5, dpi = 600)

# ---- Supplementary Figure S1: BCG sensitivity ----

bcg_interp_df <- bcg_raw %>% mutate(coverage = ifelse(year == 2019, bcg_interp_2019, coverage))
figS1_data <- bind_rows(bcg_raw %>% mutate(series = "Raw BCG"), bcg_interp_df %>% mutate(series = "Interpolated 2019"))

figS1 <- ggplot(figS1_data, aes(x = year, y = coverage, color = series)) +
  geom_line(linewidth = 1.25) + geom_point(size = 2.7) +
  geom_vline(xintercept = 2019, linetype = "dashed", linewidth = 0.7, color = cols$dark) +
  annotate("text", x = 2019.5, y = 55, label = "2019 outlier", hjust = 0, size = 5) +
  scale_color_manual(values = c("Interpolated 2019" = cols$blue, "Raw BCG" = cols$red)) +
  scale_y_continuous(limits = c(30, 103), breaks = seq(30, 100, 10)) +
  labs(x = "Year", y = "BCG coverage (%)", color = NULL) +
  theme_pub(14) + theme(legend.position = "top")

ggsave(file.path(OUTPUT_DIR, "figures", "Supplementary_Figure_S1_BCG_Sensitivity.png"), figS1, width = 13, height = 8.5, dpi = 600)

# ---- Supplementary Figure S2: Coverage heatmap ----

heatmap_vax <- c("BCG","DTP1","DTP3","HEPB3","HIB3","IPV1","MCV1","MCV2","PCV3","POL3","RCV1","ROTAC")
figS2_data <- ksa_data %>% filter(vaccine_code %in% heatmap_vax) %>% mutate(vaccine_code = factor(vaccine_code, levels = rev(heatmap_vax)))

figS2 <- ggplot(figS2_data, aes(x = year, y = vaccine_code, fill = coverage)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradientn(colors = c("#D73027","#F39C34","#F1C40F","#2ECC71"),
                       values = rescale(c(0, 50, 80, 100)), limits = c(0, 100), na.value = "grey85", name = "Coverage (%)") +
  labs(x = "Year", y = "Vaccine") + theme_pub(13) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(OUTPUT_DIR, "figures", "Supplementary_Figure_S2_Coverage_Heatmap.png"), figS2, width = 15, height = 8.5, dpi = 600)

# ---- Supplementary Figure S3: VPDP-RI ranking ----

figS3_data <- table3b_vpdp %>% filter(Vaccine != "BCG") %>%
  mutate(`VPDP-RI numeric` = as.numeric(`VPDP-RI`), Vaccine = fct_reorder(Vaccine, `VPDP-RI numeric`))

figS3 <- ggplot(figS3_data, aes(x = `VPDP-RI numeric`, y = Vaccine)) +
  geom_col(fill = cols$blue, width = 0.65) +
  geom_text(aes(label = sprintf("%.3f", `VPDP-RI numeric`)), hjust = -0.15, size = 4.8) +
  coord_cartesian(xlim = c(0.975, 1.002), clip = "off") +
  labs(x = "VPDP-RI", y = "Vaccine") + theme_pub(13)

ggsave(file.path(OUTPUT_DIR, "figures", "Supplementary_Figure_S3_VPDP_RI_Ranking.png"), figS3, width = 12, height = 8.5, dpi = 600)

# ==============================================================================
# 10. WORD TABLE EXPORTS
# ==============================================================================

main_doc <- read_docx()
main_doc <- body_add_par(main_doc, "Manuscript-ready main tables", style = "heading 1")
main_doc <- add_to_doc(main_doc, "Table 1. Vaccine-disability relevance classification", table1_classification, fs = 7)
main_doc <- add_to_doc(main_doc, "Table 2. Descriptive statistics and long-term trend analysis", table2_trends, fs = 7)
main_doc <- add_to_doc(main_doc, "Table 3A. COVID-19 coverage shock and recovery", table3a_main,
                       note = "BCG excluded from main COVID-19 shock table due to 2019 outlier; see Supplementary Tables S11A-S11D.", fs = 7)
main_doc <- add_to_doc(main_doc, "Table 3B. VPDP-RI Resilience Index ranking", table3b_vpdp,
                       note = "BCG shown for transparency; interpret cautiously due to 2019 outlier.", fs = 7)
main_doc <- add_to_doc(main_doc, "Table 4. GCC benchmarking, 2023", table4_gcc, fs = 7)
main_doc <- add_to_doc(main_doc, "Table 5A. Latest coverage and 2030 readiness", table5a_readiness, fs = 7)
main_doc <- add_to_doc(main_doc, "Table 5B. HPV scale-up scenarios, 2023-2030", table5b_hpv, fs = 7)
print(main_doc, target = file.path(OUTPUT_DIR, "docx", "JDR_Manuscript_Ready_Main_Tables.docx"))

supp_doc <- read_docx()
supp_doc <- body_add_par(supp_doc, "Supplementary tables", style = "heading 1")
supp_doc <- add_to_doc(supp_doc, "Supplementary Table S1. Dataset screening", dataset_screening, fs = 8)
supp_doc <- add_to_doc(supp_doc, "Supplementary Table S2. Vaccine indicator availability", availability_summary, fs = 7)
supp_doc <- add_to_doc(supp_doc, "Supplementary Table S3. Era-wise mean coverage", table_s3, fs = 7)
supp_doc <- add_to_doc(supp_doc, "Supplementary Table S4. Era observation counts", table_s4, fs = 7)
supp_doc <- add_to_doc(supp_doc, "Supplementary Table S5. COVID-19 full table with BCG and ITS", table_s5, fs = 7)
supp_doc <- add_to_doc(supp_doc, "Supplementary Table S6. GCC vaccine-specific comparisons", table_s6, fs = 7)
supp_doc <- add_to_doc(supp_doc, "Supplementary Table S7. GCC pooled comparison", table_s7, fs = 8)
supp_doc <- add_to_doc(supp_doc, "Supplementary Table S8. Recent trend-based 2030 projections", recent_proj, fs = 7)
supp_doc <- add_to_doc(supp_doc, "Supplementary Table S9. HPV yearly scenarios", hpv_yearly, fs = 7)
supp_doc <- add_to_doc(supp_doc, "Supplementary Table S10. Readiness status summary", table_s10, fs = 8)
supp_doc <- add_to_doc(supp_doc, "Supplementary Table S11A. BCG outlier diagnostic", diag_table, fs = 8)
supp_doc <- add_to_doc(supp_doc, "Supplementary Table S11B. BCG trend sensitivity", bcg_trend_sens, fs = 7)
supp_doc <- add_to_doc(supp_doc, "Supplementary Table S11C. BCG era sensitivity", bcg_era_sens, fs = 7)
supp_doc <- add_to_doc(supp_doc, "Supplementary Table S11D. BCG COVID and VPDP-RI sensitivity", bcg_covid_sens, fs = 7)
print(supp_doc, target = file.path(OUTPUT_DIR, "docx", "JDR_Manuscript_Ready_Supplementary_Tables.docx"))

# ==============================================================================
# 11. SESSION INFO
# ==============================================================================

capture.output(sessionInfo(), file = file.path(OUTPUT_DIR, "session_info.txt"))

message("\nAnalysis completed successfully.")
message("Outputs saved in: ", normalizePath(OUTPUT_DIR, mustWork = FALSE))