rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
})

# =========================================================
# Main workflow
# 1. Cox analysis
# 2. Core mediator selection
# 3. SEM mediation analysis
# =========================================================

# ====================== 1. Global configuration ======================
PROJECT_ROOT <- "."
DATA_DIR <- file.path(PROJECT_ROOT, "data")
PA_DIR <- file.path(DATA_DIR, "PA.ZJFX")
OUT_DIR <- file.path(PA_DIR, "ZJ")
LABEL_DIR <- file.path(DATA_DIR, "Moxingbiaoqian")
FOLLOW_DIR <- file.path(DATA_DIR, "JIXian")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

PATHS <- list(
  follow_time = file.path(FOLLOW_DIR, "follow_time_92407.csv"),
  
  master_nmr = file.path(PA_DIR, "master_bb+bp+nmr.csv"),
  master_pro = file.path(PA_DIR, "master_pro.csv"),
  
  pa_nmr = file.path(PA_DIR, "pa_exposure_mediator_results(bb+bp+nmr).csv"),
  pa_pro = file.path(PA_DIR, "pa_exposure_mediator_results(pro).csv"),
  
  static_features = file.path(LABEL_DIR, "all_static_features.csv"),
  outcome_wear = file.path(LABEL_DIR, "all_icd_f_group_1000_10y_wear.csv"),
  
  cox_nmr = file.path(OUT_DIR, "Cox_Results_bb.NMR_using_follow_time.csv"),
  cox_nmr_sig = file.path(OUT_DIR, "Sig_Biomarkers_bb.NMR_using_follow_time.csv"),
  cox_pro = file.path(OUT_DIR, "Cox_Results_pro_using_follow_time.csv"),
  cox_pro_sig = file.path(OUT_DIR, "Sig_Biomarkers_pro_using_follow_time.csv"),
  
  matched_nmr = file.path(OUT_DIR, "Matched_Disease_PA_Factor_All.csv"),
  core_nmr = file.path(OUT_DIR, "Core_Mediators_Final_Relaxed.csv"),
  matched_pro = file.path(OUT_DIR, "Matched_Disease_PA_Protein_All.csv"),
  core_pro = file.path(OUT_DIR, "Core_Mediators_Protein_Relaxed.csv"),
  
  sem_log_nmr = file.path(OUT_DIR, "SEM_Metab_Inflamm_Log_adjCov.txt"),
  sem_res_nmr = file.path(OUT_DIR, "SEM_Inf_Met_results_adjCov.csv"),
  sem_sum_nmr = file.path(OUT_DIR, "SEM_Inf_Met_summary_adjCov.csv"),
  
  sem_log_pro = file.path(OUT_DIR, "SEM_Protein_Log_adjCov.txt"),
  sem_res_pro = file.path(OUT_DIR, "SEM_Pro_results_adjCov.csv"),
  sem_sum_pro = file.path(OUT_DIR, "SEM_Pro_summary_adjCov.csv")
)

ID_COL_NAME <- "ParticipantID"
MIN_SAMPLE_SIZE <- 20

DISEASE_KEYWORDS <- c("J44", "J45", "J43", "J84", "J47", "J96", "J18", "J22", "J69")

PA_KEYWORDS <- list(
  Light = "LIGHTOVERALLAVERAGEINSTANCE0",
  ModerateVigorous = "MODERATEVIGOROUSOVERALLAVERAGEINSTANCE0",
  Sedentary = "SEDENTARYOVERALLAVERAGEINSTANCE0",
  Sleep = "SLEEPOVERALLAVERAGEINSTANCE0"
)

# ====================== 2. Shared helper functions ======================
read_csv_auto <- function(file, ...) {
  tryCatch(
    read.csv(file, fileEncoding = "UTF-8", ...),
    error = function(e) read.csv(file, fileEncoding = "GBK", ...)
  )
}

clean_id <- function(data) {
  data %>%
    rename(Participant.ID = any_of(c("Participant.ID", "participant_id", "ID", "id", "eid"))) %>%
    mutate(Participant.ID = as.character(Participant.ID)) %>%
    drop_na(Participant.ID) %>%
    distinct(Participant.ID, .keep_all = TRUE)
}

fill_na_median <- function(x) {
  if (is.numeric(x)) replace(x, is.na(x), median(x, na.rm = TRUE)) else x
}

get_mode <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

check_columns <- function(df, required_cols, step_name) {
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    warning(paste0("[", step_name, "] Missing columns: ", paste(missing_cols, collapse = ", ")))
  } else {
    message(paste0("[", step_name, "] Column check passed"))
  }
  df
}

map_pa_short <- function(x) {
  case_when(
    str_detect(x, "Light - Overall average") ~ "LPA",
    str_detect(x, "Moderate-Vigorous - Overall average") ~ "MVPA",
    str_detect(x, "Sedentary - Overall average") ~ "SB",
    str_detect(x, "Sleep - Overall average") ~ "Sleep",
    TRUE ~ NA_character_
  )
}

normalize_nmr_marker <- function(x) {
  x %>%
    gsub(" \\| Instance 0", "", .) %>%
    tolower() %>%
    gsub("[ -/]", "_", .) %>%
    gsub("_+", "_", .)
}

normalize_protein_marker <- function(x) {
  x %>%
    str_trim() %>%
    toupper() %>%
    gsub("[^A-Z0-9]", "_", .) %>%
    gsub("_+", "_", .)
}

ultimate_clean <- function(x) {
  x %>%
    str_replace_all("[^a-zA-Z0-9_]", "") %>%
    str_replace_all("_+", "_") %>%
    str_trim() %>%
    toupper() %>%
    ifelse(. %in% c("PARTICIPANTID", "PARTICIPANT_ID", "EID", "ID", "PID"), "PARTICIPANTID", .)
}

rename_id_to_standard <- function(df) {
  names(df) <- ultimate_clean(names(df))
  if ("PARTICIPANTID" %in% names(df)) names(df)[names(df) == "PARTICIPANTID"] <- ID_COL_NAME
  if ("EID" %in% names(df)) names(df)[names(df) == "EID"] <- ID_COL_NAME
  if ("ID" %in% names(df)) names(df)[names(df) == "ID"] <- ID_COL_NAME
  if (ID_COL_NAME %in% names(df)) df[[ID_COL_NAME]] <- as.character(df[[ID_COL_NAME]])
  df
}

fuzzy_match_col <- function(target_name, col_names) {
  if (is.na(target_name) || target_name == "") return(NA_character_)
  
  clean_target <- ultimate_clean(target_name)
  clean_cols <- tibble(original = col_names, clean = ultimate_clean(col_names)) %>%
    filter(clean != "")
  
  exact_match <- clean_cols %>% filter(clean == clean_target) %>% pull(original)
  if (length(exact_match) > 0) return(exact_match[1])
  
  fuzzy_match <- clean_cols %>% filter(str_detect(clean, clean_target)) %>% pull(original)
  if (length(fuzzy_match) > 0) return(fuzzy_match[1])
  
  prefix_match <- clean_cols %>%
    filter(str_detect(clean, paste0("^", str_sub(clean_target, 1, 5)))) %>%
    pull(original)
  if (length(prefix_match) > 0) return(prefix_match[1])
  
  NA_character_
}

first_non_na <- function(...) {
  vals <- list(...)
  for (v in vals) {
    if (!is.null(v) && length(v) > 0 && !is.na(v[1]) && v[1] != "") return(v[1])
  }
  NA_character_
}

check_pa_col <- function(pa_col_raw, data) {
  if (is.null(pa_col_raw) || pa_col_raw == "") return(NULL)
  clean_pa <- ultimate_clean(pa_col_raw)
  matched_col <- names(data)[ultimate_clean(names(data)) == clean_pa]
  if (length(matched_col) > 0) return(matched_col[1])
  NULL
}

# ====================== 3. Cox analysis ======================
get_available_J_codes <- function(master_df, follow_time_df) {
  j_in_master <- colnames(master_df)[grepl("^J\\d+$", colnames(master_df))]
  j_in_master <- setdiff(j_in_master, "Participant.ID")
  
  has_time <- j_in_master[paste0("event_time_", j_in_master) %in% colnames(follow_time_df)]
  has_stat <- j_in_master[paste0("event_status_", j_in_master) %in% colnames(follow_time_df)]
  
  intersect(has_time, has_stat)
}

build_followup_for_disease <- function(disease_code, master_df, follow_time_df, max_years = 10) {
  time_col <- paste0("event_time_", disease_code)
  status_col <- paste0("event_status_", disease_code)
  
  if (!(time_col %in% names(follow_time_df)) || !(status_col %in% names(follow_time_df))) {
    return(NULL)
  }
  
  out <- follow_time_df %>%
    transmute(
      Participant.ID = as.character(Participant.ID),
      time_days = suppressWarnings(as.numeric(.data[[time_col]])),
      outcome = suppressWarnings(as.integer(.data[[status_col]]))
    ) %>%
    filter(!is.na(time_days), time_days > 0) %>%
    filter(outcome %in% c(0L, 1L)) %>%
    mutate(time = time_days / 365.25) %>%
    select(Participant.ID, time, outcome)
  
  if (disease_code %in% names(follow_time_df)) {
    lab_df <- tibble(
      Participant.ID = as.character(follow_time_df$Participant.ID),
      label = suppressWarnings(as.numeric(follow_time_df[[disease_code]]))
    )
    
    out <- out %>%
      left_join(lab_df, by = "Participant.ID") %>%
      filter(is.na(label) | label != 2) %>%
      select(-label)
  } else if (disease_code %in% names(master_df)) {
    lab_df <- master_df %>%
      transmute(
        Participant.ID = as.character(Participant.ID),
        label = suppressWarnings(as.numeric(.data[[disease_code]]))
      )
    
    out <- out %>%
      left_join(lab_df, by = "Participant.ID") %>%
      filter(is.na(label) | label != 2) %>%
      select(-label)
  }
  
  out %>%
    mutate(time = ifelse(time < 0 | time > max_years, NA_real_, time)) %>%
    drop_na(time, outcome)
}

run_cox_analysis <- function(
    master_file,
    follow_time_df,
    result_file,
    sig_file,
    marker_limit = NULL,
    min_cases = 10,
    min_total = 50,
    dataset_label = "NMR"
) {
  message("========== Running Cox analysis: ", dataset_label, " ==========")
  
  df <- read.csv(
    master_file,
    stringsAsFactors = FALSE,
    na.strings = c("NA", ""),
    fileEncoding = "UTF-8"
  ) %>%
    clean_id()
  
  respiratory_cols <- get_available_J_codes(df, follow_time_df)
  if (length(respiratory_cols) == 0) {
    stop("No eligible J outcomes with both event_time and event_status in follow_time.")
  }
  
  disease_cols_all <- colnames(df)[grepl("^[A-Za-z]+[0-9]+$", colnames(df))]
  basic_cols <- c("age", "sex", "Sex", "Body.mass.index..BMI.._x", "Time.difference_x")
  activity_cols <- c(
    "Light...Overall.average...Instance.0",
    "Moderate.Vigorous...Overall.average...Instance.0",
    "Sedentary...Overall.average...Instance.0",
    "Sleep...Overall.average...Instance.0"
  )
  
  cols_exclude <- intersect(c(disease_cols_all, basic_cols, activity_cols), colnames(df))
  
  all_candidate_markers <- df %>%
    select(-all_of(cols_exclude)) %>%
    select(where(is.numeric)) %>%
    colnames()
  
  variant_markers <- all_candidate_markers[sapply(all_candidate_markers, function(x) {
    sd(df[[x]], na.rm = TRUE) > 0.001
  })]
  
  bio_markers <- if (is.null(marker_limit)) variant_markers else head(variant_markers, marker_limit)
  
  biomarker_stats <- df %>%
    select(all_of(bio_markers)) %>%
    summarise(across(
      everything(),
      list(
        mean = ~mean(., na.rm = TRUE),
        sd = ~sd(., na.rm = TRUE),
        median = ~median(., na.rm = TRUE)
      )
    ))
  
  df_biomarker_std <- df %>%
    clean_id() %>%
    select(Participant.ID, all_of(bio_markers)) %>%
    mutate(across(all_of(bio_markers), ~{
      col_name <- cur_column()
      mean_val <- pull(biomarker_stats, paste0(col_name, "_mean"))
      sd_val <- pull(biomarker_stats, paste0(col_name, "_sd"))
      median_val <- pull(biomarker_stats, paste0(col_name, "_median"))
      
      val_num <- as.numeric(.)
      val_num[is.na(val_num)] <- median_val
      val_num[abs(val_num - mean_val) > 5 * sd_val] <- median_val
      
      if (!is.na(sd_val) && is.finite(sd_val) && sd_val > 0.001) {
        (val_num - mean_val) / sd_val
      } else {
        val_num
      }
    }))
  
  basic_info <- df %>%
    transmute(
      Participant.ID = as.character(Participant.ID),
      age = suppressWarnings(as.numeric(age)),
      Sex = as.character(Sex),
      bmi = suppressWarnings(as.numeric(Body.mass.index..BMI.._x))
    ) %>%
    clean_id() %>%
    mutate(
      sex = ifelse(is.na(Sex) | Sex == "", get_mode(Sex), Sex),
      age_std = ifelse(is.na(age), NA_real_, (age - mean(age, na.rm = TRUE)) / sd(age, na.rm = TRUE)),
      bmi_std = ifelse(is.na(bmi), NA_real_, (bmi - mean(bmi, na.rm = TRUE)) / sd(bmi, na.rm = TRUE))
    ) %>%
    select(Participant.ID, age_std, sex, bmi_std)
  
  analysis_datasets <- list()
  
  for (disease_col in respiratory_cols) {
    fu <- build_followup_for_disease(disease_col, df, follow_time_df, max_years = 10)
    if (is.null(fu)) next
    
    disease_analysis <- basic_info %>%
      left_join(fu, by = "Participant.ID") %>%
      left_join(df_biomarker_std, by = "Participant.ID") %>%
      mutate(sex = as.factor(sex)) %>%
      drop_na(time, outcome) %>%
      distinct(Participant.ID, .keep_all = TRUE)
    
    case_count <- sum(disease_analysis$outcome == 1, na.rm = TRUE)
    total_count <- nrow(disease_analysis)
    
    if (case_count >= min_cases && total_count >= min_total) {
      analysis_datasets[[disease_col]] <- disease_analysis
    }
  }
  
  if (length(analysis_datasets) == 0) {
    stop(paste0("No ", dataset_label, " outcomes met the sample size requirement."))
  }
  
  valid_markers <- map_lgl(bio_markers, function(marker) {
    all_sd <- map_dbl(analysis_datasets, function(dat) sd(dat[[marker]], na.rm = TRUE))
    any(all_sd > 0)
  })
  
  sig_biomarkers <- bio_markers[valid_markers]
  
  all_results <- data.frame(
    Disease_Name = character(),
    Variable_Name = character(),
    Variable_Type = character(),
    Coefficient = numeric(),
    HR = numeric(),
    HR_95CI_Lower = numeric(),
    HR_95CI_Upper = numeric(),
    HR_95CI = character(),
    Multivariate_P = numeric(),
    Multivariate_FDR_Q = numeric(),
    Effect_Direction = character(),
    stringsAsFactors = FALSE
  )
  
  confounders <- c("age_std", "sex", "bmi_std")
  all_variables <- c(confounders, sig_biomarkers)
  
  for (disease_name in names(analysis_datasets)) {
    dat <- analysis_datasets[[disease_name]]
    
    model_data <- dat %>%
      select(all_of(all_variables), time, outcome) %>%
      mutate(
        sex = as.factor(sex),
        across(where(is.numeric), fill_na_median)
      ) %>%
      drop_na(time, outcome)
    
    for (var in all_variables) {
      tryCatch({
        if (var %in% confounders) {
          formula <- as.formula(paste("Surv(time, outcome) ~", var))
        } else {
          formula <- as.formula(paste(
            "Surv(time, outcome) ~", var, "+", paste(confounders, collapse = " + ")
          ))
        }
        
        cox_model <- coxph(formula, data = model_data)
        cox_summary <- summary(cox_model)
        
        if (var %in% rownames(cox_summary$coefficients)) {
          coef_row <- cox_summary$coefficients[var, ]
          hr <- exp(coef_row["coef"])
          
          ci <- confint(cox_model, var)
          hr_lower <- exp(ci[1])
          hr_upper <- exp(ci[2])
          
          effect_dir <- case_when(
            is.na(hr) ~ "Model_Failed",
            hr > 1 ~ "Risk_Factor",
            hr < 1 ~ "Protective_Factor",
            TRUE ~ "No_Effect"
          )
          
          all_results <- bind_rows(all_results, data.frame(
            Disease_Name = disease_name,
            Variable_Name = var,
            Variable_Type = ifelse(var %in% confounders, "Confounder", "Biomarker"),
            Coefficient = round(coef_row["coef"], 4),
            HR = round(hr, 4),
            HR_95CI_Lower = round(hr_lower, 3),
            HR_95CI_Upper = round(hr_upper, 3),
            HR_95CI = paste0("(", round(hr_lower, 3), "-", round(hr_upper, 3), ")"),
            Multivariate_P = coef_row["Pr(>|z|)"],
            Multivariate_FDR_Q = NA,
            Effect_Direction = effect_dir,
            stringsAsFactors = FALSE
          ))
        }
      }, error = function(e) {
        all_results <- bind_rows(all_results, data.frame(
          Disease_Name = disease_name,
          Variable_Name = var,
          Variable_Type = ifelse(var %in% confounders, "Confounder", "Biomarker"),
          Coefficient = NA,
          HR = NA,
          HR_95CI_Lower = NA,
          HR_95CI_Upper = NA,
          HR_95CI = "Model_Failed",
          Multivariate_P = NA,
          Multivariate_FDR_Q = NA,
          Effect_Direction = "Model_Failed",
          stringsAsFactors = FALSE
        ))
      })
    }
  }
  
  all_results <- all_results %>%
    group_by(Disease_Name) %>%
    mutate(
      Multivariate_P_Clean = ifelse(
        is.na(Multivariate_P) | Multivariate_P < 0 | Multivariate_P > 1,
        NA,
        Multivariate_P
      ),
      Multivariate_FDR_Q = p.adjust(Multivariate_P_Clean, method = "BH")
    ) %>%
    ungroup() %>%
    select(-Multivariate_P_Clean)
  
  write.csv(all_results, result_file, row.names = FALSE, fileEncoding = "UTF-8")
  
  significant_results <- all_results %>%
    filter(!is.na(Multivariate_FDR_Q) & Multivariate_FDR_Q < 0.05)
  
  if (nrow(significant_results) > 0) {
    write.csv(significant_results, sig_file, row.names = FALSE, fileEncoding = "UTF-8")
  } else {
    write.csv(
      data.frame(Note = "No significant biomarkers found"),
      sig_file,
      row.names = FALSE,
      fileEncoding = "UTF-8"
    )
  }
  
  invisible(list(all_results = all_results, significant_results = significant_results))
}

# ====================== 4. Core mediator selection ======================
select_core_mediators_nmr <- function(pa_file, cox_file, matched_file, core_file) {
  required_cols <- c(
    "mediator_core", "pa_short", "Mediator_Group",
    "PA_Beta", "HR", "Disease_Name", "PA_FDR", "Disease_FDR"
  )
  
  pa_data <- read_csv_auto(pa_file)
  
  significant_pa_markers <- pa_data %>%
    filter(!is.na(FDR) & FDR < 0.15) %>%
    arrange(FDR) %>%
    select(Exposure_Type, mediator, beta, p, FDR, Mediator_Group) %>%
    mutate(
      mediator_core = normalize_nmr_marker(mediator),
      pa_short = map_pa_short(Exposure_Type)
    ) %>%
    drop_na(pa_short) %>%
    rename(PA_Beta = beta, PA_FDR = FDR)
  
  cox_data <- read.csv(cox_file) %>%
    mutate(marker_core = normalize_nmr_marker(gsub("\\.\\.\\.Instance\\.0", "", Variable_Name))) %>%
    filter(Variable_Type == "Biomarker") %>%
    drop_na(Disease_Name, HR, Multivariate_FDR_Q) %>%
    rename(Disease_FDR = Multivariate_FDR_Q)
  
  integrated_data <- significant_pa_markers %>%
    select(pa_short, mediator_core, Mediator_Group, PA_Beta, PA_FDR) %>%
    inner_join(
      cox_data %>% select(Disease_Name, marker_core, HR, Disease_FDR),
      by = c("mediator_core" = "marker_core"),
      relationship = "many-to-many"
    ) %>%
    distinct(Disease_Name, pa_short, mediator_core, .keep_all = TRUE) %>%
    arrange(Disease_Name, pa_short, PA_FDR, Disease_FDR)
  
  significant_mediators <- integrated_data %>%
    filter(PA_FDR < 0.15 & Disease_FDR < 0.15) %>%
    check_columns(required_cols, "NMR mediator filtering")
  
  write.csv(significant_mediators, matched_file, row.names = FALSE, fileEncoding = "UTF-8")
  
  inflammation_keywords <- c("baso", "eosino", "haem", "retic", "lymph", "neutro", "leuko", "inflamm")
  metabolite_keywords <- c("gluc", "chol", "trigly", "lipid", "hba1c", "album", "creat", "urea", "bilir", "metab", "glyc", "insul")
  
  mediators_bio <- significant_mediators %>%
    filter(
      (Mediator_Group == "Inflammation" &
         str_detect(mediator_core, paste(inflammation_keywords, collapse = "|"))) |
        (Mediator_Group == "Metabolite" &
           str_detect(mediator_core, paste(metabolite_keywords, collapse = "|")))
    ) %>%
    check_columns(required_cols, "NMR biological pathway filtering")
  
  factor_pa_count <- mediators_bio %>%
    group_by(mediator_core) %>%
    summarise(pa_occur = n_distinct(pa_short), .groups = "drop")
  
  mediators_final <- mediators_bio %>%
    left_join(factor_pa_count, by = "mediator_core") %>%
    filter(pa_occur <= 6) %>%
    select(-pa_occur) %>%
    group_by(Disease_Name, pa_short) %>%
    slice_min(order_by = PA_FDR + Disease_FDR, n = 50, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(Disease_Name, pa_short, PA_FDR, Disease_FDR)
  
  write.csv(mediators_final, core_file, row.names = FALSE, fileEncoding = "UTF-8")
  invisible(mediators_final)
}

select_core_mediators_protein <- function(pa_file, cox_file, matched_file, core_file) {
  required_cols <- c(
    "mediator_core", "pa_short", "Mediator_Group",
    "PA_Beta", "HR", "Disease_Name", "PA_FDR", "Disease_FDR"
  )
  
  pa_data <- read_csv_auto(pa_file)
  
  significant_pa_markers <- pa_data %>%
    select(Exposure_Type, mediator, beta, p, FDR) %>%
    mutate(pa_short = map_pa_short(Exposure_Type)) %>%
    drop_na(pa_short) %>%
    filter(!is.na(FDR) & (
      (pa_short == "Sleep" & FDR < 0.05) |
        (pa_short != "Sleep" & FDR < 0.15)
    )) %>%
    arrange(FDR) %>%
    mutate(
      mediator_short = str_trim(str_split(mediator, ";", simplify = TRUE)[, 1]),
      mediator_clean1 = gsub(" \\| Instance 0", "", mediator_short),
      mediator_core = normalize_protein_marker(mediator_clean1)
    ) %>%
    rename(PA_Beta = beta, PA_FDR = FDR) %>%
    mutate(Mediator_Group = "Protein")
  
  cox_data <- read.csv(cox_file) %>%
    mutate(
      marker_short = str_trim(str_split(Variable_Name, "\\.", simplify = TRUE)[, 1]),
      marker_core = normalize_protein_marker(marker_short)
    ) %>%
    filter(Variable_Type == "Biomarker") %>%
    drop_na(Disease_Name, HR, Multivariate_FDR_Q) %>%
    rename(Disease_FDR = Multivariate_FDR_Q)
  
  integrated_data <- significant_pa_markers %>%
    select(pa_short, mediator_core, Mediator_Group, PA_Beta, PA_FDR) %>%
    inner_join(
      cox_data %>% select(Disease_Name, marker_core, HR, Disease_FDR),
      by = c("mediator_core" = "marker_core"),
      relationship = "many-to-many"
    ) %>%
    distinct(Disease_Name, pa_short, mediator_core, .keep_all = TRUE) %>%
    arrange(Disease_Name, pa_short, PA_FDR, Disease_FDR)
  
  significant_mediators <- integrated_data %>%
    filter(
      (pa_short == "Sleep" & PA_FDR < 0.05 & Disease_FDR < 0.05) |
        (pa_short != "Sleep" & PA_FDR < 0.15 & Disease_FDR < 0.15)
    ) %>%
    check_columns(required_cols, "Protein mediator filtering")
  
  write.csv(significant_mediators, matched_file, row.names = FALSE, fileEncoding = "UTF-8")
  
  protein_keywords <- c(
    "ITGAV", "MYDGF", "CRP", "IL", "TNF", "ADIPO", "LEPTIN",
    "INTEGRIN", "GROWTH", "CYTOKINE", "CHEMOKINE",
    "SLEEP", "MELATONIN", "CORTISOL", "HSP", "NFKB"
  )
  
  mediators_bio <- significant_mediators %>%
    filter(str_detect(toupper(mediator_core), paste(protein_keywords, collapse = "|"))) %>%
    check_columns(required_cols, "Protein biological pathway filtering")
  
  factor_pa_count <- mediators_bio %>%
    group_by(mediator_core) %>%
    summarise(pa_occur = n_distinct(pa_short), .groups = "drop")
  
  mediators_specific <- mediators_bio %>%
    left_join(factor_pa_count, by = "mediator_core") %>%
    filter(pa_occur <= 5) %>%
    select(-pa_occur)
  
  mediators_sleep <- mediators_specific %>%
    filter(pa_short == "Sleep") %>%
    group_by(Disease_Name, pa_short) %>%
    slice_min(order_by = PA_FDR + Disease_FDR, n = 50, with_ties = FALSE) %>%
    ungroup()
  
  mediators_other <- mediators_specific %>%
    filter(pa_short != "Sleep") %>%
    group_by(Disease_Name, pa_short) %>%
    slice_min(order_by = PA_FDR + Disease_FDR, n = 30, with_ties = FALSE) %>%
    ungroup()
  
  mediators_final <- bind_rows(mediators_sleep, mediators_other) %>%
    arrange(Disease_Name, pa_short, PA_FDR, Disease_FDR)
  
  write.csv(mediators_final, core_file, row.names = FALSE, fileEncoding = "UTF-8")
  invisible(mediators_final)
}

# ====================== 5. SEM mediation analysis ======================
build_cov_df <- function(static_df) {
  age_col <- fuzzy_match_col("AGE", names(static_df))
  sex_col <- fuzzy_match_col("SEX", names(static_df))
  bmi_col <- fuzzy_match_col("BODYMASSINDEXBMI", names(static_df))
  if (is.na(bmi_col)) bmi_col <- fuzzy_match_col("BMI", names(static_df))
  
  if (is.na(age_col) || is.na(sex_col) || is.na(bmi_col)) {
    stop(paste0("static missing AGE/SEX/BMI: age=", age_col, " sex=", sex_col, " bmi=", bmi_col))
  }
  
  static_df %>%
    select(all_of(c(ID_COL_NAME, age_col, sex_col, bmi_col))) %>%
    rename(
      AGE_COV = !!sym(age_col),
      SEX_COV = !!sym(sex_col),
      BMI_COV = !!sym(bmi_col)
    ) %>%
    mutate(!!ID_COL_NAME := as.character(.data[[ID_COL_NAME]])) %>%
    distinct(.data[[ID_COL_NAME]], .keep_all = TRUE)
}

add_cov_z <- function(df) {
  if (!all(c("AGE", "BMI", "SEX") %in% names(df))) {
    stop("add_cov_z(): missing AGE/BMI/SEX columns")
  }
  
  df <- df %>%
    mutate(
      AGE = suppressWarnings(as.numeric(AGE)),
      BMI = suppressWarnings(as.numeric(BMI)),
      SEX = as.character(SEX)
    )
  
  age_med <- suppressWarnings(median(df$AGE, na.rm = TRUE))
  bmi_med <- suppressWarnings(median(df$BMI, na.rm = TRUE))
  sex_mode <- get_mode(df$SEX)
  
  df$AGE[is.na(df$AGE)] <- age_med
  df$BMI[is.na(df$BMI)] <- bmi_med
  df$SEX[is.na(df$SEX) | df$SEX == ""] <- sex_mode
  df$SEX <- as.factor(df$SEX)
  
  age_sd <- sd(df$AGE, na.rm = TRUE)
  bmi_sd <- sd(df$BMI, na.rm = TRUE)
  
  df %>%
    mutate(
      AGE_Z = ifelse(is.finite(age_sd) && age_sd > 0, (AGE - mean(AGE, na.rm = TRUE)) / age_sd, AGE),
      BMI_Z = ifelse(is.finite(bmi_sd) && bmi_sd > 0, (BMI - mean(BMI, na.rm = TRUE)) / bmi_sd, BMI)
    )
}

choose_outcome_cols <- function(outcome_pref, disease_codes) {
  candidate_cols <- names(outcome_pref)[
    sapply(names(outcome_pref), function(col) {
      if (col == ID_COL_NAME) return(FALSE)
      vals <- na.omit(outcome_pref[[col]])
      any(str_detect(col, paste(disease_codes, collapse = "|"))) &&
        length(unique(vals)) > 1 &&
        sum(vals == 1, na.rm = TRUE) > 0 &&
        sum(vals == 0, na.rm = TRUE) > 0
    })
  ]
  
  valid_outcomes <- map_dfr(disease_codes, function(code) {
    cols <- candidate_cols[str_detect(candidate_cols, code)]
    if (length(cols) == 0) return(tibble(DISEASE_COL = NA_character_, DISEASE_CODE = code))
    
    x_col <- cols[str_detect(cols, paste0(code, "\\.X$|", code, "\\.x$"))]
    y_col <- cols[str_detect(cols, paste0(code, "\\.Y$|", code, "\\.y$"))]
    raw_col <- cols[!str_detect(cols, "\\.X$|\\.Y$|\\.x$|\\.y$")]
    
    final_col <- case_when(
      length(x_col) > 0 ~ x_col[1],
      length(y_col) > 0 ~ y_col[1],
      length(raw_col) > 0 ~ raw_col[1],
      TRUE ~ NA_character_
    )
    
    tibble(DISEASE_COL = final_col, DISEASE_CODE = code)
  }) %>%
    filter(!is.na(DISEASE_COL))
  
  valid_outcomes
}

read_core_mediators_by_disease <- function(file_path, factor_type, logger = message) {
  logger(paste0("Reading core mediator file: ", file_path))
  
  df <- read_csv_auto(file_path, check.names = FALSE)
  names(df) <- ultimate_clean(names(df))
  
  pa_candidates <- names(df)[names(df) %in% c("PASHORT", "PATYPE", "PATYPERAW")]
  mediator_candidates <- names(df)[names(df) %in% c("MEDIATORCORE", "MEDIATORNAME", "VARIABLE", "FACTOR")]
  disease_candidates <- names(df)[names(df) %in% c("DISEASENAME", "DISEASE")]
  group_candidates <- names(df)[names(df) %in% c("MEDIATORGROUP", "GROUP")]
  
  pa_col <- pa_candidates[1]
  mediator_col <- mediator_candidates[1]
  disease_col <- disease_candidates[1]
  group_col <- group_candidates[1]
  
  if (is.na(pa_col) || is.na(mediator_col) || is.na(disease_col)) {
    stop(paste0("Required PA / Mediator / Disease column not found: ", file_path))
  }
  
  out <- df %>%
    rename(
      PA_TYPE_RAW = !!sym(pa_col),
      MEDIATOR_NAME = !!sym(mediator_col),
      DISEASE_NAME = !!sym(disease_col)
    )
  
  if (!is.na(group_col)) {
    out <- out %>% rename(Mediator_Group = !!sym(group_col))
  } else {
    out <- out %>% mutate(Mediator_Group = factor_type)
  }
  
  out <- out %>%
    mutate(
      PA_CATEGORY = case_when(
        str_detect(toupper(PA_TYPE_RAW), "LIGHT|LPA") ~ "Light",
        str_detect(toupper(PA_TYPE_RAW), "MODERATE|MVPA|VIGOROUS") ~ "ModerateVigorous",
        str_detect(toupper(PA_TYPE_RAW), "SEDENTARY|SB") ~ "Sedentary",
        str_detect(toupper(PA_TYPE_RAW), "SLEEP") ~ "Sleep",
        TRUE ~ PA_TYPE_RAW
      ),
      DISEASE_CODE = str_extract(DISEASE_NAME, "J\\d+")
    ) %>%
    filter(!is.na(DISEASE_CODE) & DISEASE_CODE %in% DISEASE_KEYWORDS) %>%
    distinct(DISEASE_CODE, PA_CATEGORY, MEDIATOR_NAME, Mediator_Group) %>%
    rename(PA_TYPE = PA_CATEGORY) %>%
    mutate(FACTOR_TYPE = factor_type) %>%
    drop_na(MEDIATOR_NAME, DISEASE_CODE)
  
  logger(paste0("Extracted core mediators: ", nrow(out), " (", factor_type, ")"))
  out
}

strip_bt <- function(x) gsub("`", "", x)

get_term <- function(model, term) {
  ct <- coef(summary(model))
  rn <- strip_bt(rownames(ct))
  term2 <- strip_bt(term)
  
  if (!(term2 %in% rn)) return(NULL)
  idx <- match(term2, rn)
  
  pcol <- if ("Pr(>|z|)" %in% colnames(ct)) "Pr(>|z|)" else "Pr(>|t|)"
  list(
    est = ct[idx, "Estimate"],
    se = ct[idx, "Std. Error"],
    p = ct[idx, pcol]
  )
}

run_single_factor_sem <- function(data, pa_type, mediator, disease, group, mediator_group = NA, log_info = NULL) {
  need <- c(ID_COL_NAME, pa_type, mediator, disease, "AGE_Z", "BMI_Z", "SEX")
  if (!all(need %in% names(data))) {
    if (!is.null(log_info)) log_info("No result: required columns missing")
    return(NULL)
  }
  
  data_filtered <- data %>% drop_na(all_of(need))
  if (nrow(data_filtered) < MIN_SAMPLE_SIZE) {
    if (!is.null(log_info)) log_info("No result: insufficient sample size after drop_na")
    return(NULL)
  }
  
  var_ok <- function(col) length(unique(na.omit(data_filtered[[col]]))) > 1
  if (!var_ok(pa_type) || !var_ok(mediator) || !var_ok(disease)) {
    if (!is.null(log_info)) log_info("No result: no variation in PA, mediator, or disease")
    return(NULL)
  }
  
  is_binary <- all(na.omit(data_filtered[[disease]]) %in% c(0, 1))
  
  cov_vec <- c("AGE_Z", "BMI_Z")
  if (length(unique(data_filtered$SEX)) > 1) cov_vec <- c(cov_vec, "SEX")
  cov_terms <- paste(cov_vec, collapse = " + ")
  
  f_a <- as.formula(paste0("`", mediator, "` ~ `", pa_type, "` + ", cov_terms))
  f_b <- as.formula(paste0("`", disease, "` ~ `", pa_type, "` + `", mediator, "` + ", cov_terms))
  f_c <- as.formula(paste0("`", disease, "` ~ `", pa_type, "` + ", cov_terms))
  
  tryCatch({
    m1 <- lm(f_a, data = data_filtered)
    a <- get_term(m1, pa_type)
    if (is.null(a)) {
      if (!is.null(log_info)) log_info("No result: failed to extract PA coefficient from path a")
      return(NULL)
    }
    
    m2 <- if (is_binary) glm(f_b, data = data_filtered, family = binomial("logit")) else lm(f_b, data = data_filtered)
    b <- get_term(m2, mediator)
    cp <- get_term(m2, pa_type)
    if (is.null(b) || is.null(cp)) {
      if (!is.null(log_info)) log_info("No result: failed to extract coefficients from path b or c_prime")
      return(NULL)
    }
    
    m0 <- if (is_binary) glm(f_c, data = data_filtered, family = binomial("logit")) else lm(f_c, data = data_filtered)
    c0 <- get_term(m0, pa_type)
    if (is.null(c0)) {
      if (!is.null(log_info)) log_info("No result: failed to extract PA coefficient from total effect model")
      return(NULL)
    }
    
    indirect_effect <- a$est * b$est
    total_effect <- c0$est
    indirect_ratio <- ifelse(total_effect != 0, indirect_effect / total_effect, NA_real_)
    
    sobel_z <- indirect_effect / sqrt((b$est^2) * (a$se^2) + (a$est^2) * (b$se^2))
    ie_p <- 2 * pnorm(-abs(sobel_z))
    
    signif_star <- function(p) {
      dplyr::case_when(
        is.na(p) ~ "NA",
        p < 0.001 ~ "***",
        p < 0.01 ~ "**",
        p < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    }
    
    tibble(
      Group = group,
      Mediator_Group = mediator_group,
      PA_Type = pa_type,
      Mediator = mediator,
      Disease = disease,
      Outcome_Type = ifelse(is_binary, "Binary", "Continuous"),
      Sample_Size = nrow(data_filtered),
      a_coef = a$est,
      a_p = a$p,
      b_coef = b$est,
      b_p = b$p,
      c_coef = c0$est,
      c_p = c0$p,
      c_prime_coef = cp$est,
      c_prime_p = cp$p,
      indirect_effect = indirect_effect,
      total_effect = total_effect,
      indirect_ratio = indirect_ratio,
      ie_p = ie_p,
      a_signif = signif_star(a$p),
      b_signif = signif_star(b$p),
      c_signif = signif_star(c0$p),
      c_prime_signif = signif_star(cp$p),
      ie_signif = signif_star(ie_p),
      indirect_significant = ifelse(!is.na(ie_p) & ie_p < 0.05, "Significant", "Not Significant"),
      model_fit = "Success"
    )
  }, error = function(e) {
    if (!is.null(log_info)) log_info(paste0("Model error: ", substr(e$message, 1, 140)))
    NULL
  })
}

run_sem_pipeline <- function(
    pipeline_name,
    factor_file,
    static_file,
    master_file,
    outcome_file,
    log_file,
    out_result_file,
    out_summary_file
) {
  if (file.exists(log_file)) file.remove(log_file)
  
  log_info <- function(msg) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    log_msg <- paste0("[", timestamp, "] ", msg)
    cat(log_msg, "\n")
    write(log_msg, log_file, append = TRUE)
  }
  
  log_info(paste0("========== Start ", pipeline_name, " SEM mediation analysis =========="))
  
  factors <- read_core_mediators_by_disease(factor_file, pipeline_name, logger = log_info)
  
  static <- read.csv(static_file, check.names = FALSE) %>% rename_id_to_standard()
  master_raw <- read.csv(master_file, check.names = FALSE) %>% rename_id_to_standard()
  outcome_raw <- read.csv(outcome_file, check.names = FALSE) %>% rename_id_to_standard()
  
  outcome_clean <- outcome_raw %>%
    mutate(across(-all_of(ID_COL_NAME), ~ ifelse(.x == 2, NA, .x)))
  
  common_ids <- intersect(static[[ID_COL_NAME]], intersect(master_raw[[ID_COL_NAME]], outcome_clean[[ID_COL_NAME]]))
  log_info(paste0("Common sample size: ", length(common_ids)))
  
  master_filtered <- master_raw %>% filter(.data[[ID_COL_NAME]] %in% common_ids)
  outcome_filtered <- outcome_clean %>% filter(.data[[ID_COL_NAME]] %in% common_ids)
  static_filtered <- static %>% filter(.data[[ID_COL_NAME]] %in% common_ids)
  
  cov_df <- build_cov_df(static_filtered)
  
  outcome_pref <- outcome_filtered %>%
    rename_with(~ paste0("OUT_", .x), -all_of(ID_COL_NAME))
  
  valid_outcomes <- choose_outcome_cols(outcome_pref, DISEASE_KEYWORDS)
  
  if (nrow(valid_outcomes) == 0) {
    log_info("No valid disease columns met the 0/1 variation requirement. Stop.")
    return(invisible(NULL))
  }
  
  analysis_data <- master_filtered %>%
    inner_join(outcome_pref, by = ID_COL_NAME) %>%
    left_join(cov_df, by = ID_COL_NAME) %>%
    mutate(!!ID_COL_NAME := as.character(.data[[ID_COL_NAME]]))
  
  pick_non_cov <- function(target, nms) {
    nms2 <- nms[!str_detect(nms, "_COV$") & !str_detect(nms, "^OUT_")]
    fuzzy_match_col(target, nms2)
  }
  
  age0 <- pick_non_cov("AGE", names(analysis_data))
  sex0 <- pick_non_cov("SEX", names(analysis_data))
  bmi0 <- pick_non_cov("BMI", names(analysis_data))
  if (is.na(bmi0)) bmi0 <- pick_non_cov("BODYMASSINDEXBMI", names(analysis_data))
  
  analysis_data$AGE <- if (!is.na(age0)) analysis_data[[age0]] else NA_real_
  analysis_data$SEX <- if (!is.na(sex0)) analysis_data[[sex0]] else NA_character_
  analysis_data$BMI <- if (!is.na(bmi0)) analysis_data[[bmi0]] else NA_real_
  
  analysis_data <- analysis_data %>%
    mutate(
      AGE = coalesce(suppressWarnings(as.numeric(AGE)), suppressWarnings(as.numeric(AGE_COV))),
      BMI = coalesce(suppressWarnings(as.numeric(BMI)), suppressWarnings(as.numeric(BMI_COV))),
      SEX = coalesce(as.character(SEX), as.character(SEX_COV))
    ) %>%
    select(-any_of(c("AGE_COV", "BMI_COV", "SEX_COV"))) %>%
    add_cov_z()
  
  exclude_cols <- c(
    ID_COL_NAME,
    valid_outcomes$DISEASE_COL,
    unlist(PA_KEYWORDS),
    "AGE", "BMI", "SEX", "AGE_Z", "BMI_Z"
  )
  exclude_cols <- intersect(exclude_cols, names(analysis_data))
  
  numeric_std_cols <- names(analysis_data)[sapply(analysis_data, is.numeric)]
  numeric_std_cols <- setdiff(numeric_std_cols, exclude_cols)
  
  if (length(numeric_std_cols) > 0) {
    analysis_data <- analysis_data %>%
      mutate(
        across(
          all_of(numeric_std_cols),
          ~ ifelse(sd(.x, na.rm = TRUE) == 0, NA, as.numeric(scale(.x))),
          .names = "{.col}_STD"
        )
      )
  }
  
  matched <- factors %>%
    inner_join(valid_outcomes, by = c("DISEASE_CODE")) %>%
    rowwise() %>%
    mutate(
      clean_mediator = str_replace_all(MEDIATOR_NAME, "[^a-zA-Z0-9]", ""),
      RAW_COL_NAME = first_non_na(
        fuzzy_match_col(MEDIATOR_NAME, names(master_filtered)),
        fuzzy_match_col(clean_mediator, names(master_filtered)),
        fuzzy_match_col(str_sub(clean_mediator, 1, 6), names(master_filtered)),
        fuzzy_match_col(paste0(MEDIATOR_NAME, "_STD"), names(analysis_data)),
        fuzzy_match_col(paste0(clean_mediator, "_STD"), names(analysis_data))
      ),
      col_exists = !is.na(RAW_COL_NAME),
      col_has_data = if (col_exists) sum(!is.na(analysis_data[[RAW_COL_NAME]])) > MIN_SAMPLE_SIZE else FALSE
    ) %>%
    ungroup() %>%
    filter(col_exists & col_has_data)
  
  if (nrow(matched) == 0) {
    log_info("No factor passed column matching and sample size requirements. Stop.")
    return(invisible(NULL))
  }
  
  all_results <- list()
  total_steps <- nrow(matched)
  current_step <- 0
  
  walk(1:nrow(matched), function(i) {
    row <- matched[i, ]
    current_step <<- current_step + 1
    
    pa_col_raw <- PA_KEYWORDS[[row$PA_TYPE]]
    if (is.null(pa_col_raw) || is.na(pa_col_raw)) pa_col_raw <- row$PA_TYPE
    pa_col <- check_pa_col(pa_col_raw, analysis_data)
    
    if (is.null(pa_col) || !pa_col %in% names(analysis_data)) {
      log_info(paste0("Step ", current_step, "/", total_steps, ": invalid PA column, skipped"))
      return()
    }
    
    res <- run_single_factor_sem(
      data = analysis_data,
      pa_type = pa_col,
      mediator = row$RAW_COL_NAME,
      disease = row$DISEASE_COL,
      group = pipeline_name,
      mediator_group = row$Mediator_Group,
      log_info = log_info
    )
    
    if (!is.null(res)) {
      all_results[[length(all_results) + 1]] <<- res
    }
  })
  
  if (length(all_results) == 0) {
    log_info("No valid output generated.")
    return(invisible(NULL))
  }
  
  final_res <- bind_rows(all_results)
  
  summary_tbl <- final_res %>%
    mutate(DISEASE_CODE = str_extract(Disease, "J\\d+")) %>%
    group_by(DISEASE_CODE, Mediator_Group) %>%
    summarise(
      successful_models = n(),
      unique_mediators = n_distinct(Mediator),
      significant_indirect_effects = sum(ie_p < 0.05, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(significant_rate = paste0(round(significant_indirect_effects / successful_models * 100, 1), "%"))
  
  write.csv(final_res, out_result_file, row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(summary_tbl, out_summary_file, row.names = FALSE, fileEncoding = "UTF-8")
  
  log_info(paste0("Results exported: ", out_result_file))
  log_info(paste0("Summary exported: ", out_summary_file))
  log_info(paste0("========== ", pipeline_name, " SEM analysis completed =========="))
  
  invisible(list(results = final_res, summary = summary_tbl))
}

# ====================== 6. Run full pipeline ======================
follow_time <- read.csv(
  PATHS$follow_time,
  stringsAsFactors = FALSE,
  na.strings = c("NA", "")
) %>%
  clean_id()

run_cox_analysis(
  master_file = PATHS$master_nmr,
  follow_time_df = follow_time,
  result_file = PATHS$cox_nmr,
  sig_file = PATHS$cox_nmr_sig,
  marker_limit = 159,
  min_cases = 10,
  min_total = 50,
  dataset_label = "NMR"
)

run_cox_analysis(
  master_file = PATHS$master_pro,
  follow_time_df = follow_time,
  result_file = PATHS$cox_pro,
  sig_file = PATHS$cox_pro_sig,
  marker_limit = NULL,
  min_cases = 5,
  min_total = 50,
  dataset_label = "Proteomics"
)

select_core_mediators_nmr(
  pa_file = PATHS$pa_nmr,
  cox_file = PATHS$cox_nmr,
  matched_file = PATHS$matched_nmr,
  core_file = PATHS$core_nmr
)

select_core_mediators_protein(
  pa_file = PATHS$pa_pro,
  cox_file = PATHS$cox_pro,
  matched_file = PATHS$matched_pro,
  core_file = PATHS$core_pro
)

run_sem_pipeline(
  pipeline_name = "Metab_Inflamm",
  factor_file = PATHS$core_nmr,
  static_file = PATHS$static_features,
  master_file = PATHS$master_nmr,
  outcome_file = PATHS$outcome_wear,
  log_file = PATHS$sem_log_nmr,
  out_result_file = PATHS$sem_res_nmr,
  out_summary_file = PATHS$sem_sum_nmr
)

run_sem_pipeline(
  pipeline_name = "Protein",
  factor_file = PATHS$core_pro,
  static_file = PATHS$static_features,
  master_file = PATHS$master_pro,
  outcome_file = PATHS$outcome_wear,
  log_file = PATHS$sem_log_pro,
  out_result_file = PATHS$sem_res_pro,
  out_summary_file = PATHS$sem_sum_pro
)