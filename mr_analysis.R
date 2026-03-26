suppressPackageStartupMessages({
  library(dplyr)
  library(TwoSampleMR)
  library(purrr)
  library(data.table)
})

gc()

# =========================================================
# Batch Mendelian randomization analysis
# =========================================================

# ====================== 1. Global configuration ======================
PROJECT_ROOT <- "."
DATA_DIR <- file.path(PROJECT_ROOT, "data")

GWAS_DIR <- file.path(DATA_DIR, "MR.gwas")
LD_DIR <- file.path(DATA_DIR, "MR.LD")
OUTCOME_DIR <- file.path(DATA_DIR, "MR.finna")
OUTPUT_DIR <- file.path(DATA_DIR, "MR.final_results")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ====================== 2. File lists ======================
gwas_files <- list.files(
  GWAS_DIR,
  pattern = "^assoc\\.regenie\\.merged_.+\\.txt$",
  full.names = TRUE
)

outcome_files <- list.files(
  OUTCOME_DIR,
  pattern = "^finngen_R12_J10_",
  full.names = TRUE
)

cat("Number of outcome files:", length(outcome_files), "\n")
print(head(basename(outcome_files), 10))

cat(
  "Found", length(gwas_files), "exposure files and",
  length(outcome_files), "outcome files\n"
)

# ====================== 3. Main MR function ======================
process_mr_analysis <- function(gwas_files, outcome_files) {
  total_combinations <- length(gwas_files) * length(outcome_files)
  current_progress <- 0
  
  results <- map(gwas_files, function(gwas_file) {
    file_tag <- gsub("^assoc\\.regenie\\.merged_|\\.txt$", "", basename(gwas_file))
    ld_file <- file.path(LD_DIR, paste0("plink_all_ld_clumped_", file_tag, ".clumped"))
    
    if (!file.exists(ld_file)) {
      message("LD file not found: ", ld_file)
      return(NULL)
    }
    
    map(outcome_files, function(outcome_path) {
      current_progress <<- current_progress + 1
      progress_percent <- round((current_progress / total_combinations) * 100, 1)
      
      cat(sprintf(
        "[%d/%d %.1f%%] Processing: %s vs %s\n",
        current_progress,
        total_combinations,
        progress_percent,
        file_tag,
        basename(outcome_path)
      ))
      
      outcome_tag <- sub("^finngen_R12_J10_", "", basename(outcome_path))
      result_dir <- file.path(OUTPUT_DIR, paste0(file_tag, "_vs_", outcome_tag))
      dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
      
      # ------------------ Exposure preparation ------------------
      exposure_gwas <- tryCatch({
        exp_data <- fread(gwas_file, header = TRUE, stringsAsFactors = FALSE)
        ld_data <- fread(ld_file, header = TRUE, stringsAsFactors = FALSE)
        iv_snps <- ld_data$SNP
        
        exp_data %>%
          filter(ID %in% iv_snps) %>%
          mutate(
            pval = 10^(-LOG10P),
            effect_allele = as.character(ALLELE1),
            other_allele = as.character(ALLELE0)
          ) %>%
          select(
            SNP = ID,
            effect_allele,
            other_allele,
            beta = BETA,
            se = SE,
            pval,
            eaf = A1FREQ
          )
      }, error = function(e) {
        message("Exposure read error: ", conditionMessage(e))
        NULL
      })
      
      if (is.null(exposure_gwas) || nrow(exposure_gwas) == 0) {
        message("  -> No valid exposure SNPs")
        return(NULL)
      }
      
      # ------------------ Outcome preparation ------------------
      outcome_gwas <- tryCatch({
        out_data <- fread(outcome_path, header = TRUE, data.table = FALSE)
        
        out_data %>%
          filter(rsids %in% exposure_gwas$SNP) %>%
          mutate(
            pval = as.numeric(pval),
            beta = as.numeric(beta),
            sebeta = as.numeric(sebeta),
            eaf = as.numeric(af_alt),
            effect_allele = as.character(alt),
            other_allele = as.character(ref)
          ) %>%
          select(
            SNP = rsids,
            effect_allele,
            other_allele,
            beta,
            se = sebeta,
            pval,
            eaf
          )
      }, error = function(e) {
        message("Outcome read error: ", conditionMessage(e))
        NULL
      })
      
      if (is.null(outcome_gwas) || nrow(outcome_gwas) == 0) {
        message("  -> No matched outcome SNPs")
        return(NULL)
      }
      
      write.csv(
        exposure_gwas,
        file.path(result_dir, "exposure_mr_format.csv"),
        row.names = FALSE,
        quote = FALSE
      )
      write.csv(
        outcome_gwas,
        file.path(result_dir, "outcome_mr_format.csv"),
        row.names = FALSE,
        quote = FALSE
      )
      
      # ------------------ MR analysis ------------------
      tryCatch({
        exposure_dat <- read_exposure_data(
          file.path(result_dir, "exposure_mr_format.csv"),
          sep = ",",
          snp_col = "SNP",
          effect_allele_col = "effect_allele",
          other_allele_col = "other_allele",
          eaf_col = "eaf",
          beta_col = "beta",
          se_col = "se",
          pval_col = "pval"
        )
        
        outcome_dat <- read_outcome_data(
          file.path(result_dir, "outcome_mr_format.csv"),
          sep = ",",
          snp_col = "SNP",
          effect_allele_col = "effect_allele",
          other_allele_col = "other_allele",
          eaf_col = "eaf",
          beta_col = "beta",
          se_col = "se",
          pval_col = "pval"
        )
        
        dat <- harmonise_data(exposure_dat, outcome_dat)
        
        if (nrow(dat) == 0) {
          message("  -> No valid SNPs available for MR analysis")
          return(NULL)
        }
        
        mr_results <- mr(
          dat,
          method_list = c(
            "mr_ivw",
            "mr_ivw_mre",
            "mr_weighted_median",
            "mr_egger_regression"
          )
        )
        
        heterogeneity <- mr_heterogeneity(dat)
        pleiotropy <- mr_pleiotropy_test(dat)
        
        write.csv(mr_results, file.path(result_dir, "MR_main_results.csv"), row.names = FALSE)
        write.csv(heterogeneity, file.path(result_dir, "heterogeneity_results.csv"), row.names = FALSE)
        write.csv(pleiotropy, file.path(result_dir, "pleiotropy_results.csv"), row.names = FALSE)
        write.csv(dat, file.path(result_dir, "harmonised_data.csv"), row.names = FALSE)
        
        cat("  -> MR analysis completed with", nrow(dat), "harmonised SNPs\n")
        
        list(
          exposure = file_tag,
          outcome = outcome_tag,
          mr_results = mr_results,
          heterogeneity = heterogeneity,
          pleiotropy = pleiotropy,
          harmonised_snps = nrow(dat)
        )
      }, error = function(e) {
        message("MR analysis failed: ", conditionMessage(e))
        NULL
      })
    }) %>%
      compact()
  }) %>%
    flatten() %>%
    compact()
  
  results
}

# ====================== 4. Run analysis ======================
cat("Starting MR analysis...\n")
final_results <- process_mr_analysis(gwas_files, outcome_files)

# ====================== 5. Summary output ======================
if (length(final_results) > 0) {
  cat("\nSummarising results...\n")
  
  summary_df <- map_dfr(final_results, ~{
    mr_res <- .x$mr_results
    het_res <- .x$heterogeneity
    pleio_res <- .x$pleiotropy
    
    ivw_row <- which(mr_res$method == "Inverse variance weighted")
    ivw_mre_row <- which(mr_res$method == "Inverse variance weighted (multiplicative random effects)")
    wm_row <- which(mr_res$method == "Weighted median")
    egger_row <- which(mr_res$method == "MR Egger")
    het_row <- which(het_res$method == "Inverse variance weighted")
    
    data.frame(
      Exposure = .x$exposure,
      Outcome = .x$outcome,
      NSnps = .x$harmonised_snps,
      
      IVW_beta = if (length(ivw_row) > 0) mr_res$b[ivw_row] else NA,
      IVW_se = if (length(ivw_row) > 0) mr_res$se[ivw_row] else NA,
      IVW_pval = if (length(ivw_row) > 0) mr_res$pval[ivw_row] else NA,
      
      IVW_MRE_beta = if (length(ivw_mre_row) > 0) mr_res$b[ivw_mre_row] else NA,
      IVW_MRE_se = if (length(ivw_mre_row) > 0) mr_res$se[ivw_mre_row] else NA,
      IVW_MRE_pval = if (length(ivw_mre_row) > 0) mr_res$pval[ivw_mre_row] else NA,
      
      WM_beta = if (length(wm_row) > 0) mr_res$b[wm_row] else NA,
      WM_se = if (length(wm_row) > 0) mr_res$se[wm_row] else NA,
      WM_pval = if (length(wm_row) > 0) mr_res$pval[wm_row] else NA,
      
      Egger_beta = if (length(egger_row) > 0) mr_res$b[egger_row] else NA,
      Egger_se = if (length(egger_row) > 0) mr_res$se[egger_row] else NA,
      Egger_pval = if (length(egger_row) > 0) mr_res$pval[egger_row] else NA,
      
      Q_pval = if (length(het_row) > 0) het_res$Q_pval[het_row] else NA,
      Egger_intercept = if (nrow(pleio_res) > 0) pleio_res$egger_intercept else NA,
      Egger_intercept_pval = if (nrow(pleio_res) > 0) pleio_res$pval else NA,
      
      stringsAsFactors = FALSE
    )
  })
  
  write.csv(summary_df, file.path(OUTPUT_DIR, "summary_results.csv"), row.names = FALSE)
  
  cat("\nAll analyses completed.\n")
  cat("Completed", nrow(summary_df), "MR analyses successfully\n")
  cat("Summary results saved to:", file.path(OUTPUT_DIR, "summary_results.csv"), "\n")
  
  cat("\nBrief summary:\n")
  cat("- Mean number of harmonised SNPs:", round(mean(summary_df$NSnps, na.rm = TRUE), 1), "\n")
  cat("- Number of significant IVW results (p < 0.05):", sum(summary_df$IVW_pval < 0.05, na.rm = TRUE), "\n")
  
  print(summary_df)
} else {
  cat("\nNo MR analysis completed successfully.\n")
}