params <- list(
  variants_csv = "",
  mutations_csv = "",
  coverage_dir = "",
  mutation_sheet = "",
  fun_cvrg_scr = "",
  fun_lm = "",
  fun_pool = "",
  overviewQC = "",
  output_dir = "")


## ----libraries----------------------------------------------------------------

library(dplyr)

## ----input--------------------------------------------------------------------
df_var <- read.csv(params$variants_csv,
  header = TRUE,
  check.names = FALSE)
df_mut <- read.csv(params$mutations_csv,
  header = TRUE,
  check.names = FALSE)


## ----filter_plot_frames_samplescore, warning=TRUE-----------------------------
source(params$fun_pool)

# pool the samples per day, discard locations
weights <- read.csv(params$overviewQC, header = TRUE, check.names = FALSE) %>%
  dplyr::select(c(samplename, total_reads))


## ----linear_regression, eval=RUN_MUTATION_REGRESSION--------------------------
# only run regression if there are values after filtering and at least two dates
# are left over
if (nrow(approved_mut_plot) > 0 &&
  length(unique(approved_mut_plot$dates)) > 1) {
  source(params$fun_lm)

  mutation_sheet <- params$mutation_sheet

  sigmuts_df <- read.csv(mutation_sheet, header = TRUE) %>%
    na_if("") %>%
    # split gene name of for easier matching
    mutate_all(funs(str_replace(., "^[^:]*:", "")))


  changing_muts <- parsing_mutation_plot_data(approved_mut_plot)
  mutations_sig_unfiltered <- refined_lm_model(changing_muts)
  mutations_sig <- filter_lm_res_top20(mutations_sig_unfiltered, 0.05)

  ## ----mutation_counts--------------------------------------------------------
  # get functions for counting and writing
  source(params$fun_tbls)
  count_frame <- write_mutations_count(df_mut, sigmuts_df, mutations_sig)
} else {
  # TODO: This will generate an error when read as csv
  # write empty files
  count_frame <- data.frame()
  mutations_sig_unfiltered <- data.frame()
}

## file output ---------------------------------------------
output_dir <- params$output_dir

# write to file
write.csv(count_frame, file.path(output_dir, "mutations_counts.csv"),
  na = "NA", row.names = FALSE, quote = FALSE)

# mutations_sig
write.csv(mutations_sig_unfiltered,
  file.path(output_dir, "linear_regression_results.csv"),
  na = "NA", row.names = FALSE, quote = FALSE)
