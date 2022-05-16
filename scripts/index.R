params <- list(
  variants_csv = "",
  mutations_csv = "",
  coverage_dir = "",
  mutation_sheet = "",
  mutation_coverage_threshold = "90",
  fun_cvrg_scr = "",
  fun_lm = "",
  fun_pool = "",
  overviewQC = "",
  output_dir = ""
)


## ----libraries--------------------------------------------------------------------

library(dplyr)


## ----input------------------------------------------------------------------------
df_var <- read.csv(params$variants_csv,
  header = TRUE,
  check.names = FALSE
)
df_mut <- read.csv(params$mutations_csv,
  header = TRUE,
  check.names = FALSE
)

mutation_coverage_threshold <- params$mutation_coverage_threshold %>%
  # check if value is given as fraction [0,1] or percentage [0,100]
  ifelse(. >= 0 & . <= 1,
    . * 100,
    .
  ) %>%
  as.numeric()


## ----sample_quality_score---------------------------------------------------------
# TODO given as input? or only the scripts dir? similar variant_report
source(params$fun_cvrg_scr)

coverages.df <- merge(get_genome_cov(params$coverage_dir),
  get_mutation_cov(params$coverage_dir),
  by = "samplename"
) %>%
  mutate(
    proport_muts_covered = round(
      (as.numeric(total_muts_cvrd) * 100) / as.numeric(total_num_muts), 1
    )
  )

if (!(length(coverages.df) == 0)) {
  # take only the samples with > 90% mutation coverage
  good_samples.df <- coverages.df %>%
    # TODO stringecy should be given by user, maybe even through visualization
    filter(as.numeric(proport_muts_covered) >= 90)

  bad_samples.df <- coverages.df %>%
    filter(!(samplename %in% good_samples.df$samplename))
} else {
  warning("\n No coverage values found.")
}


## ----filter_plot_frames_samplescore, warning=TRUE---------------------------------
source(params$fun_pool)

# only use good samples for processing and visualiztaion
approved_var_plot <-
  df_var %>%
  filter(samplename %in% good_samples.df$samplename) %>%
  unique()
approved_mut_plot <- df_mut %>%
  filter(samplename %in% good_samples.df$samplename)

# pool the samples per day, discard locations
weights <- read.csv(params$overviewQC, header = TRUE, check.names = FALSE) %>%
  dplyr::select(c(samplename, total_reads))

# check if filtered dataframe has actual values
VARIANTS_FOUND <- nrow(approved_var_plot) > 0
MUTATIONS_FOUND <- nrow(approved_mut_plot) > 0
# only run regression if at least two dates are left over
RUN_MUTATION_REGRESSION <- if (MUTATIONS_FOUND) {
  length(unique(approved_mut_plot$dates)) > 1
} else {
  FALSE
}
# predefine to reuse later
APPROVED_MUTATIONS_FOUND <- RUN_MUTATION_REGRESSION


if (VARIANTS_FOUND) {
  approved_var_plot_location_pooled <- pool_by_weighted_mean(
    approved_var_plot,
    weights, "day_location"
  )
} else {
  warning(paste0(
    "\nNot enough variants were found or passed the quality filters.",
    "\nThe associated plots are skipped."
  ))
}


## ----linear_regression, eval=RUN_MUTATION_REGRESSION------------------------------
if (RUN_MUTATION_REGRESSION) {
  source(params$fun_lm)

  mutation_sheet <- params$mutation_sheet
  sigmuts.df <- read.csv(mutation_sheet, header = TRUE) %>%
    na_if("") %>%
    # split gene name of for easier matching
    mutate_all(funs(str_replace(., "^[^:]*:", "")))


  changing_muts <- parsing_mutation_plot_data(approved_mut_plot)
  mutations_sig_unfiltered <- refined_lm_model(changing_muts)
  mutations_sig <-
    filter_lm_res_top20(mutations_sig_unfiltered, 0.05)

  # filter good samples for only mutations with sig pvalues for plotting
  filtered_approved_mut_plot <-
    dplyr::select(
      approved_mut_plot,
      c(
        samplename,
        dates,
        location_name,
        coordinates_lat,
        coordinates_long,
        mutations_sig$mutation
      )
    )

  # check if the filtered df has actual values besides meta data
  APPROVED_MUTATIONS_FOUND <- (length(filtered_approved_mut_plot) > 5) &&
    (nrow(filtered_approved_mut_plot) > 0)

  approved_mut_plot_location_pooled <- pool_by_mean(filtered_approved_mut_plot,
    na_handling = TRUE,
    group_fun = "day_location"
  )
}
