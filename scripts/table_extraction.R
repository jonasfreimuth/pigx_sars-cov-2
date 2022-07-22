library(magrittr)
library(stringr)
library(dplyr)

get_mutations_counts <- function(mutation_plot_data,
                                 mutation_sheet_df,
                                 sign_mut_vec) {
  #' input:
  #'   * mutation_plot_data: data_mut_plot.csv df
  #'   * mutation_sheet_df: the mutation sheet with NAs at empty cells
  #'   * sign_mut_vec: A vector of mutation strings which showed a stat.
  #'   significant increase in proportion.
  #'
  #' output:
  #'   A dataframe containing counts of mutations per sample and in total.
  #'   FIXME Be more precise about what the cols are

  .get_mutation_counts_row <- function(sample_row,
                                       mutation_sheet_df,
                                       sign_incr_muts,
                                       meta_cols_excl) {
    #' input:
    #'   * sample_row: a row of data_mut_plot.csv df
    #'   * mutation_sheet_df: the mutation sheet with NAs at empty cells
    #'   * sign_incr_muts: A vector of mutation strings which showed a stat.
    #'   significant increase in proportion.
    #'   * meta_cols_excl: A vector containing the names of columns which
    #'   contain metadata instead of mutation proportions and are to be
    #'   excluded.
    #'
    #' output:
    #'   A dataframe containing counts of mutations in this row.
    #'   FIXME Be more precise about what the cols are

    sigmut_vec_all <- mutation_sheet_df %>%
      unlist(use.names = FALSE) %>%
      unique() %>%
      na.omit()

    mutations <- sample_row %>%
      as.matrix() %>%
      t() %>%
      as_tibble() %>%
      dplyr::select(-all_of(meta_cols_excl)) %>%
      dplyr::select(-where(is.na)) %>%
      names() %>%

      # extract only nucleotide mutation part, drop AA muts
      str_extract("[A-Z0-9*_]+$")

    counts_tot_sample <- data.frame(
      sample                = sample_row["samplename"],
      total_muts            = length(mutations),
      total_sigmuts         = sum(mutations %in% sigmut_vec_all),
      # get num of muts with significant increase over time
      tracked_muts_after_lm = sum(mutations %in% sign_incr_muts)
    ) %>%

      # get number of mutations which aren't signature mutations
      mutate(non_sigmuts = total_muts - total_sigmuts)

    counts_var_sample <- lapply(
      mutation_sheet_df,
      function(var_muts, mutations) {
        var_muts <- na.omit(var_muts)
        return(sum(mutations %in% var_muts))
      }, mutations) %>%
      bind_cols() %>%
      set_names(paste0("sigmuts_", names(.)))

    return(bind_cols(counts_tot_sample, counts_var_sample))
  }

  # transform mutation_sheet to one comparable vector
  sigmut_vec_all <- mutation_sheet_df %>%
    unlist(use.names = FALSE) %>%
    unique() %>%
    na.omit()

  # get vector of significantly increasing mutations
  sign_incr_muts <- mutations_sig$mutation %>%

    # extract only nucleotide mutation part, drop AA muts
    str_extract("[A-Z0-9*_]+$")

  # create vector of metadata col names to be excluded
  meta_cols_excl <- c(
    "samplename",
    "location_name",
    "coordinates_lat",
    "coordinates_long",
    "dates"
  )

  # get names of mutations without meta data
  mutations <- mutation_plot_data %>%
    select(-all_of(meta_cols_excl)) %>%
    names() %>%

    # extract only nucleotide mutation part, drop AA muts
    str_extract("[A-Z0-9*_]+$")

  # signature mutations found across samples
  common_sigmuts_df <- mutation_plot_data %>%
    dplyr::select(dplyr::contains(sigmut_vec_all))

  if (length(mutations_sig) > 0) {
    # total counts across all samples, without duplicated counts
    # (compared to what I'd get when summing up all the columns)
    counts_tot_all <- data.frame(
      sample = "Total",
      total_muts = length(mutations),
      total_sigmuts = ncol(common_sigmuts_df),

      # FIXME Is this not just the length of sign_incr_muts?
      tracked_muts_after_lm = sum(mutations %in% sign_incr_muts)
    ) %>%

      # get number of mutations which aren't signature mutations
      mutate(non_sigmuts = total_muts - total_sigmuts)

    counts_var_all <- lapply(
      mutation_sheet_df,
      function(var_muts, mutations) {
        var_muts <- na.omit(var_muts)
        return(sum(mutations %in% var_muts))
      }, mutations) %>%
      bind_cols() %>%
      set_names(paste0("sigmuts_", names(.)))

    counts_all <- bind_cols(
      counts_tot_all,
      counts_var_all
    )

    counts_per_sample <- apply(
      mutation_plot_data,
      1,
      .get_mutation_counts_row,
      mutation_sheet_df,
      sign_incr_muts,
      meta_cols_excl
    ) %>%
      bind_rows()

    count_frame <- bind_rows(counts_all, counts_per_sample)
  } else {
    count_frame <- data.frame()
  }

  return(count_frame)
}
