library(stringr)
library(magrittr)
library(dplyr)

get_mutations_counts <- function(mutation_plot_data,
                                 mutation_sheet_df,
                                 mutations_sig) {
  #' input:
  #'   * mutation_plot_data: data_mut_plot.csv df
  #'   * mutation_sheet_df: the mutation sheet with NAs at empty cells
  #'   * mutations_sig: TODO
  #'
  #' output:
  #'   A dataframe containing counts of mutations per sample and in total.
  #'   FIXME Be more precise about what the cols are

  count_muts <- function(sample_row, mutation_sheet_df, sign_mut_vec) {
    #' function used in rowwise apply() call
    #' takes row as input, calculates mutation counts and returns a dataframe
    #'
    # transform mutation_sheet to one comparable vector
    sigmut_vec_all <- mutation_sheet_df %>%
      unlist(use.names = FALSE) %>%
      unique() %>%
      na.omit()

    # create vector of metadata col names to be excluded
    meta_cols_excl <- c(
      "samplename",
      "location_name",
      "coordinates_lat",
      "coordinates_long",
      "dates"
    )

    sample_mut_row <- sample_row[-which(names(sample_row) %in% meta_cols_excl)]
    mut_vec_sample <- names(sample_mut_row)[!is.na(sample_mut_row)]

    counts_tot_sample <- data_frame(
      sample                = sample_row["samplename"],
      total_muts            = length(mut_vec_sample),
      total_sigmuts         = sum(
        str_extract(mut_vec_sample, "[A-Z0-9*_]+$") %in% sigmut_vec_all
      ),
      tracked_muts_after_lm = sum(mut_vec_sample %in% sign_mut_vec)
    ) %>%
      mutate(non_sigmuts = total_muts - total_sigmuts)

    counts_var_sample <- lapply(
      colnames(mutation_sheet_df),
      function(var) {
        var_muts <- na.omit(mutation_sheet_df[[var]])
        count <- sum(str_extract(mut_vec_sample, "[A-Z0-9*_]+$") %in% var_muts)
        return(
          data.frame(count) %>%
            set_names(paste0("sigmuts_", var))
        )
      }) %>%
      bind_cols()


    return(bind_cols(counts_tot_sample, counts_var_sample))
  }

  # transform mutation_sheet to one comparable vector
  sigmut_vec_all <- mutation_sheet_df %>%
    unlist(use.names = FALSE) %>%
    unique() %>%
    na.omit()

  # create vector of metadata col names to be excluded
  meta_cols_excl <- c(
    "samplename",
    "location_name",
    "coordinates_lat",
    "coordinates_long",
    "dates"
  )

  # get names of mutations without meta data
  mutations <- names(mutation_plot_data %>% select(-all_of(meta_cols_excl)))

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

      # get how many of all found mutations will be tracked because of
      # significant increase over time
      tracked_muts_after_lm = ncol(
        mutation_plot_data %>%
          dplyr::select(all_of(mutations_sig$mutation))
      )) %>%
      # get number of mutations which aren't signature mutations
      mutate(non_sigmuts = total_muts - total_sigmuts)

    counts_var_all <- lapply(
      colnames(mutation_sheet_df),
      function(var) {
        var_muts <- na.omit(mutation_sheet_df[[var]])
        count <- sum(str_extract(mutations, "[A-Z0-9*_]+$") %in% var_muts)
        return(
          data.frame(count) %>%
            set_names(paste0("sigmuts_", var))
        )
      }) %>%
      bind_cols()

    counts_all <- bind_cols(counts_tot_all, counts_var_all)


    # get all the per sample counts
    counts_sample <- apply(
      mutation_plot_data,
      1,
      count_muts,
      mutation_sheet_df,
      mutations_sig$mutation
    ) %>%
      bind_rows()

    count_frame <- bind_rows(counts_all, counts_sample)
  } else {
    count_frame <- data.frame()
  }

  return(count_frame)
}
