library(dplyr)

count_muts <- function(sample_row, mutation_sheet_df, sign_incr_muts) {
  #' function used in rowwise apply() call
  #' takes row as input, calculates mutation counts and returns a dataframe
  #'
  # transform mutation_sheet to one comparable vector
  sigmut_vec_all <- mutation_sheet_df %>%
    unlist(use.names = FALSE) %>%
    unique() %>%
    na.omit()

  # transform char. vector into dataframe
  sample_row <- as_tibble(t(as.matrix(sample_row)))

  # create vector of metadata col names to be excluded
  meta_cols_excl <- c(
    "samplename",
    "location_name",
    "coordinates_lat",
    "coordinates_long",
    "dates"
  )

  mutations_ps <- sample_row %>%
    dplyr::select(-all_of(meta_cols_excl))

  counts_tot_sample <- data.frame(
    # mutations only
    sample = as.character(sample_row["samplename"]),
    # count all mutations which are not NA
    total_muts = as.numeric(rowSums(!is.na(mutations_ps))),
    # count all mutations that are signature mutations
    total_sigmuts = as.numeric(
      rowSums(
        !is.na(
          mutations_ps %>%
            dplyr::select(dplyr::contains(sigmut_vec_all))
        )
      )
    ),
    # get num of muts with significant increase over time
    tracked_muts_after_lm = as.numeric(
      rowSums(
        !is.na(
          mutations_ps %>%
            dplyr::select(dplyr::contains(sign_incr_muts))
        )
      )
    )
  ) %>% 

    # get number of mutations which aren't signature mutations
    mutate(non_sigmuts = total_muts - total_sigmuts)


  counts_var_sample <- counts_tot_sample %>%
    dplyr::select(sample)

  # get number of siganture mutation per variant
  for (var in colnames(mutation_sheet_df)) {
    counts_var_sample[, paste0("sigmuts_", var)] <- as.numeric(
      rowSums(
        !is.na(
          mutations_ps %>%
            dplyr::select(dplyr::contains(na.omit(mutation_sheet_df[[var]])))
        )
      )
    )
  }
  return(dplyr::left_join(counts_tot_sample, counts_var_sample, by = "sample"))
}

get_mutations_counts <- function(mutation_plot_data,
                                  mutation_sheet_df,
                                  mutations_sig) {
  #' takes data_mut_plot.csv df, mutation_sheet_df with NAs at empty cells,
  #' mutations_sig.df as input
  #' counts mutations and return them as a dataframe

  # transform mutation_sheet to one comparable vector
  sigmut_vec_all <- mutation_sheet_df %>%
    unlist(use.names = FALSE) %>%
    unique() %>%
    na.omit()

  # get vector of significantly increasing mutations
  sign_incr_muts <- mutations_sig$mutation

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
    names()

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
          dplyr::select(all_of(sign_incr_muts))
      )
    ) %>%

      # get number of mutations which aren't signature mutations
      mutate(non_sigmuts = total_muts - total_sigmuts)


    counts_var_all <- counts_tot_all %>%
    dplyr::select(sample)

    # get number of siganture mutation per variant
    for (var in colnames(mutation_sheet_df)) {
      sigmut_pv <- mutation_plot_data %>%
        dplyr::select(dplyr::contains(na.omit(mutation_sheet_df[[var]])))
      counts_var_all[, paste0("sigmuts_", var)] <- length(sigmut_pv)
    }

    counts_all <- left_join(
      counts_tot_all,
      counts_var_all,
      by = "sample"
    )

    counts_per_sample <- apply(
      mutation_plot_data,
      1,
      count_muts,
      mutation_sheet_df,
      sign_incr_muts
    ) %>%
    bind_rows()

    count_frame <- bind_rows(counts_all, counts_per_sample)
  } else {
    count_frame <- data.frame()
  }

  return(count_frame)
}
