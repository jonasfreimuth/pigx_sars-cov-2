library(dplyr)

count_muts <- function(sample_row, mutation_sheet_df) {
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
  mutations_ps <- sample_row[
    (which(
      names(sample_row) %in% "coordinates_long"
    ) + 1):length(names(sample_row))
  ]
  count_frame <- data.frame(
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
            dplyr::select(dplyr::contains(mutations_sig$mutation))
        )
      )
    )
  )
  # get number of mutations which aren't signature mutations
  count_frame <- count_frame %>%
    mutate(non_sigmuts = total_muts - total_sigmuts)
  # get number of siganture mutation per variant
  for (var in colnames(mutation_sheet_df)) {
    count_frame[, paste0("sigmuts_", var)] <- as.numeric(
      rowSums(
        !is.na(
          mutations_ps %>%
            dplyr::select(dplyr::contains(na.omit(mutation_sheet_df[[var]])))
        )
      )
    )
  }
  return(count_frame)
}

write_mutations_count <- function(mutation_plot_data,
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

  # total counts across all samples, without duplicated counts
  # (compared to what I'd get when summing up all the columns)
  if (length(mutations_sig) > 0) {
    count_frame <- data.frame(
      sample = "Total",
      total_muts = length(mutations),
      total_sigmuts = length(common_sigmuts_df),
      # get how many of all found mutations will be tracked because of
      # significant increase over time
      tracked_muts_after_lm = length(
        mutation_plot_data %>%
          dplyr::select(dplyr::contains(mutations_sig$mutation))
      )
    )
    # get number of mutations which aren't signature mutations
    count_frame <- count_frame %>%
      mutate(non_sigmuts = total_muts - total_sigmuts)

    # get number of siganture mutation per variant
    for (var in colnames(mutation_sheet_df)) {
      sigmut_pv <- mutation_plot_data %>%
        dplyr::select(dplyr::contains(na.omit(mutation_sheet_df[[var]])))
      count_frame[, paste0("sigmuts_", var)] <- length(sigmut_pv)
    }

    counts_per_sample <- do.call(
      bind_rows,
      apply(mutation_plot_data, 1, count_muts, mutation_sheet_df)
    )
    count_frame <- bind_rows(count_frame, counts_per_sample)
  } else {
    count_frame <- data.frame()
  }

  return(count_frame)
}
