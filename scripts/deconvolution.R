## ----setup, include = FALSE, warning = FALSE----------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)

library(knitr)
library(dplyr)
library(ggplot2)
library(DT)
library(tidyr)
library(qpcR)
library(stringr)
library(magrittr)
library(base64url)

## command line arguments
args <- commandArgs(trailingOnly = TRUE)

# give default parameters
if (length(args) == 0) {
  args <- c(
    sample_name = "",
    output_dir = "",
    vep_file = "",
    snv_file = "",
    sample_sheet = "",
    mutation_sheet = "",
    deconvolution_functions = ""
  )
}

params <- list(
  sample_name = args[[1]],
  output_dir = args[[2]],
  vep_file = args[[3]],
  snv_file = args[[4]],
  sample_sheet = args[[5]],
  mutation_sheet = args[[6]],
  deconvolution_functions = args[[7]]
)


## function loading
source(params$deconvolution_functions)


## ----print_input_settings, echo = FALSE---------------------------------------
sample_name         <- params$sample_name
sample_sheet        <- data.table::fread(params$sample_sheet)
mutation_sheet      <- params$mutation_sheet

csv_output_dir      <- file.path(params$output_dir, "variants")
mutation_output_dir <- file.path(params$output_dir, "mutations")
date                <- as.character(sample_sheet[name == sample_name]$date)
location_name       <- as.character(sample_sheet[name == sample_name]$location_name)
coordinates_lat     <- as.character(sample_sheet[name == sample_name]$coordinates_lat)
coordinates_long    <- as.character(sample_sheet[name == sample_name]$coordinates_long)


## ----process_signature_mutations, include = FALSE-----------------------------
# Read signature data
sig_mutations.df <- read.csv(mutation_sheet, header = TRUE)

if ("source" %in% colnames(sig_mutations.df)) {
  sig_mutations.df <- sig_mutations.df[, -(which(names(sig_mutations.df) %in%
    "source"))]
}

sig_mutations.df <- sig_mutations.df %>%
  dplyr::na_if("") %>%
  tidyr::pivot_longer(everything(), values_drop_na = TRUE)

vepfile.df <- read.table(params$vep_file, sep = ",", header = TRUE)

vepfile.df <- na_if(vepfile.df, "-")

# deduplicate dataframe
# reasoning see description of "dedupeMuts"
# FIXME: I was not sure, how I can include this in my function but should be
# possible, will come back to it later
dupes <- duplicated(sig_mutations.df$value, fromLast = TRUE)
allDupes <- dupes | duplicated(sig_mutations.df$value, fromLast = FALSE)

if (any(allDupes)) {
  sigmuts_deduped <- sig_mutations.df[!dupes, ]
}

# FIXME: I think this should be done with some apply() function?
for (mut in sig_mutations.df$value) {
  sigmuts_deduped <- dedupeMuts(mut, sig_mutations.df, sigmuts_deduped)
}

sigmuts_deduped_no_gene <- sigmuts_deduped %>%
  rowwise() %>%
  mutate(mutation = str_split(value, ":")[[1]][2]) %>%
  dplyr::select(name, mutation)


## ----match_snvs_to_signature_mutations, include = FALSE-----------------------
variant_protein_mut <- get_protein_mut(params$vep_file)
# match the variant names according to the signature mutations, if there is no
# signature mutation no name will be given
# variant characterizing is done by NT mutations
variant_protein_mut <- dplyr::left_join(variant_protein_mut,
  sigmuts_deduped_no_gene,
  by = c("gene_mut" = "mutation")
)


## ----merge_vep_with_lofreq_info, include = FALSE------------------------------
# get the SNV frequency values and coverage information for the mutations from
# the LoFreq output
lofreq.info <- as_data_frame(parse_snv_csv(params$snv_file))
vep.info <- variant_protein_mut

complete.df <- dplyr::left_join(lofreq.info, vep.info,
  by = c("gene_mut" = "gene_mut"), copy = TRUE
) %>%
  rowwise() %>%
  mutate(gene_mut_collapsed = paste(genes, gene_mut, sep = ":"))

# filter for mutations which are signature mutations
match.df <- complete.df %>%
  filter(!is.na(name))

# filter for everything that is not a signature mutation
nomatch.df <- complete.df %>%
  filter(is.na(name))

write.csv(match.df,
  file.path(mutation_output_dir,
  paste0(sample_name,
    "_sigmuts.csv")))

write.csv(nomatch.df,
  file.path(mutation_output_dir,
  paste0(sample_name,
    "_non_sigmuts.csv")))

# Tables are displayed here in report


## ----getting_unique_muts_bulk, include = FALSE--------------------------------
# get  NT mutations only, input for the signature matrix
mutations.vector <- match.df$gene_mut_collapsed
# get bulk frequency values, will be input for the deconvolution function
bulk_freq.vector <- as.numeric(match.df$freq)

# only execute the deconvolution when at least one signature mutation was found
executeDeconvolution <- length(mutations.vector) > 0


if (executeDeconvolution) {
  ## ----creating_signature_matrix, include = FALSE-------------------------------
  # create an empty data frame add a column for the Wildtype
  # Wildtype in this case means the reference version of SARS-Cov-2
  # for the deconvolution to work we need the "wild type" frequencies too.
  # The matrix from above got mirrored, wild type mutations are simulated the
  # following: e.g. T210I (mutation) -> T210T ("wild type")
  msig_simple <- createSigMatrix(mutations.vector, mutation_sheet)

  # When multiple columns look like the same, the deconvolution will not work,
  # because the function can't distinguish between those columns. The workaround
  # for now is to identify those equal columns and merge them into one, returning
  # also a vector with the information about which of the columns were merged.
  msig_simple <- cbind(muts = mutations.vector, msig_simple)
  msig_transposed <- dedupeDF(msig_simple)
  msig_stable_transposed <- msig_transposed[[1]]
  msig_dedupe_transposed <- msig_transposed[[2]]

  dropped_variants <- c()

  # for every variant update the rownames with the group they are in
  # FIXME: Shorten this and similar constructs
  for (variant in rownames(msig_stable_transposed[-(rownames(msig_stable_transposed) %in% "muts"), ])) {
    grouping_res <- dedupeVariants(
      variant,
      msig_stable_transposed,
      msig_dedupe_transposed
    )

    msig_dedupe_transposed <- grouping_res[[1]]
    dropped_variants <- c(dropped_variants, grouping_res[[2]])
  }

  # transpose the data frame back to column format for additional processing
  if (length(msig_dedupe_transposed) >= 1) {
    # the 1 get's rid of the additional first row which is an transposing artifact
    msig_simple_unique <- as.data.frame(t(msig_dedupe_transposed[, -1])) %>%
      mutate(across(!c("muts"), as.numeric))
  }

  # clean the vector to know which variants has to be add with value 0 after deconvolution
  dropped_variants <- unique(dropped_variants)
  dropped_variants <- dropped_variants[!is.na(dropped_variants)]


  ## ----calculate_sigmat_weigths, include = FALSE--------------------------------
  deconv_lineages <- colnames(msig_simple_unique[, -which(names(msig_simple_unique) %in% c("muts", "WT"))])

  # create list of proportion values that will be used as weigths
  sigmut_proportion_weights <- list()


  for (lineage in deconv_lineages) {
    if (grepl(",", lineage)) {
      group <- unlist(str_split(lineage, ","))
      avrg <- sum(sig_mutations.df$name %in% group) / length(group)
      value <- sum(msig_simple_unique[lineage]) / avrg
    } else {
      value <- sum(msig_simple_unique[lineage]) / sum(sig_mutations.df$name == lineage)
    }
    sigmut_proportion_weights[lineage] <- value
  }

  sigmut_proportion_weights <- as_tibble(sigmut_proportion_weights)

  # applying weights on signature matrix
  # FIXME: there should be a way to do this vectorized
  msig_simple_unique_weighted <- msig_simple_unique
  for (lineage in deconv_lineages) {
    msig_simple_unique_weighted[lineage] <- msig_simple_unique_weighted[lineage] / as.numeric(sigmut_proportion_weights[lineage])
  }


  ## ----simulating_WT_mutations, include = FALSE---------------------------------
  # construct additional WT mutations that are not weighted
  msig_stable_all <- simulateWT(
    mutations.vector, bulk_freq.vector,
    msig_simple_unique_weighted[, -which(names(msig_simple_unique_weighted) == "muts")],
    match.df$cov
  )

  msig_stable_unique <- msig_stable_all[[1]]


  ## ----deconvolution, include = FALSE-------------------------------------------
  # this hack is necessary because otherwise the deconvolution will throw:
  # Error in x * wts: non-numeric argument to binary operator
  # also see: https://stackoverflow.com/questions/37707060/converting-data-frame-column-from-character-to-numeric/37707117
  sig <- apply(
    msig_stable_unique[, -which(names(msig_stable_unique) %in% "muts")],
    2,
    function(x) {
      as.numeric(as.character(x))
    }
  )

  bulk_all <- as.numeric(msig_stable_all[[2]])
  variant_abundance <- deconv(bulk_all, sig)


  ## ----plot, echo = FALSE-------------------------------------------------------
  # work in progress...only to show how it theoretically can look like in the
  # report
  variants <- colnames(msig_stable_unique[, -1])
  df <- data.frame(rbind(variant_abundance))

  # TODO: Replace the long column name consisting of concatenated variant
  # names with "others" and put the variant names in a label.
  # TODO: This should be done much earlier when building the data frame.
  condensed_variants_names <- unlist(
    lapply(
      variants,
      function(x) {
        if (str_detect(x, ".*,.*,.*")) {
          "others"
        } else {
          x
        }
      }
    )
  )

  colnames(df) <- condensed_variants_names

  variants_labels <- unlist(
    lapply(
      variants,
      function(x) str_replace_all(x, ",", "\n")
    )
  )

  df <- df %>%
    tidyr::pivot_longer(everything())

  # Handling of ambiguous cases and grouped variants

  # case 1: add dropped variants again with value 0 in case all of the other
  # variants add up to 1
  if (round(sum(df$value), 1) == 1) {
    for (variant in dropped_variants) {
      df <- rbind(df, c(variant, 0))
    }
  }

  # case 2: in case "others" == 0, both variants can be split up again and being
  # given the value 0 OR case 3: in case multiple vars can really not be
  # distinguished from each other they will be distributed normaly
  if (any(str_detect(variants, ","))) {
    grouped_rows <- which(str_detect(variants, ","))
    for (row in grouped_rows) {
      if (df[row, "value"] == 0) {
        grouped_variants <- unlist(str_split(df[row, "name"], ","))
        for (variant in grouped_variants) {
          # add new rows, one for each variant
          df <- rbind(df, c(variant, 0))
        }
      } else if (df[row, "value"] != 0) {
        grouped_variants <- unlist(str_split(df[row, "name"], ","))
        # normal distribution, devide deconv value by number of grouped variants
        distributed_freq_value <-
          as.numeric(as.numeric(df[row, "value"]) / length(grouped_variants))
        for (variant in grouped_variants) {
          # add new rows, one for each variant
          df <- rbind(df, c(variant, distributed_freq_value))
        }
      }
      # drop grouped row
      df <- df[-row, ]
    }
  }

  df <- transform(df, value = as.numeric(value))

  write.csv(df, file.path(
    mutation_output_dir,
    paste0(
      sample_name,
      "variant_proportion_barplot_data.csv"
    )
  ))

  # plot comes here in report


  ## ----csv_output_variant_plot, include = FALSE---------------------------------
  # prepare processed variant values to output them as a csv which will be used
  # for the plots in index.rmd those outputs are not offically declared as outputs
  # which can lead to issues - that part should be handled by a seperate file
  # (and maybe rule) get all possible variants
  all_variants <- colnames(msig_simple[, -which(names(msig_simple) %in% "muts")])
  output_variants <- file.path(csv_output_dir, "data_variant_plot.csv")

  # 1. write dataframe with this information here
  output_variant_plot <- data.frame(
    samplename = character(),
    dates = character(),
    location_name = character(),
    coordinates_lat = character(),
    coordinates_long = character()
  )

  # add columns for all possible variants to the dataframe
  for (variant in all_variants) {
    output_variant_plot[, variant] <- numeric()
  }

  meta_data <- c(
    samplename = sample_name,
    dates = date,
    location_name = location_name,
    coordinates_lat = coordinates_lat,
    coordinates_long = coordinates_long
  )

  output_variant_plot <- bind_rows(output_variant_plot, meta_data)

  # get rownumber for current sample
  sample_row <- which(grepl(sample_name, output_variant_plot$samplename))

  # write mutation frequency values to df
  for (i in all_variants) {
    if (i %in% df$name) {
      # check if variant already has a column
      if (i %in% colnames(output_variant_plot)) {
        output_variant_plot[sample_row, ][i] <- df$value[df$name == i]
        output_variant_plot <- output_variant_plot %>%
          mutate(others = 1 - rowSums(across(all_of(all_variants)), na.rm = TRUE))
      }
    }
  }

  # 2. check if file exists already
  if (file.exists(output_variants)) {
    previous_df <- read.table(
      output_variants,
      sep = "\t",
      header = TRUE,
      colClasses = "character",
      check.names = FALSE
    )
    # convert numeric values to character
    output_variant_plot <-
      as.data.frame(lapply(output_variant_plot, as.character), check.names = FALSE)
    # merge with adding cols and rows
    output_variant_plot <- full_join(previous_df,
      output_variant_plot,
      by = colnames(previous_df),
      copy = TRUE
    )
  }

  # 3. write to output file
  write.table(output_variant_plot, output_variants,
    sep = "\t",
    na = "NA", row.names = FALSE, quote = FALSE
  )


  ## ----csv_output_mutation_plot, include = FALSE--------------------------------
  # prepare processed mutation values to output them as a csv which will be used
  # for the plots in index.rmd those outputs are not officially declared as
  # outputs which can lead to issues - that part should be handled by a seperate
  # file (and maybe rule)
  # get all possible mutations

  # one aa mutation can have different codon mutations reported with
  # different freqs- for the summary table they have to be summed up
  # (process see line 1872 of documentation)
  complete.df <- complete.df %>%
    group_by(across(c(-freq, -gene_mut, -gene_mut_collapsed, AA_mut))) %>%
    summarise(
      freq = sum(as.numeric(freq)),
      gene_mut = paste(gene_mut, collapse = ",")
    ) %>%
    rowwise() %>%
    mutate(AA_mut = replace(AA_mut, is.na(AA_mut), "\\:\\")) %>%
    # 211006 this exclusion is necessary because this mutation has a wrong entry
    # in VEP which gives two AA_muts instead of probably 1 deletion
    filter(!(gene_mut %in% "G13477A")) %>%
    ungroup()

  # report the gene, translated_AA_mut and NT mut accordingly
  # easier to spot translation inconsitentcies that way
  all_mutations <- paste(complete.df$AA_mut[!is.na(complete.df$AA_mut)],
    complete.df$gene_mut,
    sep = "::"
  )

  output_mutations <- file.path(
    mutation_output_dir,
    paste0(sample_name, "_mutations.csv")
  )

  # 1. write dataframe with this information here
  output_mutation_frame <- data.frame(
    samplename = character(),
    dates = character(),
    location_name = character(),
    coordinates_lat = character(),
    coordinates_long = character()
  )
  # add columns for all possible mutations to the dataframe
  for (mutation in all_mutations) {
    output_mutation_frame[, mutation] <- numeric()
  }
  meta_data <- c(
    samplename = sample_name,
    dates = date,
    location_name = location_name,
    coordinates_lat = coordinates_lat,
    coordinates_long = coordinates_long
  )

  output_mutation_frame <- bind_rows(output_mutation_frame, meta_data)

  # write mutation frequency values to df
  for (i in all_mutations) {
    i_NT <- str_split(i, "::")[[1]][2]
    if (i_NT %in% complete.df$gene_mut) { # split gene name to match with AA mut
      # check if variant already has a column
      if (i %in% colnames(output_mutation_frame)) {
        output_mutation_frame[, i] <- complete.df$freq[which(
          complete.df$gene_mut == i_NT
        )]
      }
    }
  }
  colnames(output_mutation_frame) <- as.character(colnames(output_mutation_frame))
  output_mutation_frame <- output_mutation_frame %>%
    dplyr::select(-contains("NA", ignore.case = FALSE))

  # convert numeric values to character
  output_mutation_frame <- as.data.frame(lapply(
    output_mutation_frame,
    as.character
  ),
  check.names = FALSE
  )

  # 3. write to output file
  write.table(output_mutation_frame, output_mutations,
    sep = "\t",
    row.names = FALSE, quote = FALSE
  )
}
