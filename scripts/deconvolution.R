## ----setup, include = FALSE, warning = FALSE----------------------------------

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
    sample_name = "Test0",
    output_dir = "/home/jfreige/proj/pigx_sars-cov-2/tests/output",
    vep_file = "/home/jfreige/proj/pigx_sars-cov-2/tests/output/variants/Test0_vep_sarscov2_parsed.txt",
    snv_file = "/home/jfreige/proj/pigx_sars-cov-2/tests/output/variants/Test0_snv.csv",
    sample_sheet = "/home/jfreige/proj/pigx_sars-cov-2/tests/sample_sheet.csv",
    mutation_sheet = "/home/jfreige/proj/pigx_sars-cov-2/tests/sample_data/mutation_sheet_211006_covidCG_NT_location.csv",
    deconvolution_functions = "/home/jfreige/proj/pigx_sars-cov-2/scripts/deconvolution_funs.R",
    mutation_depth_threshold = "100"
  )
}

params <- list(
  sample_name = args[[1]],
  output_dir = args[[2]],
  vep_file = args[[3]],
  snv_file = args[[4]],
  sample_sheet = args[[5]],
  mutation_sheet = args[[6]],
  deconvolution_functions = args[[7]],
  mutation_depth_threshold = args[[8]]
)

# pretty print parameters for easier debugging
cat("Script running with parameters:\n\n")
par_vec <- c()
for (i in seq_along(params)) {
  par_vec[i] <- paste0(names(params)[i], " = \"", params[[i]], "\"")
}
cat(paste(par_vec, collapse = ",\n"))
cat("\n\n")

## function loading
source(params$deconvolution_functions)

# set script wide separator used for variants that contain the same signature
# muations
var_sep <- ","

## ----printInputSettings, echo = FALSE-----------------------------------------
sample_name         <- params$sample_name
sample_sheet        <- data.table::fread(params$sample_sheet)
mutation_sheet      <- params$mutation_sheet

variants_output_dir <- file.path(params$output_dir, "variants")
mutation_output_dir <- file.path(params$output_dir, "mutations")


sample_entry <- sample_sheet %>%
  filter(name == sample_name)

date                <- as.character(sample_entry$date)
location_name       <- as.character(sample_entry$location_name)
coordinates_lat     <- as.character(sample_entry$coordinates_lat)
coordinates_long    <- as.character(sample_entry$coordinates_long)


sigmut_output_file <- file.path(
  mutation_output_dir,
  paste0(
    sample_name,
    "_sigmuts.csv"
  )
)

non_sigmut_output_file <- file.path(
  mutation_output_dir,
  paste0(
    sample_name,
    "_non_sigmuts.csv"
  )
)

mutation_output_file <- file.path(
  mutation_output_dir,
  paste0(
    sample_name,
    "_mutations.csv"
  )
)

variants_file <- file.path(
  variants_output_dir,
  paste0(
    sample_name,
    "_variants.csv"
  )
)

variants_with_meta_file <- file.path(
  variants_output_dir,
  paste0(sample_name, "_variants_with_meta.csv")
)


## ----process_signature_mutations, include = FALSE-----------------------------
# Read signature data
sigmut_df <- read.csv(mutation_sheet, header = TRUE) %>%

  # preprocess signature data
  # TODO Is the first select necessary? The sample sheet seems to have no
  # such column
  dplyr::select(-matches("source")) %>%
  dplyr::na_if("") %>%
  tidyr::pivot_longer(everything(), values_drop_na = TRUE) %>%
  dplyr::select(variant = name, mutation = value)

vep_output_df <- read.csv(params$vep_file, sep = ",", header = TRUE) %>%
  dplyr::na_if("-")

# group variants with the same muation and concat their names to preserve the
# information
sigmuts_deduped <- sigmut_df %>%
  group_by(mutation) %>%
  summarise(variant = paste(variant, collapse = var_sep))

# remove mutation type info / extract info on base changes
sigmuts_deduped_no_gene <- sigmuts_deduped %>%
  mutate(mutation = str_extract(mutation, "(?<=:)[[:alnum:]]+"))

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
# get the SNV frequency values and read depth information for the mutations from
# the LoFreq output
lofreq_info <- parse_snv_csv(params$snv_file)

complete_df <- dplyr::left_join(
  lofreq_info,
  variant_protein_mut,
  by = "gene_mut", copy = TRUE
  ) %>%
  mutate(gene_mut_collapsed = paste(genes, gene_mut, sep = ":"))

complete_dep_filtered_df <- complete_df %>%
 filter(as.numeric(dep) > as.numeric(params$mutation_depth_threshold))

# filter for mutations which are signature mutations
match_df <- complete_dep_filtered_df %>%
  filter(!is.na(variant))

# filter for everything that is not a signature mutation
nomatch_df <- complete_dep_filtered_df %>%
  filter(is.na(variant))

cat("Writing signature mutation file to ", sigmut_output_file, "...\n")
write.csv(
  match_df,
  sigmut_output_file
)

cat("Writing non signature mutation file to ", non_sigmut_output_file, "...\n")
write.csv(
  nomatch_df,
  non_sigmut_output_file
)

# Tables are displayed here in report

## ----echo = FALSE-------------------------------------------------------------
# get  NT mutations only, input for the signature matrix
mutations_vec <- match_df$gene_mut_collapsed

# only execute the deconvolution when at least one signature mutation was found
execute_deconvolution <- length(mutations_vec) > 0


if (execute_deconvolution) {
  ## ----creating_signature_matrix, include = FALSE-----------------------------
  # create an empty data frame add a column for the Wildtype
  # Wildtype in this case means the reference version of SARS-Cov-2
  # for the deconvolution to work we need the "wild type" frequencies too.
  # The matrix from above got mirrored, wild type mutations are simulated the
  # following: e.g. T210I (mutation) -> T210T ("wild type")
  msig_simple <- create_sig_matrix(mutations_vec, mutation_sheet)

  # deduplicate matrix

  variant_names <- colnames(msig_simple)

  is_dupe <- duplicated(msig_simple, MARGIN = 2)
  dupe_variants <- variant_names[is_dupe]

  # coerce back to dataframe for easier processing
  msig_simple_df <- as.data.frame(msig_simple)

  # find out of which variant a dupe variant is a dupe of, generate groups
  # of variants which are duplicates of each other
  dupe_group_list <- list()

  for (dupe_var in dupe_variants) {
    if (!dupe_var %in% unique(unlist(dupe_group_list))) {
      dupe_var_col <- msig_simple_df[[dupe_var]]

      dupe_group_logi <- apply(
        msig_simple,
        2,
        function(col, dupe_col) {
          all(col == dupe_col)
        }, dupe_var_col
      )

      dupe_group_vec <- variant_names[dupe_group_logi]

      dupe_group_list[[dupe_group_vec[1]]] <- dupe_group_vec
    }
  }

  # concat dupe groups to form a new composite name for the now unique col
  dupe_group_names <- lapply(dupe_group_list, paste, collapse = var_sep) %>%
    unlist()

  # juggle names to get a named vector with names and values flipped
  # needed by dplyr::rename()
  old_names <- names(dupe_group_names)
  new_names <- dupe_group_names

  dupe_group_names <- old_names %>%
    set_names(new_names)

  # generate deduped signature matrix
  # is a col was duplicated this contains only the first col of each dupe group
  msig_deduped_df <- msig_simple_df[, !is_dupe] %>%
    rename(!! dupe_group_names)

  ## ----calculate_sigmat_weigths, include = FALSE------------------------------
  deconv_lineages <- colnames(msig_deduped_df)

  # create list of proportion values that will be used as weigths
  sigmut_proportion_weights <- list()
  for (lineage in deconv_lineages) {
    if (lineage == "Others") {
      # !! 17/02/2022 It's not yet tested how robust this behaves when one would
      # mindlessly clutter the mutationsheet
      # with lineages that are very unlikely to detect or not detected

      # n all detected mutations / n all known mutations
      # TODO Find out why
      value <- nrow(msig_deduped_df) / nrow(sigmuts_deduped)
    } else if (grepl(var_sep, lineage)) {
      # as we can not weight the variants separately by their n detected sig
      # muts / n known sig muts (for this group), we use n detected sigmuts
      # (that is still accurate) / group average n known signature mutations
      group <- unlist(str_split(lineage, var_sep))
      avrg <- sum(sigmut_df$variant %in% group) / length(group)
      value <- sum(msig_deduped_df[lineage]) / avrg
    } else {
      # n lineage signature mutations detected in sample /
      # n known lineage signature mutations (provided in the mutation
      # sheet)
      value <- sum(msig_deduped_df[lineage]) /
        sum(sigmut_df$variant == lineage)
    }
    sigmut_proportion_weights[lineage] <- value
  }

  # apply weights to signature matrix
  msig_deduped_df_weighted <- msig_deduped_df %>%
    mutate(across(
      everything(), ~ .x / sigmut_proportion_weights[[cur_column()]]
    )) %>%
    replace(is.na(.), 0)


  ## ----simulating_WT_mutations, include = FALSE-------------------------------
  # construct additional WT mutations that are not weighted
  # for the deconvolution to work we need the "wild type" frequencies too. The
  # matrix from above got (elementwise) inverted, wild type mutations are
  # simulated the following: e.g. T210I (mutation) -> T210T ("wild type")
  #
  # for each mutation, generate a dummy mutation that results in no change

  # get bulk frequency values, will be input for the deconvolution function
  bulk_freq_vec <- as.numeric(match_df$freq)

  # generate frequency values for the dummy mutations. As they represent
  # all the variants without the respective mutation, they are the remainder
  # to one of the original mutations frequencies
  others_freq_vec <- 1 - bulk_freq_vec

  # generate vector with all mutation frequencies: dummy and real
  bulk_all <- c(others_freq_vec, bulk_freq_vec)

  # make matrix with Others mutations and inverse the values and wild type
  # freqs
  msig_inverse <- msig_deduped_df_weighted %>%
    mutate(across(everything(), ~ as.numeric(!as.logical(.x)))) %>%

    # apply weights right away
    mutate(across(
      everything(),
      ~ .x / as.numeric(sigmut_proportion_weights["Others"])
    ))


  # generate combined signature matrix for variants, dummy and real
  msig_all <- rbind(msig_inverse, msig_deduped_df_weighted) %>%
    as.matrix()

  # central deconvolution step
  variant_abundance <- deconv(bulk_all, msig_all)


  ## ----plot, echo = FALSE-----------------------------------------------------
  # work in progress...only to show how it theoretically can look like in the
  # report
  variant_abundance_df <- data.frame(
    variants = deconv_lineages,
    abundance = variant_abundance
  )

  # Handling of ambiguous cases and grouped variants

  # case 1: add dropped variants again with value 0 in case all of the other
  # variants add up to 1
  if (all.equal(sum(variant_abundance_df$abundance), 1)) {
    variant_abundance_df <- variant_abundance_df %>%
      bind_rows(data.frame(
        variants = dupe_variants,
        abundance = rep(0, length(dupe_variants))
        ))
  }

  # case 2: in case "others" == 0, both variants can be split up again and being
  # given the value 0 OR
  # case 3: in case multiple vars can really not be distinguished from each
  # other they will be distributed normaly

  while (any(str_detect(variant_abundance_df$variant, var_sep))) {
    grouped_rows <- which(str_detect(variant_abundance_df$variant, var_sep))
    # fixme: this loop might be unneccessary, since only the first row should
    # been picked, everything else will be handled by the while loop
    for (row in grouped_rows) {
      if (variant_abundance_df[row, "abundance"] == 0) {
        grouped_variants <- unlist(
          str_split(variant_abundance_df[row, "variant"], var_sep)
        )
        for (variant in grouped_variants) {
          # add new rows, one for each variant
          variant_abundance_df <- rbind(variant_abundance_df, c(variant, 0))
        }
      } else if (variant_abundance_df[row, "abundance"] != 0) {
        grouped_variants <- unlist(
          str_split(variant_abundance_df[row, "variant"], var_sep)
        )
        # normal distribution, devide deconv value by number of grouped variants
        distributed_freq_value <-
          as.numeric(as.numeric(variant_abundance_df[row, "abundance"]) /
            length(grouped_variants))
        for (variant in grouped_variants) {
          # add new rows, one for each variant
          variant_abundance_df <- rbind(
            variant_abundance_df,
              c(variant, distributed_freq_value)
            )
        }
      }
      # drop grouped row
      variant_abundance_df <- variant_abundance_df[-row, ]
    }
  }

  variant_abundance_df <- transform(
    variant_abundance_df, abundance = as.numeric(abundance)
  )


  cat("Writing variant abundance file to ", variants_file, "...\n")
  write.csv(variant_abundance_df, variants_file)

  # plot comes here in report
} else {
  # write dummy variants file
  # TODO: do this as a proper emty table with the correct col names
  file.create(variants_file)
}

# TODO: check if the else of the above if is handled correctly

## ----csv_output_variant_plot, include = F-------------------------------------
# prepare processed variant values to output them as a csv which will be used
# for the plots in index.rmd those outputs are not offically declared as outputs
# which can lead to issues - that part should be handled by a seperate
# file (and maybe rule)
output_variant_plot <- data.frame(
  samplename = character(),
  dates = character(),
  location_name = character(),
  coordinates_lat = character(),
  coordinates_long = character()
)
if (!execute_deconvolution) {
  # if no signatur mutation found write empty output file
  # TODO: sombody should check whether this empty file with header is enough,
  # or a more sensible default is required
  write.csv(output_variant_plot, variants_with_meta_file,
    na = "NA", row.names = FALSE, quote = FALSE
  )
} else {
  # get all possible variants
  all_variants <- colnames(msig_simple)
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

  # write mutation frequency values to variant_abundance_df
  for (i in all_variants) {
    if (i %in% variant_abundance_df$variant) {
      # check if variant already has a column
      if (i %in% colnames(output_variant_plot)) {
        output_variant_plot[sample_row, ][i] <- variant_abundance_df$abundance[
          variant_abundance_df$variant == i
        ]
        output_variant_plot <- output_variant_plot %>%
          mutate(others = 1 - rowSums(
            across(all_of(all_variants)), na.rm = TRUE
            ))
      }
    }
  }

  # 3. write to output file
  write.csv(output_variant_plot, variants_with_meta_file,
    na = "NA", row.names = FALSE, quote = FALSE
  )
}

if (execute_deconvolution) {
  ## ----csv_output_mutation_plot, include = FALSE------------------------------
  # prepare processed mutation values to output them as a csv which will be used
  # for the plots in index.rmd those outputs are not officially declared as
  # outputs which can lead to issues - that part should be handled by a seperate
  # file (and maybe rule)
  # get all possible mutations

  # one aa mutation can have different codon mutations reported with
  # different freqs- for the summary table they have to be summed up
  # (process see line 1872 of documentation)
  complete_df <- complete_df %>%
    group_by(across(c(-freq, -gene_mut, -gene_mut_collapsed, AA_mut))) %>%
    summarise(
      freq = sum(as.numeric(freq)),
      gene_mut = paste(gene_mut, collapse = var_sep)
    ) %>%
    rowwise() %>%
    mutate(AA_mut = replace(AA_mut, is.na(AA_mut), "\\:\\")) %>%
    # 211006 this exclusion is necessary because this mutation has a wrong entry
    # in VEP which gives two AA_muts instead of probably 1 deletion
    filter(!(gene_mut %in% "G13477A")) %>%
    ungroup()

  # report the gene, translated_AA_mut and NT mut accordingly
  # easier to spot translation inconsitentcies that way
  all_mutations <- paste(complete_df$AA_mut[!is.na(complete_df$AA_mut)],
    complete_df$gene_mut,
    sep = "::"
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

  # write mutation frequency values to variant_abundance_df
  for (i in all_mutations) {
    i_nt <- str_split(i, "::")[[1]][2]
    if (i_nt %in% complete_df$gene_mut) { # split gene name to match with AA mut
      # check if variant already has a column
      if (i %in% colnames(output_mutation_frame)) {
        output_mutation_frame[, i] <- complete_df$freq[which(
          complete_df$gene_mut == i_nt
        )]
      }
    }
  }


  colnames(output_mutation_frame) <- as.character(colnames(
    output_mutation_frame
  ))

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
  cat("Writing mutation file to ", mutation_output_file, "...\n")
  write.csv(output_mutation_frame, mutation_output_file,
    row.names = FALSE, quote = FALSE
  )
} else {

  cat("Writing dummy mutation file to ", mutation_output_file, "...\n")

  writeLines(
    "Deconvolution not run, this is a dummy file.",
    mutation_output_file
  )
}
