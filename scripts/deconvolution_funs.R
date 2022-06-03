library("stringr")
library("dplyr")

parse_snv_csv <- function(snvfile, ...) { # allele frequency from v-pipe vcf
  #' input: csv file derived from vpipe vcf when using LoFreq,
  #' parsing snv-csv file for coverage, frequency and genomic mutation information

  snvtable <- read.csv(snvfile, sep = ",", header = TRUE)

  ref <- snvtable$Ref
  pos <- snvtable$Pos
  var <- snvtable$Var

  snv_info_df <- data.frame(
    gene_pos = pos,
    gene_mut = paste0(ref, pos, var),
    freq     = snvtable$AF,
    dep      = snvtable$DP
  )

  return(snv_info_df)
}

detectable_deletions <- function(x, colnames) {
  #' deletions can span a range of position values indicated by a dash in the
  #' prot_mut_loc column. In those cases those ranges have to be split up in
  #' order to be able to detect single positions. This function should be used
  #' on a already filtered set of deletions derived from the vep-output. They
  #' are then split up and expanded with the number of rows related to the
  #' number of position. The function returns those extended lines as a
  #' dataframe which structure is matching the df structure which is returned by
  #' "get_protein_mut"

  # check that the mutations spans multiple nucleotides AND multiple Amino Acids
  if (x["gene_mut_loc.2"] != x["gene_mut_loc.3"] &
    str_detect(x["prot_mut_loc"], "-")) {

    # extract columns from input dataframe where the content of the rows won't
    # change, but where rows will be duplicated when merging back with the
    # extended deletion rows
    constant_before <- as.data.frame((cbind(
      x["gene_mut_loc.1"],
      x["gene_mut_loc.2"],
      x["gene_mut_loc.3"]
    )))

    constant_after <- as.data.frame((cbind(
      x["Conseq"],
      x["genes"]
    )))

    # split the position values denoting a range and make as many new rows as
    # positions spanned with the missing position values e.g: 1-3 (1 row) will
    # become 1,2,3 (3 rows)
    split_prot_pos <- str_split(x["prot_mut_loc"], "-")
    list_prot_pos <- list()
    for (i in 1:length(split_prot_pos)) {
      list_prot_pos[[i]] <- seq(
        as.numeric(split_prot_pos[[i]][1]),
        as.numeric(split_prot_pos[[i]][2])
      )
    }

    # TODO instead of qpcR bind_cols from dplyer
    # split groups of reference Amino acids (AAs.1) if necessary, concat the
    # extendet columns
    extendet <- as.data.frame(
      qpcR:::cbind.na(
        as.character(unlist(list_prot_pos)),
        as.character(unlist(str_split(x["AAs.1"], ""))),
        as.character(str_split(x["AAs.2"], ""))
      )
    )

    # the dash is the sign for an Amino Acid that was deleted
    extendet[is.na(extendet)] <- "-"

    # join the columns with the extended rows back with the rest of the original
    # dataframe
    full <- as.data.frame(cbind(constant_before, extendet, constant_after))
    colnames(full) <- colnames

    # match the extended rows with added positions values to the Amino acids.
    # If a group of reference Amino acids in AAs.1 is given it's split and
    # matched according to their order. The mutation in AAs.2 is only applied
    # to the last position every other position get's an dash in column AAs.2
    # e.g.: 1-3 is extended to separate rows, if there was ABC in AAs.1 and C in
    # AAs.2, the results would be A1-, B2-, C3C
    full$prot_mut_loc <- vapply(full$prot_mut_loc, paste,
      collapse = " ",
      character(1L)
    )
    full$AAs.1 <- vapply(full$AAs.1, paste, collapse = " ", character(1L))
    full$AAs.2 <- vapply(full$AAs.2, paste, collapse = " ", character(1L))

    return(full)
  }
}

get_protein_mut <- function(vepfile) {
  # parse together from vep output "Protein_position" and "Amino_acid"
  #' input: []_sarscov2_parsed.txt, parsed vcf from VEP CLI, comma separated
  #' this function is for parsing the information about mutation position in the
  #' amino acid sequence, the reference aa and the alternative mutation into the
  #' aa-mutation-notation which is later on comparable to the lists of signature
  #' mutations

  # reading in whole vep txt output
  # you should include in the vep script to parse out the #
  # in the beginning of the line or include that step here.
  # TODO: include check about correct VEP file input format
  vepfile_df <- read.csv(vepfile, sep = ",", header = TRUE)
  # parsing snv and protmut location

  # parsing gene mutation
  gene_mutation <- data.frame(
    gene_mut_loc = str_split_fixed(vepfile_df$Location, "[:-]+", n = 3),
    nucleotides = str_split_fixed(
      str_split_fixed(vepfile_df$Uploaded_variation,
        "[_-]+",
        n = 4
      )[, 4], "/",
      n = 2
    )
  )

  gene_mutation$gene_mut <- paste0(
    gene_mutation$nucleotides.1,
    gene_mutation$gene_mut_loc.2,
    gene_mutation$nucleotides.2
  )

  # parsing snv and protmut location
  locations <- data.frame(
    gene_mut_loc = str_split_fixed(vepfile_df$Location, "[:-]+", n = 3),
    gene_mut = gene_mutation$gene_mut,
    prot_mut_loc = vepfile_df$Protein_position,
    AAs = str_split_fixed(vepfile_df$Amino_acids, "/", 2),
    Conseq = vepfile_df$Consequence,
    genes = vepfile_df$SYMBOL
  )

  locations <- dplyr::na_if(locations, "")

  # delete all rows with no protein position value
  locations <- distinct(locations %>% filter(!grepl("^-", prot_mut_loc)))
  # specific B117 mutations: 21990-21993, 21764-21770, maybe also 3675-3677,
  # 69-70 - all there

  deletions_df <- locations %>%
    filter(gene_mut_loc.2 != gene_mut_loc.3 &
      Conseq == "inframe_deletion" &
      str_detect(prot_mut_loc, "-"))

  colnames <- colnames(deletions_df)

  # if there is a deletion the snv would span a couple of positions, if there is
  # not such a spanning region there are no deletions
  # ! 06/05/2021 Vic - I think, I don't know how robust this is, but it will
  # work for the sig mutations we have so far
  if (nrow(locations) >= 1 && !(any(is.na(locations[, "gene_mut_loc.3"])))) {
    deletions <- dplyr::bind_rows(apply(deletions_df, 1, detectable_deletions,
      colnames = colnames
    ))
    locations <- dplyr::bind_rows(locations, deletions)
  }

  # substitute "nothing" at alternative-aa-column (for deletions) with "-"
  locations$AAs.2[is.na(locations$AAs.2)] <- "-"

  # clean - characters
  locations$AA_mut <- paste(locations$AAs.1,
    locations$prot_mut_loc,
    locations$AAs.2,
    sep = ""
  )

  # adding gene information
  locations$AA_mut <- paste(locations$genes,
    locations$AA_mut,
    sep = ":"
  )

  return(locations)
}


create_sig_matrix <- function(mutations_vector, mutation_sheet_file) {
  #' for making the signature matrix based on the signature mutations found in
  #' the sample (given as input as a vector)for it self
  #' returns simple signature matrix as data.frame without frequency values

  # read in provided mutation sheet
  mutations_df <- read.csv(mutation_sheet_file) %>%

    # remove source col if there. TODO What problems would this col cause?
    dplyr::select(-matches("source"))

  # making a matrix with the signature mutations found in the sample
  # make binary matrix matching the mutations to the mutation-list per variant
  # to see how many characterising mutations where found by variant

  sig_mat <- lapply(colnames(mutations_df), function(variant) {
    return(as.numeric(mutations_vector %in% mutations_df[[variant]]))
  }) %>%

    dplyr::bind_cols() %>%

    magrittr::set_names(names(mutations_df)) %>%

    # add "Others" col of all "0"s to indicate possible other variants not
    # possessing any of the mutations
    # TODO ensure this is how this column is supposed to work
    mutate(Others = rep(0, length(mutations_vector))) %>%

    # coerce to matrix
    as.matrix() %>%

    magrittr::set_rownames(mutations_vector)


  return(sig_mat)
}

simulate_others <- function(
  mutations_vector, bulk_freq_vektor,
  simple_sigmat_dataframe, coverage_vektor, others_weight
  ) {
  #' for the deconvolution to work we need the "wild type" frequencies too. The
  #' matrix from above got mirrored, wild type mutations are simulated the
  #' following: e.g. T210I (mutation) -> T210T ("wild type")

  # 1. make "Others mutations"
  muts_others <- lapply(mutations_vector, function(x) {
    str_replace(
      x, regex(".$"),
      str_sub(str_split(x, ":")[[1]][2], 1, 1)
    )
  })
  muts_others_df <- data.frame(muts = unlist(muts_others))
  # 2. make frequency values, subtract the observed freqs for the real mutations
  # from 1
  bulk_others <- lapply(bulk_freq_vektor, function(x) {
    1 - x
  })

# 3. make matrix with Others mutations and inverse the values and wild type
# freqs
msig_inverse <- bind_cols(
  muts_others_df,
  as.data.frame(+(!simple_sigmat_dataframe))
)

  # 4. apply Others weight
  # fixme: it could be this can be implemented in the step above already
  msig_inverse[msig_inverse == 1] <- 1 / others_weight

  # fixme: not sure if this really is a nice way to concat those things...
  muts_all <- c(muts_others, mutations_vector)
  muts_all_df <- data.frame(muts = unlist(muts_all))

  bulk_all <- c(bulk_others, bulk_freq_vektor)
  bulk_all_df <- data.frame(freq = unlist(bulk_all))

  coverage_all <- c(coverage_vektor, coverage_vektor)
  coverage_all_df <- data.frame(cov = unlist(coverage_all))

  msig_all <- rbind(
    msig_inverse[, -which(names(msig_inverse) %in% "muts")],
    simple_sigmat_dataframe
  )

  # 4. concat the data frames
  # without bulk freq for building the signature matrix
  msig_stable <- bind_cols(muts_all_df, msig_all)

  # with bulk freq for export and overview
  msig_stable_complete <- bind_cols(
    muts_all_df, msig_all, bulk_all_df,
    coverage_all_df
  )

  return(list(msig_stable, bulk_all, msig_stable_complete))
}

# When multiple columns look like the same, the deconvolution will not work,
# because the function can't distinguish between those columns. The workaround
# for now is to identify those equal columns and merge them into one, returning
# also a vector with the information about which of the columns were merged.
# deduplicate dataframe
dedupe_df <- function(msig_stable) {
  # transpose and add mutations as first column
  msig_stable_transposed <- as.data.frame(cbind(
    variants = colnames(msig_stable),
    t(msig_stable)
  ))

  # mark duplicated columns, forward and backwards to get ALL the duplicates,
  # otherwise the first one would missing
  dupes_variants <- duplicated(
    msig_stable_transposed[
      , -which(names(msig_stable_transposed) %in% "variants")
      ],
    fromLast = TRUE
  )

  msig_dedupe_transposed <- msig_stable_transposed[!dupes_variants, ]

  return(list(msig_stable_transposed, msig_dedupe_transposed))
}

dedupe_variants <- function(variant, variants_df, depup_variants_df) {
  # get variant group per mutation pattern
  duped_variants <- c()
  row_number_variant <- which(grepl(variant, variants_df$variants))
  for (row in 1:nrow(variants_df)) {
    # TODO: what are those magic numbers?
    if (all(variants_df[row_number_variant, -1] == variants_df[row, -1])) {
      duped_variants <- c(duped_variants, variants_df[row, "variants"])
    }
  }
  group_name_variants <- paste(duped_variants, collapse = ",")
  for (row in depup_variants_df$variants) {
    if (grepl(row, group_name_variants)) {

      # if variants are getting pooled with Others they are just Others and
      # nothing else
      if (str_detect(group_name_variants, "Others")) {
        rownames(depup_variants_df)[
          rownames(depup_variants_df) == row
          ] <- "Others"
        variants_to_drop <- duped_variants[!grepl("Others", duped_variants)]
      } else {
        rownames(depup_variants_df)[
          rownames(depup_variants_df) == row
          ] <- group_name_variants
        variants_to_drop <- NA
        # TODO you can stop after this ( I think)
      }
    }
  }
  # clean the vector to know which variants has to be add with value 0 after
  # deconvolution
  variants_to_drop <- unique(variants_to_drop)[!is.na(variants_to_drop)]
  return(list(depup_variants_df, variants_to_drop))
}

deconv <- function(bulk, sig) {
  #' This function performs the deconvolution using a signature matrix for the
  #' mutations found in the sample and bulk frequency values derived by the SNV
  #' caller
  #' it was build by Altuna

  rlm_model <- suppressWarnings(MASS::rlm(sig, bulk, maxit = 100, method = "M"))


  rlm_coefficients <- rlm_model$coefficients

  rlm_coefficients <- ifelse(rlm_coefficients < 0, 0, rlm_coefficients)

  sum_of_cof <- sum(rlm_coefficients)

  # normalize so coefficients add to 1
  rlm_coefficients <- rlm_coefficients / sum_of_cof

  as.vector(rlm_coefficients)
}

deconv_debug <- function(bulk, sig) {
  #' This function performs the deconvolution using a signature matrix for the
  #' mutations found in the sample and bulk frequency values derived by the SNV
  #' caller
  #' it was build by Altuna

  rlm_model <- suppressWarnings(MASS::rlm(sig, bulk, maxit = 100, method = "M"))

  rlm_coefficients <- rlm_model$coefficients

  rlm_coefficients <- ifelse(rlm_coefficients < 0, 0, rlm_coefficients)

  sum_of_cof <- sum(rlm_coefficients)

  # normalize so coefficients add to 1
  rlm_coefficients <- rlm_coefficients / sum_of_cof

  return(list(as.vector(rlm_coefficients), rlm_model$fitted.values))
}
