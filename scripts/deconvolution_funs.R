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

simulate_others <- function(mutations_vector,
                            bulk_freq_vector,
                            simple_sigmat_dataframe,
                            coverage_vector,
                            others_weight) {
  #' for the deconvolution to work we need the "wild type" frequencies too. The
  #' matrix from above got (elementwise) inverted, wild type mutations are
  #' simulated the following: e.g. T210I (mutation) -> T210T ("wild type")
  #'
  #' for each mutation, generate a dummy mutation that results in no change

  # generate frequency values for the dummy mutations. As they represent
  # all the variants without the respective mutation, they are the remainder
  # to one of the original mutations frequencies
  inv_freq_vec <- 1 - bulk_freq_vector

  # make matrix with Others mutations and inverse the values and wild type
  # freqs
  msig_inverse <- simple_sigmat_dataframe %>%
    mutate(across(everything(), ~ as.numeric(!as.logical(.x)))) %>%

    # apply weights right away
    mutate(across(everything(), function(x) {
      x[x == 1] <- 1 / others_weight

      return(x)
      }
    ))

  bulk_all <- c(inv_freq_vec, bulk_freq_vector)

  msig_all <- rbind(
    msig_inverse[, -which(names(msig_inverse) %in% "muts")],
    simple_sigmat_dataframe
  )

  return(list(msig_all, bulk_all))
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
