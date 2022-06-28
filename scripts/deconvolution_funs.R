library("MASS")
library("stringr")
library("dplyr")
library("data.table")

parse_snv_csv <- function(snvfile) { # allele frequency from v-pipe vcf
  #' input: csv file derived from vpipe vcf when using LoFreq,
  #' parsing snv-csv file for coverage, frequency and genomic mutation information

  snvtable <- fread(snvfile, sep = ",", header = TRUE)

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
  # TODO Currently unused, see note on deletions in get_protein_mut
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
  vepfile_df <- fread(vepfile, sep = ",", header = TRUE)
  # parsing snv and protmut location


  locations <- vepfile_df %>%

    rename(
      prot_pos  = Protein_position,
      conseq    = Consequence,
      gene_name = SYMBOL) %>%

    # get general info on mutation and its position
    mutate(
      mut_nucs   = str_extract(Uploaded_variation, "[A-Z*]/[A-Z*]$"),
      mut_chrom  = str_extract(Location, "^[^:]+"),
      mut_start  = str_extract(Location, "(?<=:)[0-9]+[^-]"),
      mut_end    = str_extract(Location, "(?<=-)[0-9]+")) %>%

    # generate unique string describing mutation
    mutate(
      mut_ref = str_extract(mut_nucs, "^[A-Z*]"),
      mut_var = str_extract(mut_nucs,  "[A-Z*]$"),

      mut_str = paste0(mut_ref, mut_start, mut_var)
    ) %>%

    # get infos on mutation consequences for protein
    mutate(
      # TODO This may behave unexpectedly; in the case of no change aa_ref
      # and aa_var are the same, aa_var will not be empty / NA
      aa_ref = str_extract(Amino_acids, "^[A-Z*]"),
      aa_var = str_extract(Amino_acids,  "[A-Z*]$"),

      aa_str = paste0(aa_ref, prot_pos, aa_var)
    ) %>%

    # delete all rows with no protein position value
    # TODO Why do we delete these mutations? This just means they do not
    # affect the proteins, the mutations may still be informative...
    filter(! str_detect(prot_pos, "^-")) %>%
    # TODO What is the relevance of this information?
    # specific B117 mutations: 21990-21993, 21764-21770, maybe also 3675-3677,
    # 69-70 - all there

    # remove unneeded cols
    dplyr::select(-matches(names(vepfile_df), ignore.case = FALSE)) %>%

    # This also removes variants that essentially look the same, but differ in
    # their VEP determined Gene
    # e.g. Gene:ENSSASG00005000002 <> Gene:ENSSASG00005000003
    distinct()

    # TODO Previously deletions were handled by generating pseudo rows for each
    # one. I currently do not have access to a file with deletions, and so
    # cannot understand how to handle this properly. However this solution here
    # is still better, because due to a bug the previous solution could almost
    # never run. Now there is at least a warning...
    if (any(!is.na(locations$mut_end))) {
      warning("Found a deletion. Help.")
    }

  return(locations)
}


create_sig_matrix <- function(mutations_vector, mutation_sheet_file) {
  #' for making the signature matrix based on the signature mutations found in
  #' the sample (given as input as a vector)for it self
  #' returns simple signature matrix as data.frame without frequency values

  # read in provided mutation sheet
  mutations_df <- fread(mutation_sheet_file) %>%

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


dedupe_sigmut_mat <- function(sigmut_mat, var_sep = "_") {
  # input:
  #   * sigmut_mat:
  #     A matrix with rows corresponing to singnature mutations and cols
  #     corresponding to variants, each dimension named. Entries are 1 if the
  #     mutation of the row is a signature mutation of the variant in the col,
  #     and 0 otherwise.
  #   * var_sep:
  #     A string separating the names of variables with identical columns.
  # output:
  #   * The input matrix with identical columns and their colnames merged.
  #   * A list of character vectors, giving names of variants which have equal
  #     cols.
  variant_names <- colnames(sigmut_mat)

  is_dupe <- duplicated(sigmut_mat, MARGIN = 2)
  dupe_variants <- variant_names[is_dupe]

  # coerce back to dataframe for easier processing
  sigmut_mat_df <- as.data.frame(sigmut_mat)

  # find out of which variant a dupe variant is a dupe of, generate groups
  # of variants which are duplicates of each other
  dupe_group_list <- list()

  for (dupe_var in dupe_variants) {
    if (!dupe_var %in% unique(unlist(dupe_group_list))) {
      dupe_var_col <- sigmut_mat_df[[dupe_var]]

      dupe_group_logi <- apply(
        sigmut_mat,
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
  msig_deduped_df <- sigmut_mat_df[, !is_dupe] %>%
    rename(!!dupe_group_names) %>%
    replace(is.na(.), 0)

  return(msig_deduped_df)
}