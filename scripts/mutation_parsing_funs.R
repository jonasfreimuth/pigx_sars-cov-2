# TODO: ensure necessary packages are available

parse_snv_csv <- function(snvfile, ...) { # allele frequency from v-pipe vcf
  #' input: csv file derived from vpipe vcf when using LoFreq,
  #' parsing snv-csv file for coverage, frequency and genomic mutation information

  snvtable <- read.table(snvfile, sep = ",", header = TRUE)
  freq <- snvtable$AF
  cov <- snvtable$DP
  Ref <- snvtable$Ref
  Pos <- snvtable$Pos
  Var <- snvtable$Var

  # concat position value and nucleotides to nucleotide-mutation-notation
  snvinfo.df <- data.frame(
    Ref = Ref,
    Pos = Pos,
    Var = Var,
    gene_mut = paste0(Ref, Pos, Var),
    stringsAsFactors = FALSE
  )

  # concat nucleotide-mutation-notation with mutation frequencies and coverage
  # seperate posistion column will be used for joining data frames
  snv.info <- cbind(
    gene_pos = snvinfo.df[, "Pos"],
    gene_mut = snvinfo.df[, "gene_mut"],
    freq, cov
  )

  return(snv.info)
}

detectable_deletions <- function(x, colnames) {
  #' deletions can span a range of position values indicated by a dash in the
  #' prot_mut_loc column. In those cases those ranges have to be split up in
  #' order to be able to detect single positions. This function should be used
  #' on a already filtered set of deletions derived from the vep-output. They are
  #' then split up and expanded with the number of rows related to the number
  #' of position. The function returns those extended lines as a dataframe which
  #' structure is matching the df structure which is returned by
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
  vepfile.df <- read.table(vepfile, sep = ",", header = TRUE)
  # parsing snv and protmut location

  # parsing gene mutation
  gene_mutation <- data.frame(
    gene_mut_loc = str_split_fixed(vepfile.df$Location, "[:-]+", n = 3),
    nucleotides = str_split_fixed(
      str_split_fixed(vepfile.df$Uploaded_variation,
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
    gene_mut_loc = str_split_fixed(vepfile.df$Location, "[:-]+", n = 3),
    gene_mut = gene_mutation$gene_mut,
    prot_mut_loc = vepfile.df$Protein_position,
    AAs = str_split_fixed(vepfile.df$Amino_acids, "/", 2),
    Conseq = vepfile.df$Consequence,
    genes = vepfile.df$SYMBOL
  )

  locations <- dplyr::na_if(locations, "")

  # delete all rows with no protein position value
  locations <- distinct(locations %>% filter(!grepl("^-", prot_mut_loc)))
  # specific B117 mutations: 21990-21993, 21764-21770, maybe also 3675-3677,
  # 69-70 - all there

  # deletions.df <- locations TODO: Check effect of this assignment
  deletions.df <- locations %>%
    filter(gene_mut_loc.2 != gene_mut_loc.3 &
      Conseq == "inframe_deletion" &
      str_detect(deletions.df$prot_mut_loc, "-"))

  colnames <- colnames(deletions.df)

  # if there is a deletion the snv would span a couple of positions, if there is
  # not such a spanning region there are no deletions
  # ! 06/05/2021 Vic - I think, I don't know how robust this is, but it will
  # work for the sig mutations we have so far
  if (nrow(locations) >= 1 && !(any(is.na(locations[, "gene_mut_loc.3"])))) {
    deletions <- dplyr::bind_rows(apply(deletions.df, 1, detectable_deletions,
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


dedupeMuts <- function(mut, sigmut.df, dedup.df) {
  #' this function is a different version of "dedupe"
  #' if a signature mutation is shared by multiple variants the mutation-tables
  #' would have seperate rows for the same mutation. This function is for
  #' concatenating all the variants that share that mutation. It requires the
  #' original df with all mutations, the deduplicated one (only one row per
  #' mutation) and the mutation for which this procedure should be applied
  #' This function is not yet completely refined and should be improved further
  # get variant group per mutation
  duped_muts <- grep(mut, sigmut.df$value)
  grouped <- sigmut.df$name[duped_muts]
  groupName <- paste(grouped, collapse = ",")

  dedup.df$name[dedup.df$value == mut] <- groupName
  return(dedup.df)
}