#' @title Files2Phyloseq
#'
#' @description
#' This function is used to convert DADA2, QIIME2, Kraken2, Bracken, MetaPhlAn3 and HUMAnN3 output files into a phyloseq object.
#' It accepts multiple file formats: qiime artifact - "qza", comma-separated - "csv", tab-separated - "tsv",
#' Newick-format - "tree" , Excel - "xls, xlsx, xlsm" and text - "txt. If file format is "txt" please specify field separator.
#' Filtration of specific samples based on metadata variables is also possible.
#'
#' @param BIOM.file (Required if not providing OTU and TAXA files) A character string of the PATH to the BIOM file
#' @param OTU.file (Required if not providing BIOM file) A character string of the PATH to the ASV/OTU/Species abundance table file
#' (accepted formats: qza, csv, tsv and txt if "sep" is specified).
#' @param TAXA.file (Required if not providing BIOM file) A character string of the PATH to the taxonomic classification file
#' (accepted formats: qza, csv, tsv and txt if "sep" is specified).
#' @param META.file (Optional) A character string of the PATH to the metadata file
#' (accepted formats: xls, xlsx, xlsm, csv, tsv and txt if "sep" is specified).
#' The first column must be named "SampleID" and must contain the sample names.
#' @param TREE.file (Optional) A character string of the PATH to the phylogenetic tree file
#' (accepted formats: qza and tree).
#' @param unclass.file (Optional) A character string of the PATH to the unclassified file.
#' When using Bracken the "Unclassified" taxon is removed, to add this information extract it from the Kraken report and create a txt file
#' with the the sample names in the first column and the corresponding "Unclassified" values in the second column.
#' Use the Extract_unclassified_reads.sh function to do this.
#' (accepted formats: tsv)
#' @param mphlan.file (Optional) A character string of the PATH to the MetaPhlAn3 file.
#' (accepted formats: tsv)
#' @param humann.file (Optional) A character string of the PATH to the HUMAnN3 file.
#' (accepted formats: tsv)
#' @param sep (Optional) 	A character string of the field separator. Values on each line of the file are separated by this character.
#' Only required if the OTU, TAXA or META files are in "txt" format. DEFAULT: "\t"
#' @param str2rm (Optional) A vector of character strings defining the patterns to be removed from the names of the samples. DEFAULT: ""
#' @param factor.lev (Optional) A list containing the names of the levels found in each categorical variable.
#' This will establish in which order level names will be displayed in the plots. DEFAULT: NULL
#' @param filt.var (Optional) A character string of the metadata variable to use for the filtration of samples. DEFAULT: NULL
#' @param lev.kept (Optional) A vector containing the names of the levels from the 'filt.var' that will be retained. DEFAULT: NULL
#' @return
#' A phyloseq object in which the variables defined as factors have their levels organised in the specified order.
#'
#' @examples
#' To create a phyloseq object from the DADA2 output files:
#' Files2Phyloseq("~/DADA2/seqtab.csv", "~/DADA2/taxaId.csv", "~/DADA2/Metadata.xlsm", factor.lev = list(Type = c("Water", "Sediment"), Source = c("Market", "Wild" )))
#'
#' To create a phyloseq object from the Bracken output files:
#' Files2Phyloseq("~/DADA2/seqtab.csv", "~/DADA2/taxaId.csv", "~/DADA2/Metadata.xlsm", factor.lev = list(Type = c("Water", "Sediment"), Source = c("Market", "Wild" )))
#' @export
Files2Phyloseq <- function(BIOM.file = NULL, OTU.file = NULL, TAXA.file = NULL, META.file = NULL, TREE.file = NULL,
                           unclass.file = NULL, mphlan.file = NULL, humann.file = NULL,
                           sep = "\t", str2rm = "", factor.lev = NULL, filt.var = NULL, lev.kept = NULL){

  ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

  # Create df with unclassified data -- ONLY for Bracken output
  if (!is.null(unclass.file)){
    unclass.data <- read.table(unclass.file)
    unclass.data %<>% tidyr::pivot_wider(names_from = V1, values_from = V2)
    unclass.data %<>% as.data.frame(.)
    unclass.data[1,] %<>% as.numeric(.)
    rownames(unclass.data) <- 0
    for (i in str2rm){
      colnames(unclass.data) <- sub(i, "", colnames(unclass.data))
    }
  }

  # Import BIOM or OTU/ASV and TAXA tables
  if (is.null(mphlan.file) & is.null(humann.file)){
    if (is.null(BIOM.file)){
      # Import and edit OTU/ASV table
      ifelse(tools::file_ext(OTU.file) == "csv", {OTU <- read.csv(OTU.file, sep = ",")},
             ifelse(tools::file_ext(OTU.file) == "tsv", {OTU <- read.csv(OTU.file, sep = "\t")},
                    ifelse(tools::file_ext(OTU.file) == "txt", {OTU <- read.csv(OTU.file, sep = sep)},
                           ifelse(tools::file_ext(OTU.file) == "qza", {OTU <- qiime2R::read_qza(OTU.file)$data},
                                  cat("Incorrect OTU.file format!!\n")))))
      ifelse(!tools::file_ext(OTU.file) == "qza", rownames(OTU) <- OTU[,1], "")
      OTU %<>% .[, names(.) != "X"]
      # Clean sample names
      for (i in str2rm){
        colnames(OTU) %<>% sub(i, "", .)
      }
      if (!is.null(unclass.file)){
        un.OTU <- full_join(unclass.data, OTU)
        rownames(un.OTU) <- c(0, rownames(OTU))
        OTU <- un.OTU
      }
      otumat <- as.matrix(OTU)
      OTU = otu_table(otumat, taxa_are_rows = TRUE)

      # Import and edit TAXA table
      ifelse(tools::file_ext(TAXA.file) == "csv", {TAXA  <- read.csv(TAXA.file, sep = ",")},
             ifelse(tools::file_ext(TAXA.file) == "tsv", {TAXA  <- read.csv(TAXA.file, sep = "\t")},
                    ifelse(tools::file_ext(TAXA.file) == "txt", {TAXA  <- read.csv(TAXA.file, sep = sep)},
                           ifelse(tools::file_ext(TAXA.file) == "qza", {TAXA <- qiime2R::read_qza(TAXA.file)$data},
                                  cat("Incorrect TAXA.file format!!\n")))))
      ifelse(!tools::file_ext(TAXA.file) == "qza", rownames(TAXA) <- TAXA[,1],
             ifelse(tools::file_ext(TAXA.file) == "qza", {TAXA <- qiime2R::parse_taxonomy(TAXA)}, ""))
      TAXA <- TAXA[, names(TAXA) != "X"]
      colnames(TAXA) <- ranks
      if (!is.null(unclass.file)){
        # Add Unclassified taxon with 0 as rowname
        TAXA = rbind(matrix(rep("Unclassified", ncol(TAXA)), nrow=1, dimnames = list(c(0), colnames(TAXA))), TAXA)
      }
      taxmat <- as.matrix(TAXA)
      TAXA = tax_table(taxmat)
    } else {
      defaultW <- getOption("warn") # to hide BIOM warning message
      options(warn = -1)
      BIOM = import_biom(BIOM.file, parseFunction=parse_taxonomy_greengenes, parallel=TRUE, version=1.0)
      options(warn = defaultW)
      # Clean sample names
      for (i in str2rm){
        colnames(otu_table(BIOM)) %<>% sub(i, "", .)
      }
      colnames(BIOM@tax_table) <- ranks
      if (!is.null(unclass.file)){
        otu_mat <- as.matrix(rbind(unclass.data, otu_table(BIOM)))
        OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
        TAXA = tax_table(rbind((matrix(rep("Unclassified", ncol(tax_table(BIOM))), nrow=1, dimnames = list(c(0), ranks))),
                              tax_table(BIOM)))
      } else {
        OTU = otu_table(BIOM)
        TAXA = tax_table(BIOM)
      }
    }
  }

  # Import MetaPhlAn3 table
  if(!is.null(mphlan.file)){
    mphlan.tbl <- read.csv(mphlan.file, sep = "\t", strip.white = T, stringsAsFactors = F, skip = 1, row.names = 1)
    # Clean sample names
    for (i in c('_metaphlan3', str2rm)){
      colnames(mphlan.tbl) %<>% sub(i, "", .)
    }
    split.names <- strsplit(rownames(mphlan.tbl), split = "|", fixed = TRUE)
    # let's filter the rows that have rownames which don't reach the species rank, we keep UNKNOW
    n.ranks <- max(sapply(split.names, length)) # get number of ranks
    mphlan.tbl %<>% .[(lengths(split.names) == n.ranks | rownames(.) == "UNKNOWN"),-1]
    split.names %<>% .[lengths(.) == n.ranks]
    split.names %<>% append(list(rep("Unclassified", n.ranks)), .)
    rownames(mphlan.tbl)[rownames(mphlan.tbl) == "UNKNOWN"] <- "Unclassified"
    tax_mat <- matrix(NA,
                      ncol = n.ranks,
                      nrow = length(split.names))
    colnames(tax_mat) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain")[1:n.ranks]
    rownames(tax_mat) = rownames(mphlan.tbl)
    for (i in 1:nrow(tax_mat)){
      tax_mat[i, 1:length(split.names[[i]])] <- split.names[[i]]
    }
    tax_mat = gsub("[a-z]__", "", tax_mat)
    TAXA = tax_table(tax_mat)
    OTU = otu_table(mphlan.tbl, taxa_are_rows = TRUE)
  }

  # Import HUMAnN3 table
  if(!is.null(humann.file)){
    humann <- file2meco::humann2meco(abund_table = humann.file, db = "MetaCyc")
    humann <- file2meco::meco2phyloseq(humann)
    for (i in c('_Abundance-RPKs', '_Abundance', '_Coverage', str2rm)){
      colnames(humann@otu_table) %<>% sub(i, "", .)
    }
    OTU = otu_table(humann)
    TAXA = tax_table(humann)
  }

  # Import METADATA table
  if(!is.null(META.file)){
    ifelse(tools::file_ext(META.file) %in% c("xls", "xlsx", "xlsm"), {META  <- readxl::read_excel(META.file)},
           ifelse(tools::file_ext(META.file) == "csv", {META  <- read.csv(META.file, sep = ",")},
                  ifelse(tools::file_ext(META.file) == "tsv", {META  <- read.csv(META.file, sep = "\t")},
                         ifelse(tools::file_ext(META.file) == "txt", {META  <- read.csv(META.file, sep = sep)},
                                cat("Incorrect META.file format!!\n")))))
    META <- sample_data(META)
    names(META)[1] <- "SampleID"
    rownames(META) <- META$SampleID
  }

  # Import phylogenetic TREE
  if(!is.null(TREE.file)){
    ifelse(tools::file_ext(TREE.file) == "qza", {TREE<- qiime2R::read_qza(TREE.file)$data},
           ifelse(tools::file_ext(TREE.file) == "tree", {TREE <- phy_tree(read_tree(TREE.file))},
                  cat("Incorrect TREE.file format!!\n")))
  }

  # Create phyloseq object
  if(exists("TREE") & exists("META")){
    data_phylo <- phyloseq(OTU, TAXA, META, TREE)
  } else if(exists("META")){
    data_phylo <- phyloseq(OTU, TAXA, META)
  } else {
    data_phylo <- phyloseq(OTU, TAXA)
  }

  # Remove specified samples
  if(!is.null(filt.var) & !is.null(lev.kept)){
    sam.data <- data.frame(data_phylo@sam_data)
    samples.kept = rownames(sam.data[sam.data[,filt.var] %in% lev.kept,])
    data_phylo = prune_samples(samples.kept, data_phylo)
    data_phylo = prune_taxa(taxa_sums(data_phylo) != 0, data_phylo) # Remove taxa with 0 total abundance
  } else if (!is.null(filt.var) & is.null(lev.kept)){
    cat("Error, argument 'lev.kept' has not been provided.\n")
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }

  # Re-organise factor levels
  if (!is.null(factor.lev)){
    for  (i in names(factor.lev)){
      lev.order <- factor.lev[[i]][sort(match(unique(get_variable(data_phylo, i)), factor.lev[[i]]))]
      data_phylo@sam_data[,i] <- factor(get_variable(data_phylo, i), levels = lev.order)
    }
    rownames(data_phylo@sam_data) <- data_phylo@sam_data$SampleID # Rownames lost in the previous process
  } else if (is.null(META.file)){
    cat("No metadata has been provided.\nAre you sure you don't have any?")
  } else if (is.null(factor.lev)){
    cat("Argument 'factor.lev' has not been provided.\nAre you sure your metadata doesn't have any categorical variables?")
  }
  return(data_phylo)
}
