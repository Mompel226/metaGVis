#' @title Filter and produce summary files
#'
#' @description
#' This function is used to filter samples with insufficient read counts, undesired taxa such as order Chloroplast and family Mitochondria and low prevalence taxa.
#' Furthermore, all Phyla with not available/missing taxonomic classification are removed.
#' Filtration at the domain rank can be customised.
#'
#' @param physeq (Required) A phyloseq object.
#' @param file.dir (Optional) A character string of the PATH defining the directory where the summary files and plots are going to be created.
#' If no directory is specified files won't be created. DEFAULT: NULL
#' @param filt.domain (Optional) A character string defining the names of the domains to filter out. DEFAULT: c("Unclassified", "Eukaryota", "Viruses", "Archaea").
#' @param min.depth (Optional) A numerical value -- a threshold -- indicating the minimum number (>=) of reads a sample should have. DEFAULT: 1000
#' @param min.prev (Optional) A numerical value -- a threshold -- indicating the minimum number (>=) of samples in which each taxon should be present. DEFAULT: 2
#' @param min.abund (Optional) A numerical value -- a threshold -- indicating the minimum abundance (>=) a taxon should have. DEFAULT: 10
#' @param abund.type (Optional) A character string defining which type of abundance should be taken into account when filtering ("Total" or "Mean") DEFAULT: "Mean"
#' @param report (Optional) A logical value. If TRUE an excel summary file is created. DEFAULT: FALSE
#' @param seqrows (Optional) A numerical value indicating if rownames are amplicon DNA sequences. If TRUE, a file showing the distribution of sequence lengths between filtered and unfiltered taxa is created. DEFAULT: FALSE.
#'
#' @return
#' A filtered phyloseq object. It also writes summary files in the specified directory.
#'
#' @examples
#' data("GlobalPatterns")
#' colnames(tax_table(GlobalPatterns)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#' To filter out order Chloroplast, family Mitochondria and domains Eukaryota, Archaea, Viruses and Unclassified:
#' data_phylo <- FilterAndSummarise(GlobalPatterns, min.depth = 0, min.prev = 0, min.abund = 0)
#'
#' To perform the previous filtration but also obtain the summary files and carry out an abundance filtration based on total abundance:
#' data_phylo <- FilterAndSummarise(GlobalPatterns, "~/output/summary_files", min.depth = 1000, min.prev = 3, min.abund = 20, abund.type = "Total")
#'
#' @export
FilterAndSummarise <- function(physeq,
                               file.dir = NULL,
                               filt.domain = c("Unclassified", "Eukaryota", "Viruses", "Archaea"),
                               min.depth = 1000,
                               min.prev = 2,
                               min.abund = 10,
                               abund.type = "Total",
                               report = FALSE,
                               seqrows = FALSE){

  # Inspect distribution of sequence lengths before filtering
  if (isTRUE(seqrows)){ dist.seq.length.1 <- table(nchar(dada2::getSequences(t(as.matrix(otu_table(physeq)))))) }
  # Remove Unclassified, Eukaryota, Viruses, Chloroplasts, and Mitochondria from the data.
  # The following also ensures that features with NA phylum annotation are also removed.
  oldMA <- as(tax_table(physeq), "matrix")
  oldDF <- data.frame(oldMA)
  newDF <- subset(oldDF, !Domain %in% filt.domain)
  newMA <- as(newDF, "matrix")
  tax_table(physeq) <- tax_table(newMA)
  physeq %<>% subset_taxa((Order !="Chloroplast") | is.na(Order))
  physeq %<>% subset_taxa((Family != "Mitochondria") | is.na(Family))
  physeq %<>% subset_taxa(!is.na(Phylum))
  clean.physeq <- physeq
  # Remove samples with less than 'min.depth' reads
  physeq = prune_samples(sample_sums(physeq) >= min.depth, physeq)
  prevdf <- calc.prev(physeq)
  # Filter taxa (species/ASVs/OTUs) with less than 'min.abund' reads (on average or in total between all samples)
  filt.t <- (prevdf$Prevalence >= min.prev & prevdf[,abund.type] >= min.abund)
  taxatokeep <- rownames(prevdf)[filt.t]
  physeq <- prune_taxa(taxatokeep, physeq)
  # Remove taxa whose abundance is zero across all samples
  physeq = prune_taxa(taxa_sums(physeq) != 0, physeq)
  # Inspect distribution of sequence lengths after filtering
  if (isTRUE(seqrows)){ dist.seq.length.2 <- table(nchar(dada2::getSequences(t(as.matrix(otu_table(physeq)))))) }
  # Create prevalence data frame with filtered data
  filt.prevdf <- calc.prev(physeq)
  # Subset phyla to create summary plots
  if (!is.null(file.dir)){
    # Create file directory
    ifelse(!dir.exists("DataReports"), dir.create("DataReports"), "")
    for (i in c("before", "after")){
      if (i == "before"){
        phylo = clean.physeq
        df = prevdf
      } else if (i == "after"){
        phylo = physeq
        df = filt.prevdf
      }
      df.P = subset(df, Phylum %in% get_taxa_unique(phylo, "Phylum"))
      plot <- ggplot(df.P, aes(Total, Prevalence / nsamples(phylo), color = Phylum)) +
        facet_wrap(~ Phylum) +
        geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +
        geom_point(size = 2, alpha = 0.7) +
        scale_x_log10() +
        labs(title = paste0("Phylum prevalence ",i, " abundance based filtration."),
             subtitle = "No chloroplast, mitochondria, unclassified phyla and specified domains.",
             x = "Total Abundance (CPM)",
             y = "Prevalence [fraction] (the dotted line indicates a 0.05 prevalence)") +
        theme(legend.position = "none")
      assign(i, plot)
    }
    Combined.plots <- ggpubr::ggarrange(before, after, nrow = 2, ncol = 1, heights = c(1, 0.85))
    # Save plot
    ggsave("DataReports/Phylum_prevalence.png", Combined.plots, height = 12, width = 16, units = "in", dpi = 400)
  }
  # output
  if (!is.null(file.dir) & isTRUE(report)){
    # Convert report into Excel file
    xlsx::write.xlsx(prevdf, paste0(file.dir, "/DataReports/DataSummary.xlsx"), sheetName = "TaxaPrevalence",
                     col.names = TRUE, row.names = TRUE, append = FALSE)
    if (isTRUE(seqrows)){
      xlsx::write.xlsx(dist.seq.length.1, paste0(file.dir, "/DataReports/DataSummary.xlsx"), sheetName = "SeqDistBeforeFilt",
                       col.names = TRUE, row.names = FALSE, append = TRUE)
      xlsx::write.xlsx(dist.seq.length.2, paste0(file.dir, "/DataReports/DataSummary.xlsx"), sheetName = "SeqDistAfterFilt",
                       col.names = TRUE, row.names = FALSE, append = TRUE)
    }
  }
  return(physeq)
}
################################################################################################################################################################################################################################
#' @title Calculate feature prevalence.
#'
#' @description
#' Function to create a data frame showing the prevalence of each feature across all samples plus extra infromation.
#'
#' @param physeq (Required) A phyloseq object.
#'
#' @return
#' A data frame.
#'
#' @examples
#' df <- calc.prev(phyloseq.object)
#'
#' @export
calc.prev <- function(physeq){
  # Due to different library sizes, we have to normalise the data to have a better look
  # Convert to CPM
  physeq.cpm <- physeq
  physeq.cpm@otu_table <- otu_table(
    apply(X = otu_table(physeq.cpm),
          MARGIN = 2,
          FUN = function(x){ (x/max(sum(x), 1e-32))*1e6 }), taxa_are_rows = TRUE
  )
  # The following R code was adapted from:
  # https://bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/Microbiome/MicrobiomeWorkflowII.html#filtering
  # Compute prevalence of each feature, store as data.frame
  prev <- apply(X = otu_table(physeq.cpm),
                  MARGIN = 1,
                  FUN = function(x){sum(x > 0)})
  # Add taxonomy and total read counts to this data.frame
  df <- as.data.frame(otu_table(physeq.cpm))
  prevdf <- data.frame(Prevalence = prev,
                       Total = base::rowSums(df),
                       Mean = apply(df, 1, function(x){ mean(x[x!=0])}), #Mean of non-zero values
                       tax_table(physeq.cpm),
                       otu_table(physeq.cpm))
  # Order taxa by Mean Abundance
  prevdf %<>% dplyr::arrange(dplyr::desc(Mean))
  return(prevdf)
}
