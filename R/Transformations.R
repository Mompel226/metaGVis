###############################################################################################################################################################################
#' @title Hellinger transformation of phyloseq object.
#'
#' @description
#' This function performs the Hellinger transformation on the abundance data of a phyloseq object.
#' This transformation should be used with raw abundance data.
#'
#' @param physeq (Required) A phyloseq object.
#'
#' @return
#' A phyloseq object with transformed abundance data.
#'
#' @examples
#' data("GlobalPatterns")
#' transformed.data <- hellinger(GlobalPatterns)
#'
#' @export
hellinger <- function(physeq){
  x <- physeq@otu_table
  xt <- vegan::decostand(x, method = "hellinger", MARGIN = 2)
  physeq@otu_table <- otu_table(xt, taxa_are_rows = TRUE)
  return(physeq)
}
###############################################################################################################################################################################
#' @title Chord transformation of phyloseq object (pseudocount of 1 added only to values = 0).
#'
#' @description
#' This function performs the Chord transformation on the abundance data of a phyloseq object.
#' The pseudocount used to avoid 0s when log transforming is 1, however, it is only added to values equal to 0.
#' This transformation should be used with raw abundance data.
#'
#' @param physeq (Required) A phyloseq object.
#'
#' @return
#' A phyloseq object with transformed abundance data.
#'
#' @examples
#' data("GlobalPatterns")
#' transformed.data <- logchord_0_plus_1(GlobalPatterns)
#'
#' @export
logchord_0_plus_1 <- function(physeq){
  x <- physeq@otu_table
  x[ x == 0 ] <- 1
  x <- log(x)
  xt <- vegan::decostand(x, method = "norm", MARGIN = 2)
  physeq@otu_table <- otu_table(xt, taxa_are_rows = TRUE)
  return(physeq)
}
###############################################################################################################################################################################
#' @title Chord transformation of phyloseq object (pseudocount of 1 added to all values).
#'
#' @description
#' This function performs the Chord transformation on the abundance data of a phyloseq object.
#' The pseudocount used to avoid 0s when log transforming is 1.
#' This transformation should be used with raw abundance data.
#'
#' @param physeq (Required) A phyloseq object.
#'
#' @return
#' A phyloseq object with transformed abundance data.
#'
#' @examples
#' data("GlobalPatterns")
#' transformed.data <- logchord_1p(GlobalPatterns)
#'
#' @export
logchord_1p <- function(physeq){
  x <- physeq@otu_table
  x <- log1p(x)
  xt <- vegan::decostand(x, method = "norm", MARGIN = 2)
  physeq@otu_table <- otu_table(xt, taxa_are_rows = TRUE)
  return(physeq)
}
###############################################################################################################################################################################
#' @title CLR (centered log-ratio) transformation of phyloseq object (pseudocount equal to minimum_value/2).
#'
#' @description
#' This function performs the Chord transformation on the abundance data of a phyloseq object.
#' The pseudocount used to avoid 0s when log transforming is equal to the minimum abundance value divided by 2.
#' This transformation should be used with raw abundance data (internal Total Sum Scaling transformation -- relative abundances).
#'
#' @param physeq (Required) A phyloseq object.
#'
#' @return
#' A phyloseq object with transformed abundance data.
#'
#' @examples
#' data("GlobalPatterns")
#' transformed.data <- clr_min_percent(GlobalPatterns)
#'
#' @export
clr_min_percent <- function(physeq){
  physeq <- microbiome::transform(physeq, transform = "clr", target = "sample")
  return(physeq)
}
###############################################################################################################################################################################
#' @title CLR (centered log-ratio) transformation of phyloseq object (pseudocount of 1 added only to values = 0).
#'
#' @description
#' This function performs the Chord transformation on the abundance data of a phyloseq object.
#' The pseudocount used to avoid 0s when log transforming is 1, however, it is only added to values equal to 0.
#' This transformation should be used with raw abundance data (internal Total Sum Scaling transformation -- relative abundances).
#'
#' @param physeq (Required) A phyloseq object.
#'
#' @return
#' A phyloseq object with transformed abundance data.
#'
#' @examples
#' data("GlobalPatterns")
#' transformed.data <- clr_min_percent(GlobalPatterns)
#'
#' @export
clr_0_plus_1 <- function(physeq){
  ps <- physeq
  x <- physeq@otu_table
  xt <- apply(x, 2, function(x){x/sum(x)}) # Transform to relative abundance
  # The clr function from the compositions package only performs log on non-zero values!!!
  ps@otu_table <- otu_table(t(compositions::clr(t(xt))), taxa_are_rows = TRUE)
  colnames(ps@otu_table) <- colnames(physeq@otu_table)
  return(ps)
}
###############################################################################################################################################################################
#' @title CLR (centered log-ratio) transformation of phyloseq object (pseudocount equal to random numbers between 0.1*min_value and min_value).
#'
#' @description
#' This function performs the Chord transformation on the abundance data of a phyloseq object.
#' The pseudocounts used to avoid 0s when log transforming are random numbers between the minimum abundance value (min_value) and 0.1*min_value.
#' This transformation should be used with raw abundance data (internal Total Sum Scaling transformation -- relative abundances).
#'
#' @param physeq (Required) A phyloseq object.
#'
#' @return
#' A phyloseq object with transformed abundance data.
#'
#' @references
#' Lubbe S. et al.  (2021) Comparison of zero replacement strategies for compositional data with large numbers of zeros. Chemometr. Intell. Lab., 210, 104248.
#'
#' @examples
#' data("GlobalPatterns")
#' transformed.data <- clr_runif(GlobalPatterns)
#'
#' @export
clr_runif <- function(physeq){
  ps <- physeq
  x <- physeq@otu_table
  xt <- apply(x, 2, function(x){x/max(sum(x), 1e-32)}) # Transform to relative abundance
  if (any(xt == 0)) {
    v <- as.vector(xt)
    minval <- min(v[v > 0]) # Extract the minimum value superior to 0
    Nof0 <- sum(v == 0) # Obtain the number of 0s
    xt[xt == 0] <- runif(Nof0, min = 0.1*minval, max=minval) # Obtain random numbers between min and max, and replace 0s
  }
  d <- apply(xt, 2, function(x) {log(x) - mean(log(x))}) # CLR transformation
  ps@otu_table <- otu_table(d, taxa_are_rows = TRUE) # Change otu_table from original phyloseq object
  return(ps)
}
