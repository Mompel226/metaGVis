#' @title Unclassified ranks
#'
#' @description
#' This function is used to rename unclassified/n.a taxa using the information of their higher taxonomic ranks.
#' In case of having an unclassified genera, it will replace it by something like: Unidentified (Family: ...).
#' If the Family rank is also unclassified then it will replace it by Unidentified (Order: ...) if using the short mode
#' or Unidentified (Order: ..., Family: N.A) if using the long mode.
#' It can also agglomerate species at the specified taxonomic rank and show the number of agglomerated species: Unidentified (Order: ...)(nº)
#' The main input of this function is a phyloseq object and the output can be a phyloseq object with a modified tax_table
#' or a named vector that uses the tax_table rownames as names (not available when using the agglomerate option).
#'
#' @param physeq A phyloseq object
#' @param rank A character string defining the taxonomic rank in which the names will be modified
#' @param word A character string defining the name given to n.a values found in the specified taxonomic rank.
#' Default: "N.A"
#' @param format A character string defining the format in which the new names will written: "short" or "long". (for more info look @description)
#' Default: "short"
#' @param agglomerate A logical value.
#' If TRUE the species that have the same taxonomy (after renaming) at the specified taxonomic rank will be merged.
#' Default: FALSE
#' @param count A logical value.
#' If TRUE the number of species that have been agglomerated is appended at the end of the taxonomic name.
#' Default: "FALSE"
#' @param output A character string defining the format of the output: "physeq" or "taxnames". (for more info look @description)
#' Default: "physeq"
#'
#' @return
#' If the "output" parameter is set to "physeq" --> the function returns a phyloseq object with a modified tax_table
#' If the "output" parameter is set to "taxnames" --> the function returns a named vector that uses the tax_table rownames as names.
#' The "taxnames" output is not available when using the agglomerate option.
#'
#' @examples
#' Obtaining a non-agglomerated phyloseq object:
#' unclassified_ranks(physeq=phyloseq.object, rank="Genus", format="short, output="physeq")
#'
#' Obtaining an agglomerated phyloseq object:
#' unclassified_ranks(physeq=phyloseq.object, rank="Genus", format="short, count = TRUE, agglomerate = TRUE, output="physeq")
#' @export
unclassified_ranks <- function(physeq, rank, word, format = "short", agglomerate = FALSE, count = FALSE, output = "physeq"){

  stop_quietly <- function() {
    on.exit(options(show.error.messages = FALSE))
    stop()
  }

  if(class(physeq) != "phyloseq"){
    cat("Error, argument 'physeq' is not of class phyloseq.\n")
    stop_quietly()
  }

  if(!rank %in% c ("Species", "Genus", "Family", "Order", "Class", "Phylum", "Domain")){
    cat("Specified rank is incorrect. Check spelling or rank names.\nTaxonomic names in this function are defined as:\nDomain-Phylum-Class-Order-Family-Genus-Species\n")
    stop_quietly()
  }

  if (!format %in% c("short", "long")){
    cat("Specified format option is incorrect.\nThe only options available are: 'short' or 'long'")
    stop_quietly()
  }

  if(rank == "Species"){
    # Merge genus and species names if species name doesn't start with uppercase
    ifelse(!grepl("[[:upper:]]", substr(physeq@tax_table[,"Species"], 1, 1)),
           physeq@tax_table[,"Species"]<-paste0(physeq@tax_table[,"Genus"], " ", physeq@tax_table[,"Species"]), "")
    if (output == "physeq"){
      return(physeq)
    } else {
      taxnames = physeq@tax_table[,rank]
      return(taxnames)
    }
  } else if(rank == "Genus"){
    upper_rank <- c ("Family", "Order", "Class", "Phylum", "Domain")
  } else if (rank == "Family"){
    upper_rank <- c ("Order", "Class", "Phylum", "Domain", "NA")
  } else if (rank == "Order"){
    upper_rank <- c ("Class", "Phylum", "Domain", "NA", "NA")
  } else if (rank == "Class"){
    upper_rank <- c ("Phylum", "Domain", "NA", "NA", "NA")
  } else if (rank == "Phylum"){
    upper_rank <- c ("Domain", "NA", "NA", "NA", "NA")
  }

  # Some tools use n.a for unclassifed ranks other use "unclassified" or similar strings.
  # This step allows to define a string character for n.a values
  if(missing(word)){
    word = "N.A"
    for (i in colnames(physeq@tax_table)){
      physeq@tax_table[,i][is.na(physeq@tax_table[,i])] <- word
    }
  }

  if (rank != "Domain"){
    # Check if any taxnames are "unclassified" and change them
    taxnames = physeq@tax_table[,rank]
    if(any(taxnames == word)){
      NAs_upperclass1 <- physeq@tax_table[,upper_rank[1]][physeq@tax_table[,rank] == word]
      if(format == "short"){
        name = paste0("Unidentified (", upper_rank[1], ": ", NAs_upperclass1, ")")
      } else {
        name = paste0("Unidentified ", rank, " (",upper_rank[1], ": ", NAs_upperclass1, ")")
      }
      physeq@tax_table[,rank][physeq@tax_table[,rank] == word] <- name
      taxnames = physeq@tax_table[,rank]
    }

    # If by any chance two levels are identical - e.g.Genus "unclassified" and same Family, then add Order:
    if(format != "short"){
      identifier = paste0("Unidentified ", rank, " (", upper_rank[1], ": N.A)")
    } else {
      identifier = paste0("Unidentified (", upper_rank[1], ": N.A)")
    }
    if (!any(grep("Domain: N.A", identifier))){
      if(any(taxnames == identifier)){
        NAs_upperclass2 <- physeq@tax_table[,upper_rank[2]][physeq@tax_table[,rank] == identifier]
        NAs_upperclass1v2 <- physeq@tax_table[,upper_rank[1]][physeq@tax_table[,rank] == identifier]
        if(format == "short"){
          name = paste0("Unidentified (", upper_rank[2], ": ", NAs_upperclass2, ")")
        } else {
          name = paste0("Unidentified ", rank, " (", upper_rank[2], ": ", NAs_upperclass2, "; ", upper_rank[1], ": ", NAs_upperclass1v2, ")")
        }
        physeq@tax_table[,rank][physeq@tax_table[,rank] == identifier] <- name
        taxnames = physeq@tax_table[,rank]
      }
    }

    # If by any chance three levels are identical - e.g.Genus "unclassified" and same for Family and Order, then add Class:
    if(format != "short"){
      identifier = paste0("Unidentified ", rank, " (", upper_rank[2], ": N.A; ", upper_rank[1], ": N.A)")
    } else {
      identifier = paste0("Unidentified (", upper_rank[2], ": N.A)")
    }
    if (!any(grep("Domain: N.A", identifier))){
      if(any(taxnames == identifier)){
        NAs_upperclass3 <- physeq@tax_table[,upper_rank[3]][physeq@tax_table[,rank] == identifier]
        NAs_upperclass1v3 <- physeq@tax_table[,upper_rank[1]][physeq@tax_table[,rank] == identifier]
        NAs_upperclass2v1 <- physeq@tax_table[,upper_rank[2]][physeq@tax_table[,rank] == identifier]
        if(format == "short"){
          name = paste0("Unidentified (", upper_rank[3], ": ", NAs_upperclass3, ")")
        } else {
          name = paste0("Unidentified ", rank, " (", upper_rank[3], ": ", NAs_upperclass3, "; ", upper_rank[2], ": ", NAs_upperclass2v1, "; ", upper_rank[1], ": ", NAs_upperclass1v3, ")")
        }
        physeq@tax_table[,rank][physeq@tax_table[,rank] == identifier] <- name
        taxnames = physeq@tax_table[,rank]
      }
    }

    # If by any chance four levels are identical - e.g.Genus "unclassified" and same for Family, Order and Class, then add Phylum:
    if(format != "short"){
      identifier = paste0("Unidentified ", rank, " (", upper_rank[3], ": N.A; ", upper_rank[2], ": N.A; ", upper_rank[1], ": N.A)")
    } else {
      identifier = paste0("Unidentified (", upper_rank[3], ": N.A)")
    }
    if (!any(grep("Domain: N.A", identifier))){
      if(any(taxnames == identifier)){
        NAs_upperclass4 <- physeq@tax_table[,upper_rank[4]][physeq@tax_table[,rank] == identifier]
        NAs_upperclass1v4 <- physeq@tax_table[,upper_rank[1]][physeq@tax_table[,rank] == identifier]
        NAs_upperclass2v2 <- physeq@tax_table[,upper_rank[2]][physeq@tax_table[,rank] == identifier]
        NAs_upperclass3v1 <- physeq@tax_table[,upper_rank[3]][physeq@tax_table[,rank] == identifier]
        if(format == "short"){
          name = paste0("Unidentified (", upper_rank[4], ": ", NAs_upperclass4, ")")
        } else {
          name = paste0("Unidentified ", rank, " (", upper_rank[4], ": ", NAs_upperclass4, "; ", upper_rank[3], ": ", NAs_upperclass3v1, "; ", upper_rank[2], ": ", NAs_upperclass2v2, "; ", upper_rank[1], ": ", NAs_upperclass1v4, ")")
        }
        physeq@tax_table[,rank][physeq@tax_table[,rank] == identifier] <- name
        taxnames = physeq@tax_table[,rank]
      }
    }

    # If by any chance five levels are identical - e.g.Genus "unclassified" and same for Family, Order, Class and Phylum, then add Domain:
    if(format != "short"){
      identifier = paste0("Unidentified ", rank, " (", upper_rank[4], ": N.A; ", upper_rank[3], ": N.A; ", upper_rank[2], ": N.A; ", upper_rank[1], ": N.A)")
    } else {
      identifier = paste0("Unidentified (", upper_rank[4], ": N.A)")
    }
    if (!any(grep("Domain: N.A", identifier))){
      if(any(taxnames == identifier)){
        NAs_upperclass5 <- physeq@tax_table[,upper_rank[5]][physeq@tax_table[,rank] == identifier]
        NAs_upperclass1v5 <- physeq@tax_table[,upper_rank[1]][physeq@tax_table[,rank] == identifier]
        NAs_upperclass2v3 <- physeq@tax_table[,upper_rank[2]][physeq@tax_table[,rank] == identifier]
        NAs_upperclass3v2 <- physeq@tax_table[,upper_rank[3]][physeq@tax_table[,rank] == identifier]
        NAs_upperclass4v1 <- physeq@tax_table[,upper_rank[4]][physeq@tax_table[,rank] == identifier]
        if(format == "short"){
          name = paste0("Unidentified (", upper_rank[5], ": ", NAs_upperclass5, ")")
        } else {
          name = paste0("Unidentified ", rank, " (", upper_rank[5], ": ", NAs_upperclass5, "; ", upper_rank[4], ": ", NAs_upperclass4v1, "; ", upper_rank[3], ": ", NAs_upperclass3v2, "; ", upper_rank[2], ": ", NAs_upperclass2v3, "; ", upper_rank[1], ": ", NAs_upperclass1v5, ")")
        }
        physeq@tax_table[,rank][physeq@tax_table[,rank] == identifier] <- name
        taxnames = physeq@tax_table[,rank]
      }
    }
  }
  # Agglomerate
  if (isTRUE(agglomerate)){
    agg.physeq <- phyloseq::tax_glom(physeq, taxrank = rank, NArm=FALSE)
    if (isTRUE(count) & rank != "Species"){
      n.rank  <- which(rank_names(physeq) %in% rank)
      # Concatenate taxa names up to one taxrank column in case some names are identical in the specified rank but not in the upper one.
      if(rank != "Domain"){
        concat.taxa <- physeq@tax_table[,c(n.rank-1, n.rank)]
        ASV.preval <- table(paste0(concat.taxa[,1],";", concat.taxa[,2]))
      } else {
        concat.taxa <- physeq@tax_table[,c(n.rank)]
        ASV.preval <- table(paste0(concat.taxa[,1]))
      }
      # table of the counts
      name_with_counts <- as.vector(paste0(rownames(ASV.preval)," (",ASV.preval,")"))
      name_with_counts <- sort(gsub("^.*;", "", name_with_counts)) # remove upper rank
      ord.data <- agg.physeq@tax_table[,rank][order(agg.physeq@tax_table[,rank])]
      # make sure taxa names match
      for (n in make.unique(name_with_counts)){
        k <- sub("\\s+[^ ]+$", "", n)
        ord.data[grepl(gsub("[()]", "", k), gsub("[()]", "", ord.data))][1] <- n # remove parenthesis for grepl to work
      }
      agg.physeq@tax_table[,rank][order(agg.physeq@tax_table[,rank])] <- ord.data
    }
  }

  if (output == "physeq" & !isTRUE(agglomerate)){
    return(physeq)
  } else if (output == "physeq" & isTRUE(agglomerate)){
    return(agg.physeq)
  } else if (output == "taxnames" & isTRUE(agglomerate)){
    cat("The output format 'taxnames' is not available when agglomerate = TRUE\n")
    stop()
  } else if (output == "taxnames" & !isTRUE(agglomerate)){
    return(taxnames)
  }
}
