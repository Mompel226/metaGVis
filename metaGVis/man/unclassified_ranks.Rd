% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unclassified_ranks.R
\name{unclassified_ranks}
\alias{unclassified_ranks}
\title{Unclassified ranks}
\usage{
unclassified_ranks(
  physeq,
  rank,
  word,
  format = "short",
  agglomerate = FALSE,
  count = FALSE,
  output = "physeq"
)
}
\arguments{
\item{physeq}{A phyloseq object}

\item{rank}{A character string defining the taxonomic rank in which the names will be modified}

\item{word}{A character string defining the name given to n.a values found in the specified taxonomic rank.
Default: "N.A"}

\item{format}{A character string defining the format in which the new names will written: "short" or "long". (for more info look @description)
Default: "short"}

\item{agglomerate}{A logical value.
If TRUE the species that have the same taxonomy (after renaming) at the specified taxonomic rank will be merged.
Default: FALSE}

\item{count}{A logical value.
If TRUE the number of species that have been agglomerated is appended at the end of the taxonomic name.
Default: "FALSE"}

\item{output}{A character string defining the format of the output: "physeq" or "taxnames". (for more info look @description)
Default: "physeq"}
}
\value{
If the "output" parameter is set to "physeq" --> the function returns a phyloseq object with a modified tax_table
If the "output" parameter is set to "taxnames" --> the function returns a named vector that uses the tax_table rownames as names.
The "taxnames" output is not available when using the agglomerate option.
}
\description{
This function is used to rename unclassified/n.a taxa using the information their higher taxonomic ranks.
In case of having an unclassified genera, it will replace it by something like: Unidentified (Family: ...).
If the Family rank is also unclassified then it will replace it by Unidentified (Order: ...) if using the short mode
or Unidentified (Order: ..., Family: N.A) if using the long mode.
It can also agglomerate species at the specified taxonomic rank and show the number of agglomerated species: Unidentified (Order: ...)(nº)
The main input of this function is a phyloseq object and the output can be a phyloseq object with a modified tax_table
or a named vector that uses the tax_table rownames as names (not available when using the agglomerate option).
}
\examples{
Obtaining a non-agglomerated phyloseq object:
unclassified_ranks(physeq=phyloseq.object, rank="Genus", format="short, output="physeq")

Obtaining an agglomerated phyloseq object:
unclassified_ranks(physeq=phyloseq.object, rank="Genus", format="short, count = TRUE, agglomerate = TRUE, output="physeq")
}
