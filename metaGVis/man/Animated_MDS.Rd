% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3D_plots.R
\name{Animated_MDS}
\alias{Animated_MDS}
\title{3D plots}
\usage{
Animated_MDS(
  plotly_object,
  title,
  path_to_conda = "/Users/daniel/anaconda3/bin/"
)
}
\arguments{
\item{plotly_object}{A phyloseq object.}

\item{title}{A character string defining the title of the plot.}

\item{path_to_conda}{A character string defining the PATH to the bin folder of the Conda package manager}
}
\value{
A gif file.
}
\description{
This function is used to create an animation from a 3D ordination plot created with plotly.
}
\examples{
Animated_MDS(plotly_object = plot_3D, title = "PCoA_bray", path_to_conda = "/Users/John/anaconda3/bin/")

}
