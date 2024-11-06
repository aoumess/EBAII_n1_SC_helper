\name{assess_covar}
\alias{assess_covar}
\title{Assess the weight of cell metadata covariates on the components of a dimension reduction components}
\usage{
  assess_covar(mat = NULL, covar.df = NULL, markers = NULL, ndim = 10L, center = TRUE, scale = TRUE)
}
\arguments{
  \item{mat}{An feature-by-cell expression \code{matrix}.}
  \item{covar.df}{A \code{data.frame} containing the covariates to assess.}
  \item{markers}{A vector of feature names corresponding to feature markersto assess.}
  \item{ndim}{Number of dimensions in which to assess covariates.}
  \item{center}{Center the matrix ?}
  \item{scale}{Scale the matrix ?}
}
\description{
  This function assesses the weight of any cell/barcode metadata covariate(s) on the first \code{ndim} components of a dimension reduction space, through two different tests, depending on the covariate data type : Spearman correlation test for continuous data (numerics), Kruskal-Wallis for categorical (fators).
}
\value{
  A heatmap displaying the weight of each Spearman correlation and/or K-W test, from 0 to 1.
}
\author{
  Bastien Job
}

