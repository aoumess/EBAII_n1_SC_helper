\name{SoupX_auto}
\alias{SoupX_auto}
\title{Easy way to run SoupX in automatic mode}
\usage{
  SoupX_auto(sobj = NULL, assay = 'RNA', scmat_raw = NULL, soupQuantile = 0.9, contaminationRange = c(.01, .8), soupRange = c(0,100), return_object = FALSE, doPlot = FALSE) }
\arguments{
  \item{sobj}{A \code{Seurat} object to describe.}
  \item{assay}{The \code{Assay} to extract data from.}
  \item{scmat_raw}{An integer matrix (or sparse matrix) of raw counts, corresponding to unfiltered single cell data (ie, before empty droplets removal).}
  \item{soupQuantile}{See \code{?SoupX::autoEstCont}.}
  \item{contaminationRange}{See \code{?SoupX::autoEstCont}.}
  \item{soupRange}{Range of raw counts in which the estimation of the soup will be performed. See \code{?SoupX::estimateSoup}.}
  \item{return_object}{Set to \code{TRUE} to return the modified \code{sobj} Seurat object with soup removed.}
  \item{doPlot}{Set to \code{TRUE} to get the SoupX plot. See \code{?SoupX::autoEstCont}.}
}
\description{
  This function allows to easily perform soup quantification and/or removal using \code{SoupX}, while keeping control on its main parameters.
}
\value{
  By default, only the \code{rho} value is returned (the soup fraction estimation). When using \code{return_object = TRUE}, the \code{Seurat} object from which soup counts have been removed is also returned.
}
\details{
  The use of a \code{scmat_raw} corresponds to the most efficient way to use SoupX, so recommended when possible.
}
\author{
  Bastien Job
}

