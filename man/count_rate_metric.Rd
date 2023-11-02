\name{count_rate_metric}
\alias{count_rate_metric}
\title{Easy way to peform raw counts rate metrics on a Seurat object}
\usage{
  count_rate_metric(sobj = NULL, assay = 'RNA', features = NULL)
}
\arguments{
  \item{sobj}{A \code{Seurat} object to describe.}
  \item{assay}{The \code{Assay} to extract data from.}
  \item{features}{ Vector of features names to use to perform their count rate per barcode.}
}
\description{
  This function allows to easily compute count rate metrics at the cell/barcode level for a set of \code{features}.
}
\value{
  A vector of count rates for the defined \code{features}, of same length as the number of cells/barcodes in the Seurat object \code{sobj}.
}
\details{
  This is mainly used to compute the \%mito, \%ribo, \%stress (...) metrics.
}
\author{
  Bastien Job
}

