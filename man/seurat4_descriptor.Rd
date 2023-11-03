\name{seurat4_descriptor}
\alias{seurat4_descriptor}
\title{Quick description of a Seurat object content}
\usage{
  seurat4_descriptor(sobj = NULL, describe = 'all', 
  sparse_level = TRUE, max_levels = 100)
}
\arguments{
  \item{sobj}{A \code{Seurat} object to describe.}
  \item{describe}{The type of entry/slot from the \code{Seurat} object to describe. Can be a combination of : 'assay' (to describe only \code{Assays}), 'dimred' (to describe only \code{Reductions}), or 'coldata' (to describe only barcode descriptors), or 'all' (to get the whole three types at once).}
  \item{sparse_level}{Compute the sparsity level of integer count matrices (can be disabled as a bit long for very large matrices).}
  \item{max_levels}{For barcode annotations, any column will be considered by default as a factor. Any column that contains at least the \code{max_levels} amount of levels will be processed as continuous data. This is an attempt to handle different data types automatically.}
}
\description{
  This function performs a detailed by quick description of the content of a Seurat (version <= 4).
}
\value{
  Nothing will be returned, all description will be written on the standard output.
}
\details{
  The default assay ("active" assay, as per Seurat's vocabulary), will be tagged with the '[ACTIVE]' word in its description.
}
\author{
  Bastien Job
}

