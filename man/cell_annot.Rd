\name{cell_annot}
\alias{cell_annot}
\title{Performs automatic cell type annotation of a Seurat object, at the cell or cluster level, using SingleR and celldex databases}
\usage{
  cell_annot(sobj = NULL, assay = 'RNA', slot = 'data', celldex_setname = 'ImmGenData', group_by = NULL, list_celldex = FALSE)
}
\arguments{
  \item{sobj}{A \code{Seurat} object to describe.}
  \item{assay}{The \code{Assay} to extract data from}
  \item{slot}{The \code{slot} from the requested \code{Assay} to use. Should be one of 'data', 'scale.data'.}
  \item{celldex_setname}{Name of the \code{celldex} database to use for annotation.}
  \item{group_by}{Name of one metadata column factor to perform the annotation at a group level (by exemple, clusters)}.
  \item{list_celldex}{Set to \code{TRUE} to get the list of available \code{celldex} databases .}
}
\description{
  This function performs an automatic annotation of cells/barcodes from a Seurat object according to a \code{celldex} annotation database.
}
\details{
  Any slot from any assay can be used.
  Annotation is by default performed at the cell/barode level, but can be set according to groups of cells/barcodes using \code{group_by}.
  Use \code{list_celldex = TRUE} to get the list of available \code{celldex} databases.
}
\value{
  A factor vector of celldex labels, corresponding to cells/barcodes in the same order as in the Seurat object.
}
\author{
  Bastien Job
}

