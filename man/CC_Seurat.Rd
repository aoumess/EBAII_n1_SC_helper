\name{CC_Seurat}
\alias{CC_Seurat}
\title{Compute cell-cyle state/score prediction, using the Seurat method}
\usage{
  CC_Seurat(sobj = NULL, assay = 'RNA', seurat_cc_genes = NULL, SmG2M = TRUE, nbin = 24, my_seed = 1337)
}
\arguments{
  \item{sobj}{A \code{Seurat} object to describe.}
  \item{assay}{The \code{Assay} to extract data from.}
  \item{seurat_cc_genes}{A list of two named vectors ('s.genes', 'g2m.genes'), containing the Seurat cell cycle genes. See \code{?Seurat::CellCycleScoring()}.}
  \item{SmG2M}{Compute the "S minus G2M" score.}
  \item{nbin}{Number of bins. See \code{?Seurat::AddModuleScore()}.}
  \item{my_seed}{Seed used for any step using RNG.}
}
\description{
  This function computes the prediction of cell cycle phase through scores from Seurat object, using the Seurat method. See \code{?Seurat::CellCycleScoring()}.
}
\value{
  A modified Seurat object with additional barcode metadata annotation : 'CC_Seurat_S.Score' (S-phase score), 'CC_Seurat_G2M.Score' (G2/M-phase score), 'CC_Seurat_SmG2M.Score' (difference between the S and G2/M scores) and 'CC_Seurat_Phase' (the cell cyle phase attribute).
}
\author{
  Bastien Job
}

