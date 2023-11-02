\name{dimred_covar_cor}
\alias{dimred_covar_cor}
\title{Assess the weight of cell metadata covariates on the components of a dimension reduction components}
\usage{
  dimred_covar_cor(sobj = NULL, dimred = 'pca', ncomp = 10, return_p = FALSE)
}
\arguments{
  \item{sobj}{A \code{Seurat} object to describe.}
  \item{dimred}{The name of the dimred object to assess.}
  \item{ncomp}{The number of \code{dimred} components to assess/plot.}
  \item{return_p}{Returns the correlation test p-values matrix.}
}
\description{
  This function assesses the weight of any cell/barcode metadata covariate on the first \code{ndim} components of a dimension reduction space, through two different tests, depending on the covariate data type : Spearman correlation test for contiunous data, Kruskal-Wallis for categorical. See \code{?pcaExplorer::correlatePCs()}.
}
\value{
  A heatmap displaying the -log10(p.adj) value of each Spearman correlation and/or K-W test, from 0 to 1000 (capped). When \code{return_p} is set to \code{TRUE}, the FDR-adjusted p-value matrix is returned.
}
\details{
  P-values FDR-adjustment is added to the \code{pcaExplorer::correlatePCs} output, using the Benjamini-Hochberg method.
  Infinite -log10(p.adj) values (from p.adj=0) are set to 999.
}
\author{
  Bastien Job
}

