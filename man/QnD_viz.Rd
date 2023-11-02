\name{QnD_viz}
\alias{QnD_viz}
\title{Quick and dirty visualization of a Seurat dataset as a UMAP/TSNE}
\usage{
  QnD_viz(sobj = NULL, assay = 'RNA', slot = 'counts', dimred = NULL, reduction_method1 = 'pca', reduction_method2 = 'umap', nfeatures = 2000, vtr = NULL, ncomp = 30, my_seed = 1337, group_by = NULL, features = NULL, pt_size = 1, dark_theme = TRUE, return_object = FALSE)
}
\arguments{
  \item{sobj}{A \code{Seurat} object to describe.}
  \item{assay}{The \code{Assay} to extract data from.}
  \item{slot}{The \code{slot} from the requested \code{Assay} to extract data from. Has to be set to  \code{NULL} if the  \code{dimred} parameter is used.}
  \item{dimred}{The  \code{Reduction} to use if it does exist in the  \code{Seurat} object. Can't be used if  \code{slot} is not \code{NULL}.}
  \item{reduction_method1}{Type of dimension reduction method to use first. Only\code{pca} is currently supported. Only used in \code{slot} mode, not when \code{dimred} is used.}
  \item{reduction_method2}{Type of dimension reduction method to use in second step. Can be one of \code{umap} or \code{tsne}.}
  \item{nfeatures}{Number of 'High Variation Features' to use when scaling of data is performed. Only used in \code{slot} mode, not when \code{dimred} is used.}
  \item{vtr}{Name of covariates to regress at the scaling step. Corresponds to barcode annotations. Only used in \code{slot} mode, not when \code{dimred} is used.}
  \item{ncomp}{Number of components from the 1st dimension reduction to perform the 2nd one. Note that the number of components generated for the 1st dimension reduction will be \code{2 * ncomp}. Only used in \code{slot} mode, not when \code{dimred} is used.}
  \item{my_seed}{Seed used for any step using RNG.}
  \item{group_by}{Name of one (or more) metadata columns to group (color) barcodes. See \code{?Seurat::DimPlot()}.}
  \item{features}{Vector of feaures to plot. It can be a feature (ie, a gene name), the name of a metadata column for continous data. See \code{?Seurat::FeaturePlot()}.}
  \item{pt_size}{Size of a dot in plots.}
  \item{dark_theme}{Set to \code{TRUE} to use Seurat's dark theme for plots. Recommended for better contrast. See \code{?Seurat::DarkTheme()}.}
  \item{return_object}{Set to \code{TRUE} to return the modified \code{sobj} Seurat object with all additional entries computed to generate the visualization.}
}
\description{
  This function performs a visualization of a Seurat object content as UMAP or TSNE, from any kind of data : raw counts, normalized, or scaled from an assay (= 'slot' mode), or from a 1st step dimension reduction (= 'dimred' mode). When starting from raw or normalized counts, control is given for the amount of features to use and regression of covariates, at the scaling step. For both 'modes', control is given in the amount of 1st dimred components to use for the second, using \code{ndim}. For the output plot, features of \code{DimPlot} can be used if \code{group_by} is used, or \code{FeaturePlot} if \code{features} is used. The modified \code{sobj} can be returned when using \code{return_object = TRUE}.
}
\value{
  By default, only the UMAP/TSNE plot. When using \code{return_object = TRUE}, the modified object is also returned.
}
\author{
  Bastien Job
}

