# SC.helper

A set of functions to help the analysis of single cell RNAseq data from a Seurat object.

Mainly targeted at the EBAII n1 single-cell RNAseq training course.

May help whoever dares to try it.

-  **SeuratObject_descriptor** : Describes a summary of the content of a SeuratObject (v4/v5) : assays and layers, cells meta.data, dimension reductions

-  **QnD_viz** : Creates a DimPlot or FeaturePlot from a Seurat object at any step of the analysis (from a Seurat object with only a raw unfiltered counts sparsematrix to a 2nd level dimension reduction), by applying the missing steps to the UMAP generation using Seurat defaults.

-  **SoupX_auto** : Wrapper to SoupX, applies the method to either a single filtered matrix or the unfiltered+filter pair. Using either the automatic mode or the fraction counts mode.

-  **CC_Seurat** : One-liner to run the Seurat implementation of cell-cyle phase status estimation and scores.

-  **CC_Cyclone** : One-liner to run the scran's Cyclone implementation of cell-cyle phase status estimation and scores.

-  **assess_covar** : Function to assess and plot the weight of covariates (and marker genes exrpession as control) on a matrix, through a PCA reduction. Helps at identifying adverse covariates to regress while maintaining the expected expression of cell markers.

-  **cell_annot** : Function to ease the automatic annotation of cell types using SingleR and a celldex reference database.

