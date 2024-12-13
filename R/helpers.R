## Seurat-descriptor
## TESTED ON : Seurat=4.4.0,5.1.0, sparseMatrixStats=1.13.4 , knitr=.44
## sobj           [Seurat]    A Seurat object
## max_levels     [int>0]     Maximal number of unique values to consider a barcode annotation as a factor rather than a continuous numeric vector
## describe       ['all', 'assays', 'dimred', 'coldata']    Type of entries to describe
SeuratObject_descriptor <- function(sobj = NULL, describe = 'all', sparse_level = TRUE, max_levels = 100) {
  
  ## Checks
  if (is.null(sobj)) stop('A Seurat object is required !')
  if(!is(sobj, 'Seurat')) stop('Provided sobj is not a proper Seurat object !')
  
  describe <- tolower(describe)
  desc.valids <- c('all', 'assay', 'dimred', 'coldata')
  if(!all(describe %in% desc.valids)) stop('At least one requested item type to describe is not valid. Expecting any combination of ["', paste(c(desc.valids), collapse = '", "'), '] ("all" supersedes any other).')
  if('all' %in% describe) describe <- desc.valids[-1]
  
  ## Additional checks on sobj
  if(!is(sobj, 'Seurat')) stop('Provided object is not a proper Seurat object !')
  
  retlist <- list()
  
  ## Seurat slots
  slots <- c('counts', 'data', 'scale.data')
  
  ## Getting projectname
  # projectname <- sobj@project.name
  
  ### OBJECT VERSION
  cat(paste0('OBJECT VERSION :\t', sobj@version), '\n')
  
  ### PROJECT NAME
  proj_name <- Seurat::Project(sobj)
  cat(paste0('PROJECT :\t[', proj_name, ']'), '\n')
  
  ### ASSAYS
  if('assay' %in% describe) {
    cat('\n[ASSAYS]\n')
    expassays <- Seurat::Assays(object = sobj)
    for (ea in seq_along(expassays)) {
      assay_name <- expassays[ea]
      cur_assay <- Seurat::Assays(object = sobj, slot = assay_name)
      cat(paste0('   ASSAY ', ea, ' :\t[', assay_name, ']', if(Seurat::DefaultAssay(sobj) == assay_name) ' [ACTIVE]'), '\n')
      ## SLOTS
      for (sl in seq_along(slots)) {
        slot_name <- slots[sl]
        cur_slot <- Seurat::GetAssayData(object = sobj, assay = assay_name, slot = slot_name)
        slot_dimprod <- prod(dim(cur_slot))
        cat(paste0('      SLOT/LAYER ', sl, ' :\t[', slot_name, ']\tDims:[', nrow(cur_slot), ' x ', ncol(cur_slot), if(slot_dimprod > 0) paste0(']  Range:[', paste(sprintf('%.2f', range(cur_slot, na.rm = TRUE)), collapse = '-')) else NULL, ']'), '\n')
        ## Total counts if slot is "counts"
        if(slot_name == 'counts') cat(paste0('         Counts :\t', round(sum(cur_slot))), '\n')
        ## Adding sparsity info when needed :
        if (sparse_level & is(cur_slot, 'dgCMatrix')) {
          splev <- sum(sparseMatrixStats::colCounts(x = cur_slot, value = 0)) / slot_dimprod
          cat(paste0('         Sparsity :\t', sprintf('%.5f', splev * 100), '%'), '\n')
        }
      }
    }
  }
  
  ### DIMRED
  if ('dimred' %in% describe) {
    cat('\n[DIMREDS]\n')
    dimreds <- Seurat::Reductions(object = sobj)
    for (dr in seq_along(dimreds)) {
      dr_name <- dimreds[dr]
      cur_dr <- Seurat::Reductions(object = sobj, slot = dr_name)
      cat(paste0('   DIMRED ', dr, ' : [', dimreds[dr], ']  Dims:[', nrow(cur_dr), ' x ', ncol(cur_dr), ']'), '\n')
    }
  }
  
  ### GRAPHS
  if(length(sobj@graphs) > 0) {
    cat('\n[GRAPHS]\n')
    for (gn in names(sobj@graphs)) {
      cat(paste0('  ', gn), '\n')
    }
  }
  
  ### BARCODE META
  if('coldata' %in% describe) {
    bc.df <- as.data.frame(sobj@meta.data)
    cat('\n[BARCODES METADATA]\n')
    for (b in seq_len(ncol(bc.df))) {
      b.fac <- as.factor(bc.df[,b])
      if (nlevels(b.fac) < max_levels) {
        # message('\n', colnames(bc.df)[b])
        b.tbl <- as.data.frame(table(b.fac, useNA = 'always'))
        colnames(b.tbl)[1] <- colnames(bc.df)[b]
        cat(knitr::kable(b.tbl, format = 'simple', escape = FALSE), sep = '\n')
      } else {
        cat('\n',colnames(bc.df)[b], '\n')
        print(summary(bc.df[,b]))
      }
    }
  }
}

# sobj2 <- readRDS('~/WORKSPACE/ENCADREMENT/2023/EBAII/TD_DATASETS/START_DATASETS/TDCT/TDCT_01a_EDf_Seurat.rds')
# sce <- readRDS('~/WORKSPACE/ENCADREMENT/2023/EBAII/TD_DATASETS/ALL_STEPS/TDCT/TDCT_01a_EDf.rds')


## CREATE A QUICK AND DIRTY 2D UMAP from a Seurat object
## TESTED ON : Seurat=4.4.0, sparseMatrixStats=1.13.4 
### sobj : a Seurat object
### assay : any assay name (default : 'RNA')
### slot : any slot ('counts', 'data', 'scale.data')
### reduction_method1 : any of ('pca', 'sca')
### reduction_method2 : any of ('umap', 'tsne')
### nfeatures : number of features (genes) to use at the different steps
### vtr: variables (names) to regress
### ncomp : number of first reduction dimension components used for the uMAP/tSNE reduction (note ; actually, 2 times more PCA components will be generated if slot mode is used)
### my_seed : seed for the RNG
### group_by : use a metadata factor to color the visualization
### features : either a continous metadata, or a feature name, to color the visualization
### dark_theme : use adrk background (better contrast)
QnD_viz <- function(sobj = NULL, assay = 'RNA', slot = 'counts', dimred = NULL, reduction_method1 = 'pca', reduction_method2 = 'umap', nfeatures = 2000, vtr = NULL, ncomp = 30, my_seed = 1337, group_by = NULL, features = NULL, pt_size = 1, dark_theme = TRUE, print_plot = TRUE, return_object = FALSE, return_plot = FALSE) {
  
  message('Checks ...')
  
  ### Mandatory
  if(is.null(sobj)) stop('A Seurat object is required !')
  if(!is.character(assay)) stop('Assay name should be a character (string)')
  if(ncomp <= 0) stop('[ncomp] should be a non-null positive integer (and <= N cells).')
  if(all(return_object, return_plot)) stop("Can't return Seurat object and plot object at the same time !")
  if(!any(return_object, return_plot, print_plot)) stop("No return of the Seurat object, nor the plot object, nor printing it : nothing to do !")
  ## Additional checks on sobj
  if(!is(sobj, 'Seurat')) stop('Provided sobj is not a proper Seurat object !')
  ## Handling slot/dimred modes
  if(!is.null(slot) & !is.null(dimred)) stop('Cannot choose between slot and dimred !')
  if(is.null(slot) & is.null(dimred)) {
    stop('Reduction has to be performed on a data matrix, or a dimension reduction, but both slot and dimred are set to NULL !')
  } else if (!is.null(slot)) {
    message('slot set : dimred will be ignored')
    valid.slots <- c('counts', 'data', 'scale.data')
    if(!slot %in% valid.slots) {
      stop('Unknown slot ! should be one of [', paste(valid.slots, collapse = ', '), ']')
    }
    ## Checking if requested assay exists
    cur_assays <- Seurat::Assays(sobj)
    if (!assay %in% cur_assays) stop('Requested assay [', assay, '] not found. Available assays are : [', paste(cur_assays, collapse = ', '), ']')
    
    ## Checking if slot exists
    cur_slots <- methods::slotNames(Seurat::Assays(object = sobj, slot = assay))
    ### Handling Seurat5
    if ('layers' %in% cur_slots) cur_slots <- base::names(Seurat::Assays(object = sobj, slot = assay)@layers)
    if(!slot %in% cur_slots) stop('Slot [', slot, '] not found in assay [', assay, '] !')
    
    ## Check slot population
    if (prod(dim(Seurat::GetAssayData(object = sobj, assay = assay, slot = slot))) == 0) stop('Requested slot [', slot, '] is an empty matrix !')
    
  } else {
    message('dimred set : slot will be ignored')
    
    ## Check if requested dimred exists
    if(!is.null(dimred)) {
      if (!dimred %in% Seurat::Reductions(sobj)) stop('Provided dimred not found !')
      # if (ncol(sobj@reductions[[dimred]]) < ncomp) stop('The provided dimred contains fewer dimensions than requested by ncomp !')
    }
  }
  
  ## Check group_by
  if (!is.null(group_by)) {
    if (!group_by %in% colnames(sobj@meta.data)) stop('group_by set but not found in sobj metadata !')
    if(is.numeric(sobj@meta.data[[group_by]])) stop('group_by should not be numeric !')
  }
  
  ## Check VTR
  if (!is.null(vtr)) {
    if (!all(vtr %in% colnames(sobj@meta.data))) stop('Not all vtr names found in Seurat object meta data !')
  }
  
  ## SLOT CASE
  if (!is.null(slot)) {
    ## Raw counts case : Normalize
    if (slot == 'counts') {
      message('NormalizeData')
      sobj <- Seurat::NormalizeData(object = sobj, assay = assay, verbose = FALSE)
      slot <- 'data'
    }
    ## Normalized, unscaled case : HVGs & scale
    if(slot == 'data') {
      message('FindVariableFeatures')
      sobj <- Seurat::FindVariableFeatures(object = sobj, assay = assay, nfeatures = nfeatures, verbose = FALSE)
      message('ScaleData')
      sobj <- Seurat::ScaleData(object = sobj, assay = assay, verbose = FALSE, vars.to.regress = vtr)
      slot <- 'scale.data'
    } 
    
    ## DIMRED 1
    if(reduction_method1 == 'pca') {
      ## PCA
      message('PCA')
      sobj <- Seurat::RunPCA(object = sobj, assay = assay, features = rownames(Seurat::GetAssayData(object = sobj, assay = assay, slot = slot)), npcs = ncomp*2, seed.use = my_seed, verbose = FALSE)
      dimred <- 'pca'
    } else stop('Unexpected reduction_method1 !')
  }
  
  if (!dimred %in% c('umap', 'tsne')) {
    ## Graph
    message('FindNeighbors')
    sobj <- Seurat::FindNeighbors(object = sobj, assay = assay, graph.name = paste0(assay, c('_nn','_snn')), reduction = dimred, dims = 1:ncomp, features = nfeatures, verbose = FALSE)
    
    ## DIMRED 2
    if(reduction_method2 == 'umap') {
      ### uMAP
      message('UMAP')
      suppressMessages(sobj <- Seurat::RunUMAP(object = sobj, assay = assay, reduction = dimred, seed.use = my_seed, dims = 1:ncomp, verbose = FALSE))
    } else if(reduction_method2 == 'tsne') {
      ### TSNE
      message('TSNE')
      suppressMessages(sobj <- Seurat::RunTSNE(object = sobj, assay = assay, graph.name = paste0(assay, '_snn'), reduction = dimred, seed.use = my_seed, dims = 1:ncomp, verbose = FALSE))
    }
    dimred <- reduction_method2
  }
  
  ## Plot
  
  ### Load patchwork &operator
  `&` <- patchwork:::`&.gg`
  ### Handle plot
  if (any(print_plot, return_plot)) {
    if(!is.null(group_by)) {
      pl <- Seurat::DimPlot(
        object = sobj, 
        dims = c(1,2), 
        reduction = dimred, 
        seed = my_seed, 
        group.by = group_by, 
        pt.size = pt_size, 
        shuffle = TRUE) & if(dark_theme) Seurat::DarkTheme()
      # if (dark_theme) print(Seurat::LabelClusters(plot = pl, id = group_by, color = 'white')) else print(Seurat::LabelClusters(plot = pl, id = group_by, color = 'black'))
      out_plot <- if (dark_theme) Seurat::LabelClusters(plot = pl, id = group_by, color = 'white') else Seurat::LabelClusters(plot = pl, id = group_by, color = 'black')
    } else if (!is.null(features)) {
      out_plot <- Seurat::FeaturePlot(
        object = sobj, 
        dims = c(1,2), 
        reduction = dimred, 
        features = features, 
        slot = 'data', 
        pt.size = pt_size, 
        order = TRUE,
        cols = c('white', 'blue')) & if(dark_theme) Seurat::DarkTheme()
    } else {
      out_plot <- Seurat::DimPlot(object = sobj, dims = c(1,2), reduction = dimred, seed = my_seed, pt.size = pt_size, shuffle = TRUE) & if(dark_theme) Seurat::DarkTheme()
    }
    ## Print_plot
    if(print_plot) print(out_plot)
    
    ## Return plot ?
    if(return_plot) return(out_plot)
  }
  
  ## Return object when requested
  if(return_object) return(sobj)
}

## Function to perform SoupX on a Seurat object
## TESTED on SoupX=1.6.2 with igraph=1.5.1, Seurat=4.4.0
## scmat_filt : (sparse)matrix corresponding to an empty droplets -filtered count matrix
## scmat_raw : (sparse)matrix corresponding to an NON- empty droplets -filtered count matrix
## soupQuantile, contaminationRange : see ?SoupX::autoEstCont
## contaminationRange : limits of the expected soup proportion
## soupRange : estimate soup fraction from features with total reads comprised in this range (only used when scmat_raw != NULL).
## return_object : TRUE = return the "unsouped" Seurat object ; FALSE = only perform soup estimation and return the rho estimation
## doPlot : Perform the SoupX estimation plot (rho distributions)
SoupX_auto <- function(scmat_filt = NULL, scmat_raw = NULL, soupQuantile = 0.9, contaminationRange = c(.01, .8), soupRange = c(0,100), return_object = FALSE, doPlot = FALSE) {
  
  ## Checks
  if(is.null(scmat_filt)) stop('A filtered count matrix is required !')
  
  if(is.null(scmat_raw)) message('No unfiltered raw counts matrix provided. Estimation will be based on filtered matrix only.')
  
  ## If no raw matrix
  if (is.null(scmat_raw)) {
    spChanRaw <- SoupX::SoupChannel(tod = scmat_filt, toc = scmat_filt, calcSoupProfile = FALSE)
    sc_rowsum <- sparseMatrixStats::rowSums2(scmat_filt)
    spProf <- data.frame(row.names = rownames(scmat_filt), est = sc_rowsum/sum(scmat_filt), counts = sc_rowsum)
    spChan <- SoupX::setSoupProfile(spChanRaw, spProf)
  } else {
    spChan <- SoupX::SoupChannel(tod = scmat_raw, toc = scmat_filt, calcSoupProfile = FALSE)
    if (min(spChan$nDropUMIs) > max(soupRange)) stop('Minimum found counts per barcode is : ', min(spChan$nDropUMIs), ', which is smaller than the upper bound of soupRange ! Please increase soupRange max !')
    spChan <- SoupX::estimateSoup(sc = spChan, soupRange = soupRange)
  }
  ## Display Top 20 contributing genes
  if (!return_object) {
    cat('\nSoup-contributing features (Top 20) :\n')
    print(knitr::kable(head(spChan$soupProfile[order(spChan$soupProfile$est, decreasing = TRUE), ], n = 20)))
  }
  ## Quick clustering needed
  spClust <- scran::quickCluster(scmat_filt, method = "igraph")
  ## Adding clusters to the SoupChannel object
  spChan <- SoupX::setClusters(sc = spChan, clusters = spClust)
  ## Estimating soup
  sX <- SoupX::autoEstCont(sc = spChan, doPlot = doPlot, tfidfMin = 1, 
                           soupQuantile = soupQuantile, maxMarkers = 100, 
                           contaminationRange = contaminationRange, rhoMaxFDR = .2, 
                           priorRho = .05, priorRhoStdDev = .1, 
                           forceAccept = FALSE)
  
  ## Removing soup (adjusting counts)
  if(return_object) {
    cat('Counts BEFORE SoupX : ', sum(scmat_filt), '\n')
    scmat_soupx <- SoupX::adjustCounts(sX, method = 'subtraction', roundToInt = TRUE, tol = .001, pCut = .01)
    cat('Counts AFTER SoupX : ', sum(scmat_soupx), '\n')
    rm(scmat_filt)
    return(scmat_soupx)
  } else return(sX$fit$rhoEst)
}

# ## Perform scBFA normalization/dimension reduction
# ## hvgs : 
# scbfa_norm <- function(sobj = NULL, assay = 'RNA', nfeatures = 2000, ndim = 50, vtr = NULL, vtr.scale = TRUE) {
#   
#   ### Mandatory
#   if (is.null(sobj)) stop('A Seurat object is required !')
#   if (!is(sobj, 'Seurat')) stop('Provided sobj is not a proper Seurat object !')
#   
#   ## CHecking HVGs
#   if (length(sobj@assays[[assay]]@var.features) == 0) sobj <- Seurat::FindVariableFeatures(sobj, assay = assay, nfeatures = nfeatures)
#   
#   ## Naming the reduction result
#   if (is.null(red.name)) red.name <- paste0(assay, '_', reduction.method)
#   
#   ## Handling variables to regress
#   if (!is.null(vtr)) {
#     if (!all(vtr %in% colnames(sobj@meta.data))) stop('Not all vtr names found in Seurat object meta data !')
#     minimeta <- sobj@meta.data[,colnames(sobj@meta.data) %in% vtr, drop = FALSE]
#     X <- matrix(ncol = 0, nrow = nrow(sobj@meta.data))
#     for(v in vtr) {
#       if(is.character(minimeta[[v]])) minimeta[[v]] <- as.factor(minimeta[[v]])
#       if(any(is.na(minimeta[[v]]))) stop(paste0("Covariate '", v, "' contains NA value(s) !"))
#       if(is.factor(minimeta[[v]])) {
#         message(paste0("Converting '", v, "' factor into model matrix and adding to the regression..."))
#         mm <- model.matrix(~minimeta[[v]])[,-1, drop = FALSE]
#         X <- cbind(X, mm)
#       } else {
#         message(paste0("Adding '", v, "' covariate", if(vtr.scale) ' (scaled)' else NULL,  ' to the regression ...'))
#         X <- cbind(X, if(vtr.scale) scale(minimeta[[v]]) else minimeta[[v]])
#       }
#     }
#   } else {
#     X <- NULL
#   }
#   
#   bfa.res <- scBFA::scBFA(scData = as.matrix(sobj@assays[[assay]]@counts[sobj@assays[[assay]]@var.features,]), numFactors = ndim, X = X)
#   dimnames(bfa.res$ZZ) <- list(colnames(sobj@assays[[assay]]@counts), paste0('SCBFA_', 1L:max.dims))
#   dimnames(bfa.res$AA) <- list(sobj@assays[[assay]]@var.features, paste0('SCBFA_', 1L:max.dims))
#   sobj@reductions[[red.name]] <- Seurat::CreateDimReducObject(embeddings = bfa.res$ZZ, loadings = bfa.res$AA, assay = assay, stdev = matrixStats::colSds(bfa.res$ZZ), key = paste0(red.name, '_'), misc = list())
#   return(sobj)
# }

## Computes Seurat's cell cycle phase/scores prediction
CC_Seurat <- function(sobj = NULL, assay = "RNA", seurat_cc_genes = NULL, SmG2M = TRUE, nbin = 24, my_seed = 1337L) {
  
  ## Check sobj
  if(is.null(sobj)) stop('A Seurat object is required !')
  if(!is(sobj, 'Seurat')) stop('Provided sobj is not a proper Seurat object !')
  cur_assays <- Seurat::Assays(sobj)
  if (!assay %in% cur_assays) stop('Requested assay [', assay, '] not found. Available assays are : [', paste(cur_assays, collapse = ', '), ']')
  
  ## Check if data slot/layer is populated
  if (crossprod(dim(Seurat::GetAssayData(object = sobj, assay = assay, slot = 'data')))[1] == 0) {
    sobj <- Seurat::NormalizeData(object = sobj, verbose = FALSE)
  }
  
  ## Run CC
  set.seed(my_seed)
  sobj <- Seurat::CellCycleScoring(object = sobj, s.features = seurat_cc_genes$s.genes, g2m.features = seurat_cc_genes$g2m.genes, assay = assay, nbin = nbin, seed = my_seed)
  sobj$CC_Seurat_S.Score <- sobj$S.Score
  sobj$CC_Seurat_G2M.Score <- sobj$G2M.Score
  sobj$CC_Seurat_Phase <- as.factor(sobj$Phase)
  if (SmG2M) sobj$CC_Seurat_SmG2M.Score <- sobj$S.Score - sobj$G2M.Score
  sobj$S.Score <- sobj$G2M.Score <- sobj$Phase <- NULL
  return(sobj)
}

## Computes scran's Cyclone cell cycle phase/scores prediction
CC_Cyclone <- function(sobj = NULL, assay = 'RNA', cyclone_cc_pairs = NULL, SmG2M = TRUE, my_seed = 1337) {
  
  ## Check sobj
  if(is.null(sobj)) stop('A Seurat object is required !')
  if(!is(sobj, 'Seurat')) stop('Provided sobj is not a proper Seurat object !')
  cur_assays <- Seurat::Assays(sobj)
  if (!assay %in% cur_assays) stop('Requested assay [', assay, '] not found. Available assays are : [', paste(cur_assays, collapse = ', '), ']')
  
  ## Run CC
  set.seed(my_seed)
  cycres <- scran::cyclone(x = Seurat::GetAssayData(object = sobj, assay = assay, slot = 'counts'), pairs = cyclone_cc_pairs, verbose = TRUE)
  sobj$CC_Cyclone_nG1.Score <- cycres$normalized.scores$G1
  sobj$CC_Cyclone_nS.Score <- cycres$normalized.scores$S
  sobj$CC_Cyclone_nG2M.Score <- cycres$normalized.scores$G2M
  sobj$CC_Cyclone_Phase <- as.factor(cycres$phases)
  if(SmG2M) sobj$CC_Cyclone_nSmG2M.Score <- cycres$normalized.scores$S - cycres$normalized.scores$G2M
  return(sobj)
}


assess_covar <- function(mat = NULL, covar.df = NULL, markers = NULL, ndim = 10L, center = TRUE, scale = TRUE, return.data = FALSE) {
  ## Checks
  if(is.null(mat)) stop('A matrix is required !')
  if(is.null(covar.df)) stop('A data.frame with covariates is required !')
  if(!is.data.frame(covar.df)) stop('covar.df should be a data.frame !')
  if (ncol(mat) != nrow(covar.df)) stop("ncol(mat) != nrow(covar.df) !")
  if(is.null(ndim)) stop('A number of dimensions to use is required !')
  if(!is.integer(ndim)) stop('ndim must be an integer > 0 !')
  if(ndim > ncol(mat)) stop('ndim must be <= ndim(mat) !')
  
  ## Check markers
  if(!is.null(markers)) markers <- markers[markers %in% rownames(mat)]
  
  ## Split factors and continuous data
  conti.names <- factor.names <- c()
  for (x in colnames(covar.df)) {
    if(is.numeric(covar.df[1,x])) conti.names <- c(conti.names, x) else factor.names <- c(factor.names, x)
  }
  col.names <- c(factor.names, conti.names, markers)
  col.types <- c(rep('factor', length(factor.names)), rep('continuous', length(conti.names)), rep('marker', length(markers)))
  
  ## Clean factors
  if (length(factor.names) > 0) {
    for (ccn in factor.names) {
      covar.df[[ccn]] <- as.factor(covar.df[[ccn]])
      covar.df[[ccn]] <- droplevels(covar.df[[ccn]])
    }
  }
  
  ## Remove factors with single level
  # covar.df[,vapply(seq_len(ncol(covar.df)), nlevels)]
  fac_sing <- vapply(covar.df, function(x) { if (!is.factor(x)) FALSE else if (nlevels(x) > 1) FALSE else TRUE }, TRUE)
  if (any(fac_sing)) {
    covar.df <- covar.df[,!fac_sing]
    factor.names <- factor.names[!factor.names %in% colnames(covar.df)[fac_sing]]
  }
  
  ## Convert conti to Zscores
  if (length(conti.names) > 0) {
    for (x in conti.names) {
      covar.df[[x]] <- (covar.df[[x]] - mean(covar.df[[x]], na.rm = TRUE)) / sd(covar.df[[x]], na.rm = TRUE)
    }
  }
  
  ## Convert markers to Zscores
  if (length(markers) > 0) {
    covar.df <- cbind(covar.df, t(mat[markers,]))
    for (x in markers) {
      covar.df[[x]] <- (covar.df[[x]] - mean(covar.df[[x]], na.rm = TRUE)) / sd(covar.df[[x]], na.rm = TRUE)
    }
  }
  
  ## Center / scale the matrix ?
  if (any(c(center, scale))) mat <- base::scale(x = mat, center = center, scale = scale)
  ## Dimension reduction
  norm.red <- irlba::irlba(A = t(mat), nv = min(ndim+1, ncol(mat)))$u
  
  rm(mat)
  
  ## Setting output matrix
  bc.mat <- matrix(NA, ncol = length(col.names), nrow = ndim, dimnames = list(paste0("PCA", seq_len(ndim)), col.names))
  ## Filling matrix
  for (cn in seq_along(col.names)) {
    # message(col.names[cn])
    if (col.names[cn] %in% c(conti.names, markers)) {
      cv2cor <- covar.df[[col.names[cn]]]
      nona <- !is.na(cv2cor)
      bc.mat[, cn] <-  abs(cor(x = cv2cor[nona], y = norm.red[nona,seq_len(ndim)], method = 'spearman'))
    } else if (col.names[cn] %in% factor.names) {
      b2kw <- covar.df[[col.names[cn]]]
      nona <- !is.na(b2kw)
      for (si in seq_len(ndim)) {
        k_test <- try(k_res <- kruskal.test(x = norm.red[nona,si], g = as.factor(b2kw[nona])), silent = TRUE)
        if (!is(k_test, class2 = 'try-error')) {
          bc.mat[si,cn] <- k_res$statistic / nrow(norm.red)
        } else bc.mat[si,cn] <- 0
      }
    }
  }
  ## Heatmap
  color.palette = c("white", "orangered3")
  myRamp.col <- circlize::colorRamp2(c(0, 1), color.palette)
  BC.hm <- ComplexHeatmap::Heatmap(matrix = bc.mat,
                                   name = 'Weight',
                                   col = myRamp.col,
                                   na_col = 'grey75',
                                   cluster_rows = FALSE,
                                   cluster_columns = FALSE,
                                   rect_gp = grid::gpar(col = "darkgrey", lwd=0.5),
                                   column_title = 'Batch factors and covariates weight on dataset',
                                   row_title = 'PCA dimensions',
                                   column_split = col.types,
                                   top_annotation = ComplexHeatmap::HeatmapAnnotation(Type = col.types, col = list(Type = setNames(object = c('lightblue','pink', 'green3'), nm = c('factor', 'continuous', 'marker')))))
  print(BC.hm)
  
  ## Return the weight matrix ?
  if(return.data) return(bc.mat)
}


## Cell annotation using SINGLER (by default with ImmGenData db from celldex)
cell_annot <- function(sobj = NULL, assay = 'RNA', slot = 'data', celldex_setname = 'ImmGenData', group_by = NULL, list_celldex = FALSE) {
  
  ## Get list of celldex databases
  library(celldex)
  cdx <- ls('package:celldex')
  detach(name = 'package:celldex')
  
  ## Return the celldex list
  if(list_celldex) {
    message(paste0(cdx, collapse = '\n'))
  } else {
    
    ## Check sobj
    if(is.null(sobj)) stop('A Seurat object is required !')
    if(!is(sobj, 'Seurat')) stop('Provided sobj is not a proper Seurat object !')
    ## Checking if assay exists
    cur_assays <- Seurat::Assays(sobj)
    if (!assay %in% cur_assays) stop('Requested assay [', assay, '] not found. Available assays are : [', paste(cur_assays, collapse = ', '), ']')
    ## Checking if slot exists
    cur_slots <- methods::slotNames(Seurat::Assays(object = sobj, slot = assay))
    if(!slot %in% cur_slots) stop('Slot [', slot, '] not found in assay [', assay, '] !')
    ## Check requested celldex db
    if (!celldex_setname %in% cdx) stop('Unsupported celldex database !')
    
    ## Check requested assay
    okassays <- Seurat::Assays(object = sobj)
    if(!assay %in% okassays) stop('Assay [', assay, '] not found !')
    
    ## Check slot
    okslot <- c('data', 'scale.data')
    if(!slot %in% okslot) stop('Unsupported slot !')
    
    ## Check group_by
    if (!is.null(group_by)) {
      if (!group_by %in% colnames(sobj@meta.data)) stop('group_by [', group_by, '] not found in sobj metadata !')
    }
    
    ## Load annotation
    annotation <- celldex::ImmGenData()
    print(celldex_setname)
    library(celldex)
    annotation <- suppressMessages(do.call(match.fun(paste0(celldex_setname)), args = list()))
    detach("package:celldex")
    
    ## Get matrix 
    norm_exp_mat <- Seurat::GetAssayData(object = sobj, assay = assay, slot = slot)
    
    ## Handle group_by
    if(!is.null(group_by)) group_by <- droplevels(as.factor(sobj[[group_by, drop = TRUE]]))
    
    ## Predict
    ann_predictions <- SingleR::SingleR(
      test = norm_exp_mat
      , ref = annotation
      , clusters = group_by
      , labels = annotation$label.fine
      , assay.type.test = "logcounts"
      , assay.type.ref = "logcounts"
      , BPPARAM = BiocParallel::SerialParam()
    )
    if (is.null(group_by)) {
      return(as.factor(ann_predictions$labels))
    } else {
      levels(group_by) <- ann_predictions$labels
      return(group_by)
    }
  }
}
