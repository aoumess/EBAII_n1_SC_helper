## Seurat-descriptor
## TESTED ON : Seurat=4.4.0, sparseMatrixStats=1.13.4 , knitr=.44
## sobj           [Seurat]    A Seurat object
## max_levels     [int>0]     Maximal number of unique values to consider a barcode annotation as a factor rather than a continuous numeric vector
## describe       ['all', 'assays', 'dimred', 'coldata']    Type of entries to describe
seurat4_descriptor <- function(sobj = NULL, describe = 'all', sparse_level = TRUE, max_levels = 100) {
  
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
  
  # ## Handling merged case
  # if(length(sobj@meta.data$misc) > 1 & all(vapply(sobj@metadata$misc, is.list, TRUE)))
  #   sobj@metadata$misc <- list(projectname = paste(vapply(sobj@metadata$misc, function(x) x$projectname, "a"), collapse = '.'), id = max(vapply(sobj@metadata$misc, function(x) x$id, 1)))
  
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
        cat(paste0('      SLOT ', sl, ' :\t[', slot_name, ']\tDims:[', nrow(cur_slot), ' x ', ncol(cur_slot), if(slot_dimprod > 0) paste0(']  Range:[', paste(sprintf('%.2f', range(cur_slot, na.rm = TRUE)), collapse = '-')) else NULL, ']'), '\n')
        ## Adding sparsity info when needed :
        if (sparse_level & is(cur_slot, 'dgCMatrix')) {
          splev <- sum(sparseMatrixStats::colCounts(x = cur_slot, value = 0)) / slot_dimprod
          cat(paste0('         Sparsity :\t', sprintf('%.5f', splev * 100), '%'), '\n')
          if(slot_name == 'counts') cat(paste0('         Counts :\t', round(sum(cur_slot))), '\n')
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
QnD_viz <- function(sobj = NULL, assay = 'RNA', slot = 'counts', dimred = NULL, reduction_method1 = 'pca', reduction_method2 = 'umap', nfeatures = 2000, vtr = NULL, ncomp = 30, my_seed = 1337, group_by = NULL, features = NULL, pt_size = 1, dark_theme = TRUE, return_object = FALSE) {
  
  message('Checks ...')
  
  ### Mandatory
  if (is.null(sobj)) stop('A Seurat object is required !')
  if(!is.character(assay)) stop('Assay name should be a character (string)')
  if(ncomp <= 0) stop('[ncomp] should be a non-null positive integer (and <= N cells).')
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
  if(!is.null(group_by)) {
    pl <- Seurat::DimPlot(object = sobj, dims = c(1,2), reduction = dimred, seed = my_seed, group.by = group_by, pt.size = pt_size, shuffle = TRUE) + if(dark_theme) Seurat::DarkTheme()
    if (dark_theme) print(Seurat::LabelClusters(plot = pl, id = group_by, color = 'white')) else print(Seurat::LabelClusters(plot = pl, id = group_by, color = 'black'))
  } else if (!is.null(features)) {
    print(Seurat::FeaturePlot(object = sobj, dims = c(1,2), reduction = dimred, features = features, slot = 'data', pt.size = pt_size) + if(dark_theme) Seurat::DarkTheme())
  } else {
    print(Seurat::DimPlot(object = sobj, dims = c(1,2), reduction = dimred, seed = my_seed, pt.size = pt_size, shuffle = TRUE) + if(dark_theme) Seurat::DarkTheme())
  }
  
  ## Return object when requested
  if(return_object) return(sobj)
}

## Function to perform SoupX on a Seurat object
## TESTED on SoupX=1.6.2 with igraph=1.5.1, Seurat=4.4.0
## sobj : (Seurat object corresponding to an empty droplets -filtered dataset
## scmat_raw : (sparse)matrix corresponding to an NON- empty droplets -filtered dataset
## soupQuantile, contaminationRange : see ?SoupX::autoEstCont
## contaminationRange : limits of the expected soup proportion
## soupRange : estimate soup fraction from features with total reads comprised in this range (only used when scmat_raw != NULL).
## return_object : TRUE = return the "unsouped" Seurat object ; FALSE = only perform soup estimation and return the rho estimation
## doPlot : Perform the SoupX estimation plot (rho distributions)
SoupX_auto <- function(sobj = NULL, assay = 'RNA', scmat_raw = NULL, soupQuantile = 0.9, contaminationRange = c(.01, .8), soupRange = c(0,100), return_object = FALSE, doPlot = FALSE) {
  
  ## Checks
  if(is.null(sobj)) stop('A Seurat object is required !')
  if(!is(sobj, 'Seurat')) stop('Provided sobj is not a proper Seurat object !')
  if(is.null(scmat_raw)) message('No unfiltered raw counts matrix provided. Estimation will be based on filtered matrix only.')
  cur_assays <- Seurat::Assays(sobj)
  if (!assay %in% cur_assays) stop('Requested assay [', assay, '] not found. Available assays are : [', paste(cur_assays, collapse = ', '), ']')
  
  ## Get expression matrix
  scmat_filt <- Seurat::GetAssayData(object = sobj, assay = assay, slot = 'counts')
  
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
  message('Soup-contributing features (Top 20) :')
  print(knitr::kable(head(spChan$soupProfile[order(spChan$soupProfile$est, decreasing = TRUE), ], n = 20)))
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
  ## Display soup markers
  # message('Soup-marking features (Top 50) :')
  # print(knitr::kable(sX$fit$markersUsed[1:50, c('gene', 'geneFrequencyGlobal', 'qval')]))
  
  ## Removing soup (adjusting counts)
  if(return_object) {
    message('Counts BEFORE SoupX : ', sum(scmat_filt))
    scmat_soupx <- SoupX::adjustCounts(sX, method = 'subtraction', roundToInt = TRUE, tol = .001, pCut = .01)
    message('Counts AFTER SoupX : ', sum(scmat_soupx))
    # sobj@assays[[assay]]@counts <- scmat_soupx
    newmeta <- as.data.frame(sobj@meta.data)
    newmeta <- newmeta[,!colnames(newmeta) %in% c('orig.ident', paste0('nCount_', assay),  paste0('nFeature_', assay))]
    sobj <- Seurat::CreateSeuratObject(counts = scmat_soupx,project = sobj@project.name, meta.data = newmeta)
    rm(scmat_soupx, scmat_filt)
    return(sobj)
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

## Add a barcode-level metric based on counts (mito_rate, ) 
count_rate_metric <- function(sobj = NULL, assay = 'RNA', features = NULL) {
  ## Check sobj
  if(is.null(sobj)) stop('A Seurat object is required !')
  if(!is(sobj, 'Seurat')) stop('Provided sobj is not a proper Seurat object !')
  cur_assays <- Seurat::Assays(sobj)
  if (!assay %in% cur_assays) stop('Requested assay [', assay, '] not found. Available assays are : [', paste(cur_assays, collapse = ', '), ']')
  
  ## Get raw count expression matrix
  count_mat <- Seurat::GetAssayData(object = sobj, assay = assay, slot = 'counts')
  infeats <- rownames(count_mat) %in% features
  metric_rate <- unname(as.vector(Matrix::colSums(count_mat[infeats,]) / sobj$nCount_RNA))
  rm(count_mat)
  return(metric_rate)
}

## Computes Seurat's cell cycle phase/scores prediction
CC_Seurat <- function(sobj = NULL, assay = 'RNA', seurat_cc_genes = NULL, SmG2M = TRUE, nbin = 24, my_seed = 1337) {
  
  ## Check sobj
  if(is.null(sobj)) stop('A Seurat object is required !')
  if(!is(sobj, 'Seurat')) stop('Provided sobj is not a proper Seurat object !')
  cur_assays <- Seurat::Assays(sobj)
  if (!assay %in% cur_assays) stop('Requested assay [', assay, '] not found. Available assays are : [', paste(cur_assays, collapse = ', '), ']')
  
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
  cycres <- scran::cyclone(Seurat::as.SingleCellExperiment(sobj, assay = assay), pairs = cyclone_cc_pairs, verbose = TRUE)
  sobj$CC_Cyclone_nG1.Score <- cycres$normalized.scores$G1
  sobj$CC_Cyclone_nS.Score <- cycres$normalized.scores$S
  sobj$CC_Cyclone_nG2M.Score <- cycres$normalized.scores$G2M
  sobj$CC_Cyclone_Phase <- as.factor(cycres$phases)
  if(SmG2M) sobj$CC_Cyclone_nSmG2M.Score <- cycres$normalized.scores$S - cycres$normalized.scores$G2M
  return(sobj)
}

## Performs a "correlation analysis" of covariates with dimred components.
## . Spearman correlation for continuous covariates
## . KW for factors
## . Returns the matrix of comp x covariate p-values (off by default)
## . Plots a heatmap of -log10(p)
## sobj : Seurat object
## dimred : name of a reduction in the Seurat object
## return_p : if TRUE, returns the p-values matrix
dimred_covar_cor <- function(sobj = NULL, dimred = 'pca', ncomp = 10, return_p = FALSE) {
  
  ## Check sobj
  if(is.null(sobj)) stop('A Seurat object is required !')
  if(!is(sobj, 'Seurat')) stop('Provided sobj is not a proper Seurat object !')
  if (! 'pca' %in% Seurat::Reductions(sobj)) stop('Reduction [', dimred, '] not found in sobj !')
  ## Get reduction
  cur_red <- Seurat::Reductions(object = sobj, slot = dimred)@cell.embeddings
  if (ncomp > ncol(cur_red)) {
    ncomp <- ncol(cur_red)
    message('More dimensions requested than available : Limiting to available ones (', ncomp, ') ...')
  }
  cur_meta <- as.data.frame(sobj@meta.data)
  ## Convert logical & character to factor
  meta_types <- vapply(seq_along(colnames(cur_meta)), function(x) is(cur_meta[,x, drop = TRUE])[1], 'a')
  for (z in which(meta_types == 'logical' | meta_types == 'character')) cur_meta[,z] <- as.factor(cur_meta[,z])
  ## Drop missing levels in factors
  meta_types <- vapply(seq_along(colnames(cur_meta)), function(x) is(cur_meta[,x, drop = TRUE])[1], 'a')
  for (z in which(meta_types == 'factor')) cur_meta[,z] <- droplevels(cur_meta[,z])
  
  pCor <- pcaExplorer::correlatePCs(pcaobj = list(x =  cur_red), coldata = cur_meta, pcs = 1:ncomp)
  ## FDR adjustment
  pCor <- matrix(data = p.adjust(p = pCor, method = 'BH'), nrow = nrow(pCor), ncol = ncol(pCor), dimnames = dimnames(pCor))
  ## Heatmap
  pl <- -log10(pCor)
  pl[is.infinite(pl)] <- 999
  rownames(pl) <- colnames(cur_red)[1:ncomp]
  print(ComplexHeatmap::Heatmap(matrix = pl, name = Seurat::Project(object = sobj), cluster_rows = FALSE, cluster_columns = FALSE, col = circlize::colorRamp2(c(0,1000), c('grey95', 'blue'))))
  ## Return p-values
  if(return_p) return(pCor)
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
