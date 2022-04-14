#' "Ancestry matching" procedure to filter residual germline variants from cell line 96 tri-nucleotide contexts. 
#' 
#' @param nClusters The number of clusters for determining the "ancestry" clusters (default 13)
#' @param nCPUs The number of CPUs for parallel computing
#' @param output.file The name of the output file where 96 contexts will be written
#' @param PCA.components Data frame containing PCA components from the PCA analysis of common variants of TCGA samples and cell lines (columns = PC components, rows = samples, row.names = IDs of samples)
#' @param contexts.96.CCLs Data frame containing 96 trinucleotide contexts of cancer cell lines (CCLs) (columns = 96 contexts, rows = samples, row.names = IDs of samples). No need to provide if remove.CNVkit.deletions is used.
#' @param contexts.96.TCGA.Germline Data frame containing 96 trinucleotide contexts of TCGA germline samples (columns = 96 contexts, rows = samples, row.names = IDs of samples). No need to provide if remove.CNVkit.deletions is used.
#' @param remove.CNVkit.deletions Set to TRUE to filter variants from TCGA germline samples that are in the deleted regions of cell lines, as determined by copy number segment determined by the CNVkit software  (default FALSE). If this option is used cell line and TCGA germline variants have to be provided.
#' @param CNVkit.folder Folder with copy nubmer segments (.cns files) from CNVkit. File names have to correspond to cell line IDs (<ID>.cns)
#' @param CNVkit.log2.threshold Threshold for filtering, variants in regions with cns log2 < threshold are filtered (default -0.3).
#' @param variants.CCLs List with cell line variants. Names of list elements correspond to cell line IDs, while each element contains data frame with variants of that cell line (columns: Chr, Start, Ref, Alt, End). Chr is in form: 1, 2, 3, ..., X, Y.
#' @param variants.TCGA.Germline List with cell line variants
#'
#' @return Data frame with 96 tri-nucleotide contexts of cell lines with determined by the ancestry matching procedure
#' @export
#'
#' @examples
ancestryMatching <- function(nClusters = 13, nCPUs = 1, output.file, PCA.components, contexts.96.CCLs, contexts.96.TCGA.Germline = NULL, remove.CNVkit.deletions = FALSE, CNVkit.folder = "", CNVkit.log2.threshold = -0.3, variants.CCLs = NULL, variants.TCGA.Germline = NULL) {
  
  require(tclust)
  require(pcaPP)
  require(foreach)
  require(doParallel)
  
  #make sure that cell line and TCGA 96 contexts are equally ordered
  if(sum(colnames(contexts.96.CCLs) != colnames(contexts.96.TCGA.Germline)) > 0) {
    stop("The columns of cell line and TCGA 96 contexts need to have equal order")
  }
  
  #cluster cell line and tumor samples on the basis of PCA components
  tclus <- tclust(PCA.components, k = nClusters, alpha = 0.05, iter.max = 100, restr = "sigma")
  
  PCA.components$cluster = tclus$cluster
  PCA.components <- PCA.components[PCA.components$cluster > 0,]  #remove outliers as determined by tclust
  
  isCellLine = row.names(PCA.components) %in% row.names(contexts.96.CCLs)
  isTCGA = row.names(PCA.components) %in% row.names(contexts.96.TCGA.Germline)
  
  CCL.clusters = data.frame(ID = row.names(PCA.components)[isCellLine], cluster = PCA.components$cluster[isCellLine]) #cell lines IDs and their corresponding clusters, some cell lines may have been removed as tclust ouliers 
  TCGA.clusters = data.frame(ID = row.names(PCA.components)[isTCGA], cluster = PCA.components$cluster[isTCGA]) #TCGA IDs and their corresponding clusters, some samples may have been removed as tclust ouliers 
  
  #calculate median "ancestors" of each TCGA cluster. Note if remove.CNVkit.deletions is used, 96 contexts have to be calculated again 
  if(!remove.CNVkit.deletions) {
    medianAncestors <- list()
    for(c in unique(TCGA.clusters$cluster)) {
      medianAncestors[[c]] <- round(l1median(X = contexts.96.TCGA.Germline[row.names(contexts.96.TCGA.Germline) %in% TCGA.clusters$ID[TCGA.clusters$cluster == c],]))
      names(medianAncestors[[c]]) = colnames(contexts.96.TCGA.Germline)
    } 
  }
  
  cl <- makeCluster(nCPUs)
  registerDoParallel(cl)
  clusterExport(cl = cl, list("CCL.clusters", "medianAncestors", "contexts.96.CCLs"), envir=environment())
  if (remove.CNVkit.deletions) {
    source("Code/getMutationContext.R")
    clusterExport(cl = cl, list("variants.TCGA.Germline"))
  }
  
  contexts.96.CCLs.ancestryMatching <- foreach(i = 1:nrow(CCL.clusters), .combine = rbind) %dopar% {
    
    #cell line's cluster
    cluster = CCL.clusters$cluster[i]
    
    #96 context of a cell line
    if(!remove.CNVkit.deletions) {
      CCL.96.context <- contexts.96.CCLs[row.names(contexts.96.CCLs) == CCL.clusters$ID[i],]
    }
    
    #remove variants from TCGA germline samples that are in the deleted regions of cell lines (as determined by CNVkit). 
    #Note, for this 96 contexts of TCGA germline samples have to be recalculated (for each cell line), therefore this can be time consuming
    if (remove.CNVkit.deletions) {
      
      cns.mask <- read.csv(paste0(CNVkit.folder, CCL.clusters$ID[i], ".cns"))
      
      #cell line
      tmpVars <- variants.CCLs[[CCL.clusters$ID[i]]]
      tmpVars <- filterVariantsCNVkit(tmpVars, cns.mask)
      
      CCL.96.context = get_96_context(tmpVars)
      CCL.96.context = CCL.96.context[,order(colnames(CCL.96.context))] #sort to ensure equal ordering
      
      #TCGA sampleds
      samples = TCGA.clusters$ID[TCGA.clusters$cluster == cluster]
      
      medianAncestors <- list()
      tmpContexts <- list()
      
      for (s in samples) {
        tmpVars <- variants.TCGA.Germline[[s]]
        tmpVars <- filterVariantsCNVkit(tmpVars, cns.mask)
        tmpContexts[[s]] = get_96_context(tmpVars)
      }
      
      tmpContexts <- data.frame(matrix(unlist(tmpContexts), nrow = length(tmpContexts), byrow = T))
      colnames(tmpContexts) <- names(contextsGermlineAll[[tmp[1]]])
      
      tmpContexts = tmpContexts[,order(colnames(tmpContexts))] #sort to ensure equal ordering
      
      medianAncestors[[cluster]] <- round(l1median(X = tmpContexts))
    }
    
    CCL.96.context.subtracted = CCL.96.context - medianAncestors[[cluster]] #subtract ancestor's 96 profile
    CCL.96.context.subtracted[CCL.96.context.subtracted < 0] <- 0 #clamp negative values to 0
    
    CCL.96.context.subtracted
  }
  
  stopCluster(cl)
  
  write.csv(contexts.96.CCLs.ancestryMatching, outFile)
}


#' Filter variants that are in the specified copy number seqments (cns), as given by the CNVkit software
#'
#' @param variants Data frame containing SNVs, first column must be Chr and second column position
#' @param cns.mask Data frame with copy number seqments as given by the CNVkit software
#' @param expandBy if set, each segment is expanded by the specified number of bases (upstream and downstream) 
#' @param cns.mask.threshold The threshold for the log2 column of cns.mask used for filtering of variants (default is -0.3; i.e., variants in regions with log2 < -0.3 are filtered out) 
#'
#' @return Data frame with filtered variants  
filterVariantsCNVkit <- function(variants, cns.mask, expandBy = 0, cns.mask.threshold = -0.3) {
  return(variants[!apply(variants, 1, cns, cns.mask, expandBy, cns.mask.threshold),])
}

#' Helper function for filterVariantsCNVkit
cns <- function(x, cns.mask, expandBy = 0, cns.mask.threshold = -0.3) {
  sum(cns.mask$log2 < cns.mask.threshold & x[1] == cns.mask$chromosome & (cns.mask$start - expandBy) <= as.numeric(x[2]) & as.numeric(x[2]) <= (cns.mask$end + expandBy)) > 0
}

