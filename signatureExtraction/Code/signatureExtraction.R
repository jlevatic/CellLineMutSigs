library(dplyr)
library(nmfgpu4R)
library(data.table)
library(stringr)
nmfgpu4R.init()


#' The signature extraction based on non-negative matrix factorization. The code relies on GPU NMF, therefore, to use this Cuda have to be installed (with appropriate GPU that supports it).
#' Signatures are extracted, then compared to a set of reference signatures (i.e., PCAWG), to select the ones that best resemble reference signatures. 
#'
#' @param contexts Data frame containing contexts for signature extraction (96 tri-nucleotide contexts and optionally indel features). Rows = samples, Cols = features, row.names = sample IDs.
#' @param output.name The name of the output file.
#' @param nBootstrapReplicates The number of bootstrap replicates (default 300).
#' @param maxK Maximal number of signatures to consider, it will go from minK...maxK (default 40).
#' @param minK Minimal number of signatures to consider (default 2).
#' @param threshold.hier Cosine similarity threshold for hierarchical extraction of signatures. Samples that are reconstructed with cosine similarity >= threshold are removed and signature extraction is repeated (default 0.97). 
#' @param max.iter.hier The maximal number of iteration for hierarchical extraction to be performed (default 3). Note that, for the value 1, hierarchical extraction will not be performed.
#' @param tissueAware If TRUE, the best matching signatures will we extracted by also considering the similarity on mean signature exposures across cancer types
#' @param tissueAWare.weight The weight by which the context similarity and the exposure similarity are combined to determine the best matching signatures ( weight * contextSimilarity + (1-weight) * exposureSimilarity )
#' @param referenceSignatures Data frame with reference signatures (rows = signatures, columns = contexts) used to find the best matching signatures
#' @param referenceSignatures.exposureMeans Data frame with means per cancer type of exposures of reference signatures (rows = cancer types, columns = signatures)
#' @param similarityThreshold The minimal cosine similairty on context required to match the extracted signatures to reference signatures (default 0.85)
#' @param save.data If TRUE the all of the resuling signatures/exposures (considering different bootstrap replicates and different number of signatures) will be saved to files. 
#'
#' @return Tables with extracted signatures their corresponding exposure scores.
#' @export
#'
#' @examples
extractSignatures = function(contexts, output.name, nBootstrapReplicates = 300, maxK = 40, minK = 2, threshold.hier = 0.97, max.iter.hier = 3, referenceSignatures, referenceSignatures.exposureMeans, tissueAware = TRUE, tissueAWare.weight = 0.9, similarityThreshold = 0.85, save.data = FALSE) {
  
  cellLine.cancerTypes = data.frame(COSMIC_ID = row.names(contexts), cancer_type = contexts$cancer_type)
  
  triSig.dataM = as.matrix(contexts[,-1])
  
  ##
  ## Parameters and initializatin of some objects
  ##
  triSig.nmfHmatAllByFact = list()
  triSig.nmfWmatAllByFact = list()
  triSig.nmfFrobeniusError = list()
  triSig.nmfRMSD = list()
  triSig.dataMResamp = list() #prepared the resampled NMF input matrices, generated just once and then re-used for any possible # factors/clusters
  
  ##
  ## Generate matrices with small random pertubations of the original matrix 
  ##
  set.seed(42)
  triSig.dataMResamp = list()
  for (nIter in 1:nBootstrapReplicates) {
    cat(sprintf("Generating bootstrap matrix: nIter %d\n", nIter))
    triSig.dataMTempIter = matrix(data = NA,
                                  nrow = nrow(triSig.dataM), ncol = ncol(triSig.dataM),
                                  dimnames = list(rownames(triSig.dataM), colnames(triSig.dataM))
    )
    
    for (aRow in 1:nrow(triSig.dataM)) {
      # for each sample (row in input matrix)
      triSig.dataMTempIter[aRow,] = sampling::UPmultinomial(triSig.dataM[aRow,]) #this changes the row a bit while keeping the total sum of mutations the same
    }
    
    attr(triSig.dataMTempIter, "nIter") = nIter
    triSig.dataMResamp[[nIter]] = triSig.dataMTempIter;
  }
  rm(triSig.dataMTempIter, nIter)
  
  ##
  ## Run NMF for each matrix generated in the previous step
  ##
  set.seed(42);
  
  triSig.dataM.names = data.frame(order = 1:nrow(triSig.dataM)) #just to keep order of samples for hier sig extraction
  row.names(triSig.dataM.names) = row.names(triSig.dataM)
  
  for (nFact in minK:maxK) {
    cat(sprintf("Running NMF: nFact %d (all iters)\n", nFact))
    triSig.nmfOutputByIter <- list()
    
    for (iter in 1:nBootstrapReplicates) {
      x <- triSig.dataMResamp[[iter]]
      seed = 42 + attr(x, "nIter")
      x = x[, colSums(x) > 0]
      
      #hierarchical signature extraction, remove samples that were well reconstructed, and repeat signature extraction
      for(iter.temp in 1:max.iter.hier) {
        if(nrow(x) < 10) { #there are less than 10 samples left
          break
        }
        
        iter.name = sprintf("run%03d_hier%01d", iter, iter.temp) 
        
        triSig.nmfOutputByIter[[iter.name]] <- nmf(x, r = nFact, algoritm = "mu", seed = seed + iter.temp, maxiter = 10000)
        
        triSig.dataM.estimated = triSig.nmfOutputByIter[[iter.name]]$W %*% triSig.nmfOutputByIter[[iter.name]]$H
        
        sampleSimilarities = getCosineSim(x, triSig.dataM.estimated)
        
        temp = merge(triSig.dataM.names, triSig.nmfOutputByIter[[iter.name]]$W, by = "row.names", all.x = TRUE)
        row.names(temp) = temp$Row.names
        temp = temp[order(temp$order), ]
        temp = temp[,-c(1,2)]
        temp[is.na(temp)] <- 0
        triSig.nmfOutputByIter[[iter.name]]$W = temp #to add zeros for removed samples
        
        x = x[sampleSimilarities < threshold.hier, ]
        print(paste0("iter: ", iter, ", hier iter: ", iter.temp))
      }
      rm(x)
      #triSig.nmfOutputByIter[[iter]] <- nmf(x, r = nFact, algoritm = "nsnmf", seed = seed, maxiter = 10000, parameters=list(theta=0.9))
    }
    
    idString = sprintf("nFact=%03d", nFact)
    iter.names = names(triSig.nmfOutputByIter)
    triSig.nmfOutput = triSig.nmfOutputByIter[[iter.names[1]]]
    
    if (ncol(triSig.nmfOutput$H) < ncol(triSig.dataM)) { #some cols were possibly removed because of all zeros, add them with zeros
      triSig.nmfOutput$H <- rbind.fill(as.data.frame(triSig.dataM), as.data.frame(triSig.nmfOutput$H))[c((nrow(triSig.dataM) + 1):(nrow(triSig.dataM) + nrow(triSig.nmfOutput$H))),] #very stupid, but does the trick
      triSig.nmfOutput$H[is.na(triSig.nmfOutput$H)] <- 0
    } 
    
    rownames(triSig.nmfOutput$H) = sprintf("%s_fact%02d", iter.names[1], 1:nFact)
    colnames(triSig.nmfOutput$W) = sprintf("%s_fact%02d", iter.names[1], 1:nFact)
    triSig.nmfHmatAllByFact[[idString]] <- triSig.nmfOutput$H
    triSig.nmfWmatAllByFact[[idString]] <- t(triSig.nmfOutput$W)
    triSig.nmfFrobeniusError[[idString]] <- triSig.nmfOutput$Frobenius
    triSig.nmfRMSD[[idString]] <- triSig.nmfOutput$RMSD
    
    #start from 2 here, because we already put 1 in the list above
    for (nIter in iter.names[-1]) {
      triSig.nmfOutput = triSig.nmfOutputByIter[[nIter]]
      rownames(triSig.nmfOutput$H) = sprintf("%s_fact%02d", nIter, 1:nFact)
      colnames(triSig.nmfOutput$W) = sprintf("%s_fact%02d", nIter, 1:nFact)
      if (ncol(triSig.nmfOutput$H) < ncol(triSig.nmfHmatAllByFact[[idString]])) {
        rnames <- c(rownames(triSig.nmfHmatAllByFact[[idString]]), rownames(triSig.nmfOutput$H)) #rownames are lost after rbind.fill, we need to restore them
        triSig.nmfHmatAllByFact[[idString]] <- rbind.fill(as.data.frame(triSig.nmfHmatAllByFact[[idString]]), as.data.frame(triSig.nmfOutput$H))
        triSig.nmfHmatAllByFact[[idString]][is.na(triSig.nmfHmatAllByFact[[idString]])] <- 0
        rownames(triSig.nmfHmatAllByFact[[idString]]) <- rnames
      } else {
        triSig.nmfHmatAllByFact[[idString]] = triSig.nmfHmatAllByFact[[idString]] %>% rbind(triSig.nmfOutput$H) #  same thing that you get using the coef(NMFfit) command in R   -> loadings
      }
      triSig.nmfWmatAllByFact[[idString]] = triSig.nmfWmatAllByFact[[idString]] %>% rbind(t(triSig.nmfOutput$W)) #  same thing that you get using the basis(NMFfit) command in R  -> transformed data
      triSig.nmfFrobeniusError[[idString]] <- c(triSig.nmfFrobeniusError[[idString]], triSig.nmfOutput$Frobenius)
      triSig.nmfRMSD[[idString]] <- c(triSig.nmfRMSD[[idString]], triSig.nmfOutput$RMSD)
    }
  }
  
  if(save.data) {
    saveRDS(triSig.nmfHmatAllByFact, paste0(output.name, "_Hmatrix.Rda")) #signatures
    saveRDS(triSig.nmfWmatAllByFact, paste0(output.name, "_Wmatrix.Rda")) #exposures
    saveRDS(triSig.nmfFrobeniusError, paste0(output.name, "_FrobeniusError.Rda")) #exposures
    saveRDS(triSig.nmfRMSD, paste0(output.name, "_RMSD.Rda")) #exposures
  }
  
  #triSig.nmfHmatAllByFact = readRDS(file = "test_Hmatrix.Rda")
  #triSig.nmfWmatAllByFact = readRDS(file = "test_Wmatrix.Rda")
  #triSig.nmfFrobeniusError = readRDS(file = "test_FrobeniusError.Rda")
  #triSig.nmfRMSD = readRDS(file = "test_RMSD.Rda")
  
  ##
  ## Calculate similarities between extracted signatures and a set of reference signatures (i.e., PCAWG)
  ##
  
  #order contexts alphabetically
  referenceSignatures <- referenceSignatures[, order(colnames(referenceSignatures))]
  referenceSignaturesNames <- rownames(referenceSignatures)
  
  similarities <- rep(0, nrow(referenceSignatures))
  similaritiesExposures <- rep(0, nrow(referenceSignatures))
  bestMatchesNames <- rep("", nrow(referenceSignatures))
  bestMatchesNfact <- rep("", nrow(referenceSignatures))
  bestMatches <- data.frame(matrix(nrow = nrow(referenceSignatures), ncol = ncol(referenceSignatures)))
  colnames(bestMatches) <- colnames(referenceSignatures)
  rownames(bestMatches) <- rownames(referenceSignatures)
  
  referenceSignatures.ExposureMeans.temp = referenceSignatures.ExposureMeans
  
  similarities.all = data.frame()
  
  print("Calculating similarities")
  
  for (nFact in names(triSig.nmfHmatAllByFact)) {
    sigs <- triSig.nmfHmatAllByFact[[nFact]] 
    
    if (tissueAware) {
      referenceSignatures.ExposureMeans = referenceSignatures.ExposureMeans.temp
      exps = t(triSig.nmfWmatAllByFact[[nFact]])
      exps = merge(cellLine.cancerTypes, exps, by.y = "row.names", by.x = "COSMIC_ID")
      row.names(exps) = exps$COSMIC_ID
      exps$COSMIC_ID = NULL
      exps = aggregateMatrix(exps)
      #make sure we compare the same cancer types
      exps = exps[exps$cancer_type  %in% referenceSignatures.ExposureMeans$cancer_type,]
      referenceSignatures.ExposureMeans = referenceSignatures.ExposureMeans[referenceSignatures.ExposureMeans$cancer_type %in% exps$cancer_type, ]
      exps = exps[order(exps$cancer_type), ]
      row.names(exps) = exps$cancer_type
      exps$cancer_type = NULL
      referenceSignatures.ExposureMeans = referenceSignatures.ExposureMeans[order(referenceSignatures.ExposureMeans$cancer_type), ]
      row.names(referenceSignatures.ExposureMeans) = referenceSignatures.ExposureMeans$cancer_type
      referenceSignatures.ExposureMeans$cancer_type  = NULL
    }
    
    #calculate similarity only on shared contexts between reference signatures and extracted signatures based on SNVs
    sigs = sigs[,colnames(sigs) %in% colnames(referenceSignatures)]
    #order cols alphabeticallty
    sigs = sigs[,order(colnames(sigs))]
    
    temp.similarities = matrix(ncol = 5, nrow = nrow(sigs) * nrow(referenceSignatures))
    colnames(temp.similarities) = c("nFact", "sig", "SBS", "cosine", "cosineExp")
    count = 1
    
    for(run in 1:nrow(sigs)) {
      #print(run)
      sig <- sigs[run,]
      sig.name <- row.names(sigs)[run]
      
      for (i in 1:nrow(referenceSignatures)) {
        referenceSig <- referenceSignatures[i,]
        referenceSigName <- referenceSignaturesNames[i]
        similarity <- cos.sim.vect(as.numeric(sig), as.numeric(referenceSig))
        temp.similarities[count, 1] = nFact
        temp.similarities[count, 2] = sig.name
        temp.similarities[count, 3] = referenceSigName
        temp.similarities[count, 4] = round(similarity, digits = 2)
        if (tissueAware) {
          similarityExp = cos.sim.vect(exps[, run], referenceSignatures.ExposureMeans[, i])
          temp.similarities[count, 5] = round(similarityExp, digits = 2)
        }
        count = count + 1
      }               
    }
    
    similarities.all = rbind.data.frame(similarities.all, temp.similarities, stringsAsFactors = F)
  }
  
  similarities.all = as.data.table(similarities.all)
  setkey(similarities.all, SBS)
  similarities.all.filtered = similarities.all[similarities.all$cosine >= similarityThreshold,]  
  similarities.all.filtered = similarities.all.filtered[order(-similarities.all.filtered$cosine),]
  
  Wmat.sel.all = data.frame()
  Hmat.sel.all = data.frame()
  rnames = c()
  
  if(!tissueAware) {
    tissueAWare.weight = 1
  }
  
  #find best matching signatures
  for (sbs in unique(similarities.all.filtered$SBS)) {
    temp = similarities.all.filtered[similarities.all.filtered$SBS == sbs,]
    temp = selectBestSig(temp, tissueAWare.weight, 1 - tissueAWare.weight)
    
    selectedSig = temp[temp$selected == "Selected",]
    
    Wmat.sel = triSig.nmfWmatAllByFact[[selectedSig$nFact]]
    row = Wmat.sel[row.names(Wmat.sel) == selectedSig$sig,]
    Wmat.sel.all = rbind(Wmat.sel.all, row)
    rnames = c(rnames, sbs)
    
    Hmat.sel = triSig.nmfHmatAllByFact[[selectedSig$nFact]]
    row = Hmat.sel[row.names(Hmat.sel) == selectedSig$sig,]
    Hmat.sel.all = rbind(Hmat.sel.all, row)
  }
  
  colnames(Wmat.sel.all) = colnames(Wmat.sel)
  row.names(Wmat.sel.all) = rnames
  Wmat.sel.all = t(Wmat.sel.all)
  Wmat.sel.all = round(Wmat.sel.all, digits = 5)
  
  row.names(Hmat.sel.all) = rnames
  colnames(Hmat.sel.all) = names(row)
  Hmat.sel.all = round(Hmat.sel.all, digits = 5)
  
  #normalize exposures
  Wmat.sel.all.norm = Wmat.sel.all
  for (z in 1:nrow(Wmat.sel.all.norm)) {
    Wmat.sel.all.norm[z, -1] <- round(Wmat.sel.all.norm[z, -1] / sum(Wmat.sel.all.norm[z, -1]) * 100, digits = 5)
  }
  
  #naming of the signatures
  res.cosine = data.frame()
  for (i in 1:nrow(Hmat.sel.all)) {
    for(j in 1:nrow(referenceSignatures)) {
      cosine = cos.sim.vect(as.numeric(Hmat.sel.all[i,colnames(Hmat.sel.all) %in% colnames(referenceSignatures)]), as.numeric(referenceSignatures[j,]))
      res.cosine = rbind(res.cosine, data.frame(refSignature = referenceSignaturesNames[j], Signature = row.names(Hmat.sel.all)[i], Cosine = cosine))
    }
  }
  
  #name mapping
  #res.cosine$SignaturePrimary = vapply(strsplit(as.character(res.cosine$Signature)," "), `[`, 1, FUN.VALUE=character(1))
  res.cosine$SignaturePrimary = res.cosine$Signature
  res.cosine$Cosine = round(res.cosine$Cosine, 2)
  res.cosine$refSignature = as.character(res.cosine$refSignature)
  res.cosine$Signature = as.character(res.cosine$Signature)
  res.cosine$SignaturePrimary = as.character(res.cosine$SignaturePrimary)
  
  names_new = rep("", nrow(Hmat.sel.all))
  names_old = row.names(Hmat.sel.all) 
  
  for(i in 1:nrow(Hmat.sel.all)) {
    sig = names_old[i]
    sig.cosine = res.cosine[res.cosine$Signature == sig, ]
    primarySig = sig.cosine$SignaturePrimary[1]
    sig.cosine = sig.cosine[sig.cosine$Cosine >= similarityThreshold,]
    
    if(nrow(sig.cosine) > 0) {
      sig.cosine = sig.cosine[order(-sig.cosine$Cosine),]
      sig.cosine.other = sig.cosine[sig.cosine$refSignature != primarySig,]
      refOther = sig.cosine.other$refSignature
      if(nrow(sig.cosine.other) > 0) {
        refOther = ifelse(sig.cosine.other$Cosine < 0.95, paste0(refOther, "L"), refOther)
      }
      
      if(max(sig.cosine$Cosine) < 0.95) {
        new_name = paste0("SBS", paste(str_replace(c(paste0(primarySig,"L"), refOther), "SBS", ""), collapse = "/"))    
      } else {
        new_name = paste0("SBS", paste(str_replace(c(primarySig, refOther), "SBS", ""), collapse = "/"))
      }
    } else {
      new_name = primarySig
    }
    names_new[i] = new_name
  }
  
  row.names(Hmat.sel.all) = names_new
  colnames(Wmat.sel.all.norm) = names_new
  
  write.csv(Hmat.sel.all, paste0(output.name, "_Signatures.csv")) #signatures
  write.csv(Wmat.sel.all.norm, paste0(output.name, "_Exposures.csv")) #exposures
  
}

## helper functions
cos.sim.vect = function(x, y) {
  if (!is.vector(x) || !is.vector(y)) {
    return(NA);
  } else {
    return((crossprod(x, y) / sqrt(crossprod(x) * crossprod(y)))[1, 1]);
  }
}

#gets cosine similarities between rows of matrices A and B
getCosineSim = function(matA, matB) {
  res = rep(0, times = nrow(matA))
  for(i in 1:nrow(matA)) {
    res[i] = cos.sim.vect(as.numeric(matA[i,]), as.numeric(matB[i,]))
  }
  return(res)
}

aggregateMatrix = function(mat, col = 1) {
  mat.means = data.frame(matrix(ncol = ncol(mat), nrow = length(unique(mat[, col]))))
  colnames(mat.means) = colnames(mat)
  
  for (i in 2:ncol(mat)) {
    sig = colnames(mat)[i]
    temp = mat[, c(col, i)]
    colnames(temp) = c("cancer_type", "sig")
    temp2 = aggregate(formula = sig ~ cancer_type, data = temp, FUN = mean)
    mat.means[, i] = temp2[, 2]
    mat.means[, 1] = temp2[, 1]
  }
  return(mat.means)
}

selectBestSig = function(data, cosineSig = 0.9, cosineSpectrum = 0.1) {
  data$combination = (cosineSig * as.numeric(data$cosine) + cosineSpectrum * as.numeric(data$cosineExp))
  data$selected = ""
  bestSig = which(data$combination == max(data$combination))
  if(length(bestSig) > 1) { #in case two signatures have the same score, select first one
    bestSig = bestSig[1]
  }
  data$selected[bestSig] = "Selected"
  return(data) 
}
##