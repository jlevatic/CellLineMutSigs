##
## Example of running signature extraction
##
#Note, set the working directory to the signatureExtraction folder: setwd("yourPath/signatureExtraction)
source("Code/signatureExtraction.R")

contexts = read.csv("C:/Users/jlevatic/Desktop/cellLines_contexts.csv", check.names = F, row.names = 1)

referenceSignatures = read.csv("Data/referenceSignatures_PCAWG.csv", check.names = F, row.names = 1)
referenceSignatures.ExposureMeans = read.csv("Data/referenceSignaturesExposureMeans.csv")

#Note that the code relies on GPU NMF (nmfgpu4R R package), therefore, to use it Cuda have to be installed (with appropriate GPU that supports it)
extractSignatures(contexts = contexts, output.name = "test", nBootstrapReplicates = 10, maxK = 5, threshold.hier = 0.9, max.iter.hier = 2, referenceSignatures = referenceSignatures, referenceSignatures.exposureMeans = referenceSignatures.exposureMeans, tissueAware = T, tissueAWare.weight = 0.9, similarityThreshold = 0.85)

##
## The following code would reproduce the results from the paper 
##
# contexts = read.csv("C:/Users/jlevatic/Desktop/cellLines_contexts.csv", check.names = F, row.names = 1)
# referenceSignatures = read.csv("Data/referenceSignatures_PCAWG.csv", check.names = F, row.names = 1)
# referenceSignatures.ExposureMeans = read.csv("Data/referenceSignaturesExposureMeans.csv")
#
# extractSignatures(contexts = contexts, output.name = "signatures", nBootstrapReplicates = 300, maxK = 40, minK = 2, threshold.hier = 0.97, max.iter.hier = 3, referenceSignatures = referenceSignatures, referenceSignatures.exposureMeans = referenceSignatures.exposureMeans, tissueAware = T, tissueAWare.weight = 0.9, similarityThreshold = 0.85)

