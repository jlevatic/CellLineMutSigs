##
## Example of running ancestry matching
##

#Note, set the working directory to the ancestryMatching folder: setwd("yourPath/ancestryMatching)
source("Code/ancestryMatching.R")

PCA.components = read.csv("Data/PCAcomponents.csv", row.names = 1, check.names = F)

#small sample of the data to run the example
PCA.components.subset = rbind(PCA.components[1:100,], PCA.components[1000:1500,])

contexts.96.CCLs <- read.csv("Data/contexts.96.CCLs.csv", row.names = 1, check.names = F)
contexts.96.TCGA.Germline <- read.csv("Data/contexts.96.TCGA.Germline.csv", row.names = 1, check.names = F)

outFile = "contexts.96.CCLs.ancestryMatching.csv"

ancestryMatching(output.file = output.file, PCA.components = PCA.components.subset, contexts.96.CCLs = contexts.96.CCLs, contexts.96.TCGA.Germline = contexts.96.TCGA.Germline, nClusters = 4)

##
## The following code would reproduce the results from the paper 
##
# PCA.components.subset = rbind(PCA.components[1:100,], PCA.components[1000:1500,])
# contexts.96.CCLs <- read.csv("Data/contexts.96.CCLs.csv", row.names = 1, check.names = F)
# contexts.96.TCGA.Germline <- read.csv("Data/contexts.96.TCGA.Germline.csv", row.names = 1, check.names = F)
# variants.CCLs = readRDS("Data/variants.CCLs.RDS")
# variants.TCGA.Germline = ## Due to privacy we are not allowed to share the germline variants of TCGA patients ##
# outFile = "contexts.96.CCLs.ancestryMatching.csv"

# ancestryMatching(output.file = output.file, nClusters = 13, PCA.components = PCA.components, contexts.96.CCLs = contexts.96.CCLs, contexts.96.TCGA.Germline = contexts.96.TCGA.Germline, remove.CNVkit.deletions = T, CNVkit.folder = "Data/CNVkit_cns/", CNVkit.log2.threshold = -0.3, variants.CCLs = variants.CCLs, variants.TCGA.Germline = variants.TCGA.Germline)