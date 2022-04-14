##
## Example of running the two way association testing procedure
##

#Note, set the working directory to the twoWayAssociationTesting.R folder: setwd("yourPath/twoWayAssociationTesting.R)
source("Code/twoWayAssociationTesting.R")

## Example of running GDSC/PSCORE test on "Dabrafenib" (x-dimension) and "BRAF" (y-dimension) for CFEs features (i.e., oncogenic mutations, copy number alterations, and DNA hypermethylation)
twoWayAssociations(folder.data = "Data/", folder.output = "CFEs_PSCORE", path.to.features = "Data/CFEs_features.csv", analysis = "pscore", xaxis = "Dabrafenib", yaxis = "BRAF", nRandomizations = 100)
## Example of running GDSC/PRISM test on "5-fluorouracil" (x- and y-dimension) for CFEs
twoWayAssociations(folder.data = "Data/", folder.output = "CFEs_PRISM", path.to.features = "Data/CFEs_features.csv", analysis = "prism", xaxis = "5-fluorouracil", yaxis = "5-fluorouracil", nRandomizations = 100)
## Example of running GDSC/GDSC (same target) test on "Cisplatin" (x-dimension) and "Oxaliplatin" (y-dimension) for CFEs
twoWayAssociations(folder.data = "Data/", folder.output = "CFEss_GDSC", path.to.features = "Data/CFEs_features.csv", analysis = "same_target", xaxis = "Cisplatin", yaxis = "Oxaliplatin", nRandomizations = 100)

## Example of running GDSC/PSCORE test on "Dabrafenib" (x-dimension) and "BRAF" (y-dimension) for mutational signatures
twoWayAssociations(folder.data = "Data/", folder.output = "MutSigs_PSCORE", path.to.features = "Data/signaturesExposures_Binarized.csv", analysis = "pscore", xaxis = "Dabrafenib", yaxis = "BRAF", nRandomizations = 100)
## Example of running GDSC/PRISM test on "5-fluorouracil" (x- and y-dimension) for mutational signatures
twoWayAssociations(folder.data = "Data/", folder.output = "MutSigs_PRISM", path.to.features = "Data/signaturesExposures_Binarized.csv", analysis = "prism", xaxis = "5-fluorouracil", yaxis = "5-fluorouracil", nRandomizations = 100)
## Example of running GDSC/GDSC (same target) test on "Cisplatin" (x-dimension) and "Oxaliplatin" (y-dimension) for mutational signatures
twoWayAssociations(folder.data = "Data/", folder.output = "MutSigs_GDSC", path.to.features = "Data/signaturesExposures_Binarized.csv", analysis = "same_target", xaxis = "Cisplatin", yaxis = "Oxaliplatin", nRandomizations = 100)


##
## The following code would reproduce the results from the paper 
## Note, this is to exemplify how to run the code, this is time consuming and it is better to run it in parallel rather than in a loop as shown here.
##

## Running all GDSC/PRISM pairs for mutational signatures and CFEs
# gdsc_prism_pairs = read.csv("Data/pairs_gdsc_prism.txt", header = F)
# for(i in 1:nrow(gdsc_prism_pairs)) {
#   twoWayAssociations(folder.data = "Data/", folder.output = "MutSigs_PRISM", path.to.features = "Data/signaturesExposures_Binarized.csv", analysis = "prism", xaxis = gdsc_prism_pairs$V1[i], yaxis = gdsc_prism_pairs$V2[i], nRandomizations = 100000, eff.size.th = 0.2)
#   twoWayAssociations(folder.data = "Data/", folder.output = "CFEs_PRISM", path.to.features = "Data/CFEs_features.csv", analysis = "prism", xaxis = gdsc_prism_pairs$V1[i], yaxis = gdsc_prism_pairs$V2[i], nRandomizations = 100000, eff.size.th = 0.2)
# }

## Running all GDSC/PSCORE pairs for mutational signatures and CFEs.
# gdsc_pscore_pairs = read.csv("Data/pairs_gdsc_pscore.txt", header = F)
# for(i in 1:nrow(gdsc_pscore_pairs)) {
#   twoWayAssociations(folder.data = "Data/", folder.output = "MutSigs_PSCORE", path.to.features = "Data/signaturesExposures_Binarized.csv", analysis = "pscore", xaxis = gdsc_pscore_pairs$V1[i], yaxis = gdsc_pscore_pairs$V2[i], nRandomizations = 100000, eff.size.th = 0.2)
#   twoWayAssociations(folder.data = "Data/", folder.output = "CFEs_PSCORE", path.to.features = "Data/CFEs_features.csv", analysis = "pscore", xaxis = gdsc_pscore_pairs$V1[i], yaxis = gdsc_pscore_pairs$V2[i], nRandomizations = 100000, eff.size.th = 0.2)
# }

## Running all GDSC/GDSC (same target) pairs for mutational signatures and CFEs.
# gdsc_gdsc_pairs = read.csv("Data/pairs_gdsc_gdsc_same_target.txt", header = F)
# for(i in 1:nrow(gdsc_gdsc_pairs)) {
#   twoWayAssociations(folder.data = "Data/", folder.output = "MutSigs_SameTarget", path.to.features = "Data/signaturesExposures_Binarized.csv", analysis = "same_target", xaxis = gdsc_gdsc_pairs$V1[i], yaxis = gdsc_gdsc_pairs$V2[i], nRandomizations = 100000, eff.size.th = 0.2)
#   twoWayAssociations(folder.data = "Data/", folder.output = "CFEs_SameTarget", path.to.features = "Data/CFEs_features.csv", analysis = "same_target", xaxis = gdsc_gdsc_pairs$V1[i], yaxis = gdsc_gdsc_pairs$V2[i], nRandomizations = 100000, eff.size.th = 0.2)
# }
