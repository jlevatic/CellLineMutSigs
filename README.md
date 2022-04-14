# CellLineMutSigs
Code associated with the paper "Mutational signatures are markers of drug sensitivity of cancer cells" by Levatic et al. (2022)

## Overview

The folders (ancestryMatching, signatureExtraction, and twoWayAssociationTesting) in contain code, data, examples, and expected output of examples of custom scripts used for the purposes of the manuscript "Mutational signatures are markers of drug sensitivity of cancer cells". The scripts are based on R programming language. The main methods are implemented as R functions. Input parameters of functions are documented within the scripts.

## ancestryMatching

Contains data, script, and example for running the code for filtering of germline variants based on the "ancestry matching" method. 

Note that, due to privacy we are not allowed to share the germline variants of TCGA patients which are required to use some functionalities of the script (namely, removal of variants in deleted regions as determined by the output of the CNVkit software). The script can be used without this functionality.

### System requirements:

R (tested on version 3.6.2)

The script was tested on the following OS: Windows 10 Pro x64

### Hardware requirements:

A standard computer with enough RAM to support the in-memory operations.

### Installation guide:

There are no specific installation requirements other than installing R and the R packages the script depends on.

The script is loaded into R environment as follows: ```source("ancestryMatching.R")```

### R dependencies:

tclust (tested on version 1.4-2)
pcaPP (tested on version 1.9-74)
foreach (tested on version 1.5.0)
doParallel (tested on version 1.0.16)
Additional optional R dependencies if remove.CNVkit.deletions option is used:
GenomicRanges (tested on version 1.38.0)
Biostrings (tested on version 2.54.0)
BSgenome.Hsapiens.UCSC.hg38

#### Expected running time of example: <1min
