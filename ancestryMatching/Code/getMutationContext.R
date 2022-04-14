snvTypes <- c("ACA>A", "ACC>A", "ACG>A", "ACT>A", "CCA>A", "CCC>A", "CCG>A", "CCT>A", "GCA>A", "GCC>A", "GCG>A", "GCT>A", "TCA>A", "TCC>A", "TCG>A", "TCT>A", "ACA>T", "ACC>T", "ACG>T", "ACT>T", "CCA>T", "CCC>T", "CCG>T", "CCT>T", "GCA>T", "GCC>T", "GCG>T", "GCT>T", "TCA>T", "TCC>T", "TCG>T", "TCT>T", "ACA>G", "ACC>G", "ACG>G", "ACT>G", "CCA>G", "CCC>G", "CCG>G", "CCT>G", "GCA>G", "GCC>G", "GCG>G", "GCT>G", "TCA>G", "TCC>G", "TCG>G", "TCT>G", "AAA>C", "AAC>C", "AAG>C", "AAT>C", "CAA>C", "CAC>C", "CAG>C", "CAT>C", "GAA>C", "GAC>C", "GAG>C", "GAT>C", "TAA>C", "TAC>C", "TAG>C", "TAT>C", "AAA>T", "AAC>T", "AAG>T", "AAT>T", "CAA>T", "CAC>T", "CAG>T", "CAT>T", "GAA>T", "GAC>T", "GAG>T", "GAT>T", "TAA>T", "TAC>T", "TAG>T", "TAT>T", "AAA>G", "AAC>G", "AAG>G", "AAT>G", "CAA>G", "CAC>G", "CAG>G", "CAT>G", "GAA>G", "GAC>G", "GAG>G", "GAT>G", "TAA>G", "TAC>G", "TAG>G", "TAT>G")

get_MS <- function(chr,pos,strand,a,sep = ">", hg = "hg38"){
    require(GenomicRanges)
    require(Biostrings)
	if (hg == "hg38") {
        require(BSgenome.Hsapiens.UCSC.hg38)
    } else {
        require(BSgenome.Hsapiens.UCSC.hg19)
    }

      # build a GR object
      gr <- GRanges(
        seqnames = Rle(chr),
        ranges = IRanges(start = pos,end = pos),
        strand = Rle(strand)
      )
      seq = get_MS_GR(gr, hg = hg)
      alt_seq = DNAStringSet(a)
      # important to do reverse complementary
      reverse_mask = grepl(pattern = "[ACTG][TG][ACTG]",x = seq)
      seq[reverse_mask] <- reverseComplement(seq[reverse_mask])
      alt_seq[reverse_mask] <- reverseComplement(alt_seq[reverse_mask])
  
      seq = as.character(seq)
      alt_res = as.character(alt_seq)
  
    res = paste(seq,alt_res,sep = sep)

    return(res)
}

get_MS_GR <- function(x, extendBy = 1, hg = "hg38") {
    # x as a Granges object
    require(GenomicRanges)
    seqlevels(x) <- gsub(pattern = "chrMT",
                       replacement = "chrM",
                       x = seqlevels(x))
    tmp <- extend(x, upstream = extendBy, downstream = extendBy)

    if (hg == "hg38") {
        require(BSgenome.Hsapiens.UCSC.hg38)
        seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, tmp)
    } else {
        require(BSgenome.Hsapiens.UCSC.hg19)
        seq <- getSeq(BSgenome.Hsapiens.UCSC.hg19, tmp)
    }
    return(seq)
}

# extracted from 
# https://support.bioconductor.org/p/78652/
extend <- function(x, upstream=0, downstream=0) {
  if (any(strand(x) == "*"))
    warning("'*' ranges were treated as '+'")
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}

get_96_context = function(tmpVars) {
  mutTypesSampleCount <- table(get_MS(paste0("chr", tmpVars$Chr), tmpVars$Start, rep("+", times = nrow(tmpVars)), tmpVars$Alt, hg = "hg38"))
  contexts <- vector(length = length(snvTypes))
  names(contexts) <- snvTypes
  contexts[names(mutTypesSampleCount)] <- mutTypesSampleCount
  return(contexts)
}