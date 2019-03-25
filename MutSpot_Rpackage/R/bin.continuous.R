#' Divide continuous feature scores into specified number of bins.
#'
#' @param feature.name Continuous feature name.
#' @param feature.url URL of continuous feature.
#' @param nbins Number of bins to divide continuous feature score into.
#' @return Binned continuous feature scores.
#' @export

bin.continuous = function(feature.name, feature.url, nbins) {
  
  
  chrOrder <- c(paste("chr", 1:22, sep = ""), "chrX")
  seqi = GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[GenomeInfoDb::seqnames(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[1:23]]
  seqnames = GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens))[1:23]
  
  z = rtracklayer::import(feature.url)
  
  mutrateBins = quantile(z$score, seq(0, 1, 1 / nbins))
  mutrateBins[1] = mutrateBins[1]-1e-5
  mutrate.bin = .bincode(z$score, mutrateBins)
  bin.mean = lapply(1:nbins, function(x) mean(z$score[mutrate.bin == x]))

  feat.bins = .bincode(z$score, mutrateBins)
  out = sapply(feat.bins, function(x) bin.mean[x])
  
  z$score = unlist(out)
  z = GenomicRanges::as.data.frame(z)
  
  z = z[ ,c(1:3, 6)]
  z = split(z, f = z$seqnames)
  
  for (i in 1:length(z)) {
    
    chr = as.character(unique(z[[i]]$seqnames))
    z[[i]][which(z[[i]]$start == min(z[[i]]$start)), "start"] <- 1
    z[[i]][which(z[[i]]$end == max(z[[i]]$end)), "end"] = GenomeInfoDb::seqlengths(seqi)[chr]
    
  }
  z = do.call(rbind, z)
  rownames(z) = NULL
  
  # Adjust start coordinates to 0-based as in bed files
  z$start = z$start - 1
  
return(z)

}
