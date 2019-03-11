#' Summarize bigWig into bins.
#' 
#' @param biwWigUrl URL of bigWig file.
#' @param bins Genomic bins GRanges.
#' @param type Calculate mean or maximum score for each bin, default = "mean".
#' @return Mean/Maximum score for each bin.
#' @export

bigwig.summarize.bins = function(bigWigUrl,bins, type = 'mean') {
  
  print(paste0('>>', bigWigUrl))
  ss.cov = rtracklayer::import(as.character(bigWigUrl), as = 'RleList')
  ss.cov = ss.cov[BiocGenerics::intersect(names(ss.cov), GenomeInfoDb::seqlevels(bins))]
  # Add empty missing chromosomes
  for (c in BiocGenerics::setdiff(GenomeInfoDb::seqlevels(bins), names(ss.cov))) { ss.cov[[c]] = Rle() }
  if (type == 'mean') {
    
    # Compute binned Mean
    ss.bins = binnedMean(bins, ss.cov, 'score')
    
  } else if (type == 'max') {
    
    # Compute binned Max
    ss.bins = binnedMax(bins, ss.cov, 'score')     
    
  }
  # Force zeros if -Inf (happens if no data for genome.bin)
  S4Vectors::values(ss.bins)$score[is.infinite(ss.bins$score)] <- 0
  S4Vectors::values(ss.bins)$score[is.nan(ss.bins$score)] <- 0
  ss.bins$score
  
}
