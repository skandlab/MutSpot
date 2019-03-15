#' Calculate maximum score in each bin.
#'
#' @param bins Genomic bins GRanges.
#' @param numvar Feature scores across genome in RleList.
#' @param mcolname Colname of feature score, default = "score".
#' @return Maximum score for each bin.
#' @export

binnedMax = function(bins, numvar, mcolname) {
  
  stopifnot(is(bins, "GRanges"))
  stopifnot(is(numvar, "RleList"))
  stopifnot(identical(sort(GenomeInfoDb::seqlevels(bins)), sort(names(numvar))))
  bins.per.chrom <- GenomicRanges::split(IRanges::ranges(bins), as.character(GenomicRanges::seqnames(bins)))
  means.list <- lapply(GenomeInfoDb::seqlevels(bins),
                       function(seqname) {
                         
                         views <- IRanges::Views(numvar[[seqname]], bins.per.chrom[[seqname]])
                         IRanges::viewMaxs(views)
                         
                       })
  new.mcol <- suppressWarnings(BiocGenerics::unsplit(means.list, as.factor(as.character(GenomicRanges::seqnames(bins)))))
  S4Vectors::mcols(bins)[[mcolname]] <- new.mcol
  bins
  
}
