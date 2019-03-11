#' Remove Region A from Region B for 1 chromosome.
#' 
#' @param rois Region B, GRanges object.
#' @param sub Region A, GRangesList/GRanges object.
#' @param chr Chromosome name.
#' @return Filtered Region B for 1 chromosome.
#' @export

subtract.regions.from.roi.chr <- function(rois, sub,chr = 'chr10') {
  
  if (class(rois) == "GRangesList") {
    
    roi.chr = BiocGenerics::unlist(rois[as.character(GenomicRanges::seqnames(rois)) == chr]); sub.chr = BiocGenerics::unlist(sub[as.character(GenomicRanges::seqnames(sub)) == chr]);
    GenomicRanges::elementMetadata(roi.chr) <- NULL; GenomicRanges::elementMetadata(sub.chr) <- NULL
    
    # Return if one list is empty
    if (length(roi.chr) == 0 || length(sub.chr) == 0) {
      
      nms = names(roi.chr)
      names(roi.chr) <- NULL
      rl = split(roi.chr, nms)
      return(rl)
      
    }
    
    r1 = c(roi.chr, sub.chr)
    r2 = IRanges::disjoin(r1)
    
    ## Problem here is that setdiff merges adjacent ranges ...
    #r3 = setdiff(r2,sub.chr)
    
    r3 = r2[-S4Vectors::queryHits(IRanges::findOverlaps(r2, sub.chr))]
    rovl = IRanges::findOverlaps(rois, r3)
    r4 = tapply(r3[S4Vectors::subjectHits(rovl)], S4Vectors::queryHits(rovl), function(xs) {IRanges::reduce(BiocGenerics::unlist(xs))})
    r4.names = tapply(names(rois[S4Vectors::queryHits(rovl)]), S4Vectors::queryHits(rovl), function(x) x[1])
    names(r4) = r4.names
    rl = GenomicRanges::GRangesList(r4);
    rl
    
  } else if (class(rois)=="GRanges") {
    
    roi.chr = rois[GenomicRanges::seqnames(rois) == chr]; sub.chr = sub[GenomicRanges::seqnames(sub) == chr];
    GenomicRanges::elementMetadata(roi.chr) <- NULL; GenomicRanges::elementMetadata(sub.chr) <- NULL
    # return if one list is empty
    if (length(roi.chr) == 0 || length(sub.chr) == 0) {
      
      return(roi.chr)
      
    }
    
    r1 = c(roi.chr,sub.chr)
    r2 = GenomicRanges::disjoin(r1)
    
    ## Problem here is that setdiff merges adjacent ranges ...
    #r3 = setdiff(r2,sub.chr)
    
    r3 = r2[-S4Vectors::queryHits(IRanges::findOverlaps(r2, sub.chr))]
    r3
    
  }
}
