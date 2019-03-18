#' Remove Region A from Region B for 1 chromosome.
#' 
#' @param rois Region B, GRanges object.
#' @param sub Region A, GRangesList/GRanges object.
#' @param chr Chromosome name.
#' @return Filtered Region B for 1 chromosome.
#' @export

subtract.regions.from.roi.chr <- function(rois, sub,chr = 'chr10') {
  
  if (class(rois) == "GRangesList") {
    
    roi.chr = BiocGenerics::unlist(rois[as.character(seqnames(rois)) == chr]); sub.chr = BiocGenerics::unlist(sub[as.character(seqnames(sub)) == chr]);
    elementMetadata(roi.chr) <- NULL; elementMetadata(sub.chr) <- NULL
    
    # Return if one list is empty
    if (length(roi.chr) == 0 || length(sub.chr) == 0) {
      
      nms = names(roi.chr)
      names(roi.chr) <- NULL
      rl = split(roi.chr, nms)
      return(rl)
      
    }
    
    r1 = c(roi.chr, sub.chr)
    r2 = disjoin(r1)
    
    ## Problem here is that setdiff merges adjacent ranges ...
    #r3 = setdiff(r2,sub.chr)
    
    r3 = r2[queryHits(findOverlaps(r2, sub.chr))]
    rovl = findOverlaps(rois, r3)
    r4 = tapply(r3[subjectHits(rovl)], queryHits(rovl), function(xs) { reduce(BiocGenerics::unlist(xs)) })
    r4.names = tapply(names(rois[queryHits(rovl)]), queryHits(rovl), function(x) x[1])
    names(r4) = r4.names
    rl = GRangesList(r4);
    rl
    
  } else if (class(rois)=="GRanges") {
    
    roi.chr = rois[seqnames(rois) == chr]; sub.chr = sub[seqnames(sub) == chr];
    elementMetadata(roi.chr) <- NULL; elementMetadata(sub.chr) <- NULL
    
    # Return if one list is empty
    if (length(roi.chr) == 0 || length(sub.chr) == 0) {
      
      return(roi.chr)
      
    }
    
    r1 = c(roi.chr,sub.chr)
    r2 = disjoin(r1)
    
    ## Problem here is that setdiff merges adjacent ranges ...
    #r3 = setdiff(r2,sub.chr)
    
    r3 = r2[-queryHits(findOverlaps(r2, sub.chr))]
    r3
    
  }
}
