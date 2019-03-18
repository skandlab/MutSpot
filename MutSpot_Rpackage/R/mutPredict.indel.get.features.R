#' Mutation hotspot recurrence prediction for indel.
#' 
#' @param roi Each indel hotspot in GRanges.
#' @param continuous.features.indel Full scores list of indel continuous features selected for model.
#' @param discrete.features.indel Full scores list of indel discrete features selected for model.
#' @return Feature matrix.
#' @export

mutPredict.indel.get.features <- function(roi, continuous.features.indel, discrete.features.indel) {
  
  # Get a list of sites for each region
  sites = IRanges::tile(roi, width = 1)
  sites = BiocGenerics::unlist(sites)  
  
  # Extract DNA sequence for +/-10bp around each site in roi
  seq = GenomicRanges::GRanges(GenomeInfoDb::seqnames(sites), IRanges::IRanges(IRanges::ranges(sites)@start - 5, IRanges::ranges(sites)@start + 5))
  seq = IRanges::Views(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, seq)
  seq = as.data.frame(seq)
  
  # Assign continuous epigenetic scores to each site in roi
  if (!is.null(continuous.features.indel)) {
    
    for (i in names(continuous.features.indel)) {
      
      df = data.frame(numeric(length(sites)))
      colnames(df) = i
      GenomicRanges::values(sites) = data.frame(GenomicRanges::values(sites), df, check.names = FALSE)
      names(GenomicRanges::mcols(sites))[length(names(GenomicRanges::mcols(sites)))] = i
      ovl = IRanges::findOverlaps(sites, continuous.features.indel[[i]])
      if (length(ovl) != 0) {
        
        GenomicRanges::values(sites)[S4Vectors::queryHits(ovl), i] = continuous.features.indel[[i]][S4Vectors::subjectHits(ovl)]$score
        
      }
      rm(ovl)
      rm(df)
      
    }
    
  }
  
  # Assign discrete epigenetic scores to each site
  if (!is.null(discrete.features.indel)) {
    
    for (i in names(discrete.features.indel)) {
      
      df = data.frame(numeric(length(sites)))
      colnames(df) = i
      GenomicRanges::values(sites) = data.frame(GenomicRanges::values(sites), df, check.names = FALSE)
      names(GenomicRanges::mcols(sites))[length(names(GenomicRanges::mcols(sites)))] = i
      ovl = IRanges::findOverlaps(sites, discrete.features.indel[[i]])
      if (length(ovl) != 0) {
        
        GenomicRanges::values(sites)[S4Vectors::queryHits(ovl), i] = 1
        
      }
      rm(ovl)
      
    }
    
  }
  
  roi.features = as.data.frame(GenomicRanges::mcols(sites))
  colnames(roi.features) = c(names(continuous.features.indel), names(discrete.features.indel))
  out = roi.features

  # Assign polyA/T/C/G scores to each site
  feat.polymer = mutPredict.indel.find.polymer(seq)
  rownames(feat.polymer) = NULL
  out = cbind(roi.features, feat.polymer)
  
  return(out)
  
}
