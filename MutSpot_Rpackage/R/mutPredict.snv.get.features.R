#' Mutation hotspot recurrence prediction for SNV.
#' 
#' @param roi Each SNV hotspot in GRanges.
#' @param continuous.features.snv Full scores list of SNV continuous features selected for model.
#' @param discrete.features.snv Full scores list of SNV discrete features selected for model.
#' @param motifs Selected nucleotide contexts to be extracted.
#' @return Feature matrix.
#' @export

mutPredict.snv.get.features <- function(roi, continuous.features.snv, discrete.features.snv, motifs) {
  
  # Get a list of sites for each roi
  sites = as.list(GenomicRanges::tile(roi, width = 1))
  sites = do.call(c, sites)
  
  # Extract DNA sequence for +/-10bp around each site in roi
  seq = GenomicRanges::GRanges(GenomicRanges::seqnames(sites), IRanges::IRanges(IRanges::ranges(sites)@start - 5, IRanges::ranges(sites)@start + 5))
  seq = IRanges::Views(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, seq)
  seq = as.data.frame(seq)   
  
  # Assign continuous epigenetic scores to each site in roi
  if (!is.null(continuous.features.snv)) {
    
    for (i in names(continuous.features.snv)) {
      
      df = data.frame(numeric(length(sites)))
      colnames(df) = i
      GenomicRanges::values(sites) = data.frame(GenomicRanges::values(sites), df, check.names = FALSE)
      names(GenomicRanges::mcols(sites))[length(names(GenomicRanges::mcols(sites)))] = i
      ovl = IRanges::findOverlaps(sites, continuous.features.snv[[i]])
      if (length(ovl) != 0) {
        
        GenomicRanges::values(sites)[S4Vectors::queryHits(ovl), i] = continuous.features.snv[[i]][S4Vectors::subjectHits(ovl)]$score
        
        }
      rm(ovl)
      rm(df)
      
    }
    
    }
  
  # Assign discrete epigenetic scores to each site
  if (!is.null(discrete.features.snv)) {
    
    for (i in names(discrete.features.snv)) {
      
      df = data.frame(numeric(length(sites)))
      colnames(df) = i
      GenomicRanges::values(sites) = data.frame(GenomicRanges::values(sites), df, check.names = FALSE)
      names(GenomicRanges::mcols(sites))[length(names(GenomicRanges::mcols(sites)))] = i
      ovl = IRanges::findOverlaps(sites, discrete.features.snv[[i]])
      if (length(ovl) != 0) {
        
        GenomicRanges::values(sites)[S4Vectors::queryHits(ovl), i] = 1
        
        }
      rm(ovl)
      
    }
    
    }
  
  roi.features = as.data.frame(GenomicRanges::mcols(sites))
  colnames(roi.features) = c(names(continuous.features.snv), names(discrete.features.snv))
  out = roi.features
  
  # Assign nucleotide context scores to each site
  if (!is.null(motifs)) {
    
    roi.motif <- mutPredict.snv.find.motif(seq, motifs)
    out = cbind(roi.features, roi.motif)
    
    }
  
  return(out)
  
}
