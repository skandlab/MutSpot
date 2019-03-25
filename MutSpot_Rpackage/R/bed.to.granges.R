#' Convert BED file to GRanges object.
#' 
#' @param bed.file BED file to be converted.
#' @return GRanges object.
#' @export

bed.to.granges <- function(bed.file) {
  
  z = read.delim(bed.file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  z$V2 = z$V2 + 1
  if (ncol(z) == 3) {
    
  z = with(z, GenomicRanges::GRanges(V1, IRanges::IRanges(V2, V3)))
  
  } else {
    
    z = with (z,GenomicRanges::GRanges(V1, IRanges::IRanges(V2, V3), score = V4))
    
  }

  return(z)
  
} 
