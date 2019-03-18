#' Remove Region A from Region B.
#'
#' @param rois Region B, GRanges object.
#' @param sub Region A, GRangesList/GRanges object.
#' @param cores Number of cores, default = 1.
#' @return Filtered Region B.
#' @export

subtract.regions.from.roi <- function(rois, sub, cores = 1) {
  
  # In parallel for all chromosomes
  chromosomes <- unique(as.character(GenomeInfoDb::seqnames(BiocGenerics::unlist(rois))))
  rois.subtracted <- parallel::mclapply(chromosomes, function(chr) {subtract.regions.from.roi.chr(rois, sub, chr)},mc.cores = cores, mc.preschedule = FALSE)
  rois.subtracted <- do.call(c, rois.subtracted)
  GenomeInfoDb::seqlevels(rois.subtracted) = GenomeInfoDb::seqlevels(rois)
  rois.subtracted
  
}
