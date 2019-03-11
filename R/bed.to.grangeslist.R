#' Converts bed file to grangeslist.
#'
#' @param fname Name of bed file.
#' @return GRangeslist of the input file.
#' @export

bed.to.grangeslist = function(fname) {
  roit <- read.table(fname,header=F)[,1:4]
  colnames(roit) <- c('chr','start','end','id')
  # create GRanges object 
  roi.gr <- with(roit, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end), id=as.character(id)))
  # split by ID into GRangesList
  roi.grl <- GenomicRanges::split(roi.gr,roi.gr$id)
}
