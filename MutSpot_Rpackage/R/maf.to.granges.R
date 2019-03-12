#' Convert MAF file to GRanges object.
#' 
#' @param fname.maf MAF file.
#' @return GRanges object containing information from MAF file.
#' @examples
#' maf.to.granges("SNV_mutations.MAF")
#' maf.to.granges("Indel_mutations.MAF")
#' @export

maf.to.granges <- function(fname.maf) {
  
  print(">> Reading compact MAF ...")
  maf.t = read.table(fname.maf, stringsAsFactors = FALSE, header = FALSE)
  colnames(maf.t)[1:6] = c('chr', 'start', 'end', 'reference.allele', 'tumor.allele', 'sampleID')
  # Create MAF GRanges object
  maf <- with(maf.t, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end), ral = reference.allele, tal = tumor.allele, sid = sampleID))
  
}
