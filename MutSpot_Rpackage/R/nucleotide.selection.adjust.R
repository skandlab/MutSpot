#' Select nucleotide context for SNVs after adjusting threshold.
#'
#' @param nucleotide.stabs.file RDS file containing frequency of nucleotide contexts selected through lasso.
#' @param new.cutoff Updated frequency cutoff/threshold to determine nucleotide contexts used in prediction model, ranges from 0.5 to 1.
#' @return Updated set of nucleotide contexts that passed the new threshold.
#' @export

nucleotide.selection.adjust = function(nucleotide.stabs.file, new.cutoff) {
  
  stabs = readRDS(nucleotide.stabs.file)
  sel = as.character(stabs[which(stabs$f >= new.cutoff), "feature"])
  
  return(sel)
  
}
