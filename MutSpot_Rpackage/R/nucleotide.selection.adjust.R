#' Select nucleotide context for SNVs after adjusting threshold.
#'
#' @param nucleotide.stabs.file RDS file containing frequency of nucleotide contexts selected through lasso.
#' @param new.cutoff Updated frequency cutoff/threshold to determine nucleotide contexts used in prediction model, ranges from 0.5 to 1.
#' @param top.nucleotide Top number of nucleotide contexts to select, default = NULL.
#' @return Updated set of nucleotide contexts that passed the new threshold.
#' @export

nucleotide.selection.adjust = function(nucleotide.stabs.file, new.cutoff, top.nucleotide = NULL) {
  
  stabs = readRDS(nucleotide.stabs.file)[[1]]
  if (!is.null(top.nucleotide)) {
    
  if (sum(stabs$f >= new.cutoff) > top.nucleotide) {
    
  stabs.select = stabs[which(stabs$f >= new.cutoff), ]
  stabs.select$feature = as.character(stabs.select$feature)
  rownames(stabs.select) = NULL
  
  stabs.coef = readRDS(nucleotide.stabs.file)[[2]]
  stabs.coef$feature = as.character(stabs.coef$feature)
  rownames(stabs.coef) = NULL
  stabs.select = merge(stabs.select, stabs.coef, by = "feature")
  stabs.select$feat.coef = abs(stabs.select$feat.coef)
  stabs.select = stabs.select[order(stabs.select$f, stabs.select$feat.coef, decreasing = TRUE), ]
  stabs.select = stabs.select[1:top.nucleotide, ]
  
  sel = as.character(stabs.select$feature)
  
  } else {
    
    sel = as.character(stabs[which(stabs$f >= new.cutoff), "feature"])
    
  }
    
    } else {
    
      sel = as.character(stabs[which(stabs$f >= new.cutoff), "feature"])
    
    }
  
  return(sel)
  
}
