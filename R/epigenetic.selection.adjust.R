#' Select epigentic features through LASSO.
#'
#' @param feature.stabs.snv.file RDS file containing frequency of SNV epigenetic features selected through lasso.
#' @param feature.stabs.indel.file RDS file containing frequency of indel epigenetic features selected through lasso.
#' @param continuous.features.selected.snv.url.file Text file containing URLs of selected continuous SNV epigentic features.
#' @param discrete.features.selected.snv.url.file Text file containing URLs of selected discrete SNV epigenetic features.
#' @param continuous.features.selected.indel.url.file Text file containing URLs of selected continuous indel epigenetic features.
#' @param discrete.features.selected.indel.url.file Text file containing URLs of selected discrete indel epigenetic features.
#' @param new.cutoff.snv Updated frequency cutoff/threshold to determine SNV epigenetic features used in prediction model, ranges from 0.5 to 1.
#' @param new.cutoff.indel Updated frequency cutoff/threshold to determine indel epigenetic features used in prediction model, ranges from 0.5 to 1.
#' @return Updated set of SNV/indel continuous and discrete features that passed the new threshold.
#' @export

epigenetic.selection.adjust = function(feature.stabs.snv.file, feature.stabs.indel.file, continuous.features.selected.snv.url.file, discrete.features.selected.snv.url.file, continuous.features.selected.indel.url.file, discrete.features.selected.indel.url.file, new.cutoff.snv, new.cutoff.indel) {
  
  # If SNV mutations available, else skip this
  if (!is.null(feature.stabs.snv.file) & !is.null(new.cutoff.snv)) {
    
    stabs = readRDS(feature.stabs.snv.file)
    sel = as.character(stabs[which(stabs$f >= new.cutoff.snv), "feature"])
    
    # If continuous SNV features selected before, else return NULL
    if (!is.null(continuous.features.selected.snv.url.file)) {
      
      continuous.snv.features = read.delim(continuous.features.selected.snv.url.file, header = FALSE, stringsAsFactors = FALSE)
      continuous.snv.features = continuous.snv.features[which(continuous.snv.features[ ,1] %in% sel), ] 
      
    } else {
        
      continuous.snv.features = NULL
      
    }
    
    # If discrete SNV features selected before, else return NULL
    if (!is.null(discrete.features.selected.snv.url.file)) {
      
      discrete.snv.features = read.delim(discrete.features.selected.snv.url.file, header = FALSE, stringsAsFactors = FALSE)
      discrete.snv.features = discrete.snv.features[which(discrete.snv.features[ ,1] %in% sel), ]
      
    } else {
        
      discrete.snv.features = NULL
      
    }
    
  } else {
    
    continuous.snv.features = NULL
    discrete.snv.features = NULL
    
  }
  
  # If indel mutations available, else skip this
  if (!is.null(feature.stabs.indel.file) & !is.null(new.cutoff.indel)) {
    
    stabs = readRDS(feature.stabs.indel.file)
    sel = as.character(stabs[which(stabs$f >= new.cutoff.indel), "feature"])
    
    # If continuous indel features selected before, else return NULL
    if (!is.null(continuous.features.selected.indel.url.file)) {
      
      continuous.indel.features = read.delim(continuous.features.selected.indel.url.file, header = FALSE, stringsAsFactors = FALSE)
      continuous.indel.features = continuous.indel.features[which(continuous.indel.features[ ,1] %in% sel), ]
      
    } else {
        
      continuous.indel.features = NULL
      
    }
    
    # If discrete features selected before, else return NULL
    if (!is.null(discrete.features.selected.indel.url.file)) {
      
      discrete.indel.features = read.delim(discrete.features.selected.indel.url.file, header = FALSE, stringsAsFactors = FALSE)
      discrete.indel.features = discrete.indel.features[which(discrete.indel.features[ ,1] %in% sel), ]
      
    } else {
        
      discrete.indel.features = NULL
      
    }
    
  } else {
    
    continuous.indel.features = NULL
    discrete.indel.features = NULL
    
  }
  
  return(list(continuous.snv.features, discrete.snv.features, continuous.indel.features, discrete.indel.features))
  
}
