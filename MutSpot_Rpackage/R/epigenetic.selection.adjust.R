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
#' @param top.features Number of top genomic features to select, default = NULL.
#' @param features.sds RDS list containing standard deviations of each feature.
#' @return Updated set of SNV/indel continuous and discrete features that passed the new threshold.
#' @export

epigenetic.selection.adjust = function(feature.stabs.snv.file, feature.stabs.indel.file, continuous.features.selected.snv.url.file, discrete.features.selected.snv.url.file, continuous.features.selected.indel.url.file, discrete.features.selected.indel.url.file, new.cutoff.snv, new.cutoff.indel, top.features = NULL, features.sds) {
  
  # If SNV mutations available, else skip this
  if (!is.null(feature.stabs.snv.file) & !is.null(new.cutoff.snv)) {
    
    stabs = readRDS(feature.stabs.snv.file)[[1]]
    if (!is.null(top.features)) {
      
      if (sum(stabs$f >= new.cutoff.snv) > top.features) {
      
    # stabs.select = stabs[which(stabs$f >= new.cutoff.snv), ]
    stabs.select = stabs
    stabs.select$feature = as.character(stabs.select$feature)
    rownames(stabs.select) = NULL
    
    stabs.coef = readRDS(feature.stabs.snv.file)[[2]]
    stabs.coef$feature = as.character(stabs.coef$feature)
    rownames(stabs.coef) = NULL
    stabs.select = merge(stabs.select, stabs.coef, by = "feature")
    
    stabs.sds = readRDS(features.sds)[[1]]
    stabs.sds = as.data.frame(stabs.sds)
    stabs.sds$feature = as.character(rownames(stabs.sds))
    rownames(stabs.sds) = NULL
    stabs.select = merge(stabs.select, stabs.sds, by = "feature")
    
    stabs.select$feat.coef2 = stabs.select$feat.coef * stabs.select$stabs.sds
    stabs.select$feat.coef2 = abs(stabs.select$feat.coef2)
    stabs.select = stabs.select[order(stabs.select$f, stabs.select$feat.coef2, decreasing = TRUE), ]
    stabs.select = stabs.select[1:top.features, ]
    
    sel = as.character(stabs.select$feature)
    
      } else {
        
        sel = as.character(stabs[which(stabs$f >= new.cutoff.snv), "feature"])
        
      }
      
    } else {
      
      sel = as.character(stabs[which(stabs$f >= new.cutoff.snv), "feature"])
        
    }
    
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
    
    stabs = readRDS(feature.stabs.indel.file)[[1]]
    if (!is.null(top.features)) {
      
      if (sum(stabs$f >= new.cutoff.indel) > top.features) {
        
        stabs.select = stabs[which(stabs$f >= new.cutoff.indel), ]
        stabs.select$feature = as.character(stabs.select$feature)
        rownames(stabs.select) = NULL
        
        stabs.coef = readRDS(feature.stabs.indel.file)[[2]]
        stabs.coef$feature = as.character(stabs.coef$feature)
        rownames(stabs.coef) = NULL
        stabs.select = merge(stabs.select, stabs.coef, by = "feature")
        
        stabs.sds = readRDS(features.sds)[[2]]
        stabs.sds = as.data.frame(stabs.sds)
        stabs.sds$feature = as.character(rownames(stabs.sds))
        rownames(stabs.sds) = NULL
        stabs.select = merge(stabs.select, stabs.sds, by = "feature")
        
        stabs.select$feat.coef2 = stabs.select$feat.coef * stabs.select$stabs.sds
        stabs.select$feat.coef2 = abs(stabs.select$feat.coef2)
        stabs.select = stabs.select[order(stabs.select$f, stabs.select$feat.coef2, decreasing = TRUE), ]
        stabs.select = stabs.select[1:top.features, ]
        
        sel = as.character(stabs.select$feature)
        
      } else {
        
        sel = as.character(stabs[which(stabs$f >= new.cutoff.indel), "feature"])
        
      }
      
    } else {
      
      sel = as.character(stabs[which(stabs$f >= new.cutoff.indel), "feature"])
      
    }
    
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
