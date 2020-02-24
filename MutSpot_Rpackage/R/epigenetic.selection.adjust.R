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
#' @param genomic.features.snv Text file containing URLs of potential continuous and discrete SNV epigenetic features to select from, default = NULL.
#' @param genomic.features.indel Text file containing URLs of potential continuous and discrete indel epigenetic features to select from, default = NULL.
#' @param genomic.features Text file containing URLs of potential continuous and discrete SNV and indel epigenetic features to select from, default = NULL.
#' @param feature.dir Directory containing binned feature bed files.
#' @return Updated set of SNV/indel continuous and discrete features that passed the new threshold.
#' @export

epigenetic.selection.adjust = function(feature.stabs.snv.file, feature.stabs.indel.file, continuous.features.selected.snv.url.file, discrete.features.selected.snv.url.file, continuous.features.selected.indel.url.file, discrete.features.selected.indel.url.file, new.cutoff.snv, new.cutoff.indel, top.features = NULL, features.sds, genomic.features.snv = genomic.features.snv, genomic.features.indel = genomic.features.indel, genomic.features = genomic.features,feature.dir=feature.dir) {
  
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
    
    if (!is.null(genomic.features)) {
      
      genomic.features.snv = genomic.features
      
    }
    
    features = read.delim(genomic.features.snv, stringsAsFactors = FALSE)
    if (sum(features$feature_type == 1) != 0) {
      
    full.continuous.snv = features[which(features$feature_type == 1), ]
    full.continuous.snv$file_path = paste(feature.dir, full.continuous.snv$feature_name, ".bed", sep="")
    full.continuous.snv = rbind(full.continuous.snv, c("local_mutrate", paste(feature.dir, "localmutrate_snv.bed", sep = ""), 1, 10))
    
    } else {
      
      full.continuous.snv = NULL
      
    }
    
    if (sum(features$feature_type != 1) != 0) {
      
    full.discrete.snv = features[which(features$feature_type != 1), ]
    
    } else {
      
      full.discrete.snv = NULL
      
    }
    
    # If continuous SNV features selected before, else return NULL
    if (sum(sel %in% full.continuous.snv$feature_name) > 0) {
      
      continuous.snv.features = full.continuous.snv[which(full.continuous.snv$feature_name %in% sel), 1:2]

    } else {
        
      continuous.snv.features = NULL
      
    }
    
    # If discrete SNV features selected before, else return NULL
    if (sum(sel %in% full.discrete.snv$feature_name) > 0) {
      
      discrete.snv.features = full.discrete.snv[which(full.discrete.snv$feature_name %in% sel), 1:2]

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
    
    if (!is.null(genomic.features)) {
      
      genomic.features.indel = genomic.features
      
    }
    
    features = read.delim(genomic.features.indel, stringsAsFactors = FALSE)
    if (sum(features$feature_type == 1) != 0) {
      
      full.continuous.indel = features[which(features$feature_type == 1), ]
      full.continuous.indel$file_path = paste(feature.dir, full.continuous.indel$feature_name, ".bed", sep="")
      full.continuous.indel = rbind(full.continuous.indel, c("local_mutrate", paste(feature.dir, "localmutrate_indel.bed", sep = ""), 1, 10))
      
    } else {
      
      full.continuous.indel = NULL
      
    }
    
    if (sum(features$feature_type != 1) != 0) {
      
    full.discrete.indel = features[which(features$feature_type != 1), ]
    
    } else {
      
      full.discrete.indel = NULL
      
    }
    
    # If continuous indel features selected before, else return NULL
    if (sum(sel %in% full.continuous.indel$feature_name) > 0) {
      
      continuous.indel.features = full.continuous.indel[which(full.continuous.indel$feature_name %in% sel), 1:2]

    } else {
        
      continuous.indel.features = NULL
      
    }
    
    # If discrete features selected before, else return NULL
    if (sum(sel %in% full.discrete.indel$feature_name) > 0) {
      
      discrete.indel.features = full.discrete.indel[which(full.discrete.indel$feature_name %in% sel), 1:2]

    } else {
        
      discrete.indel.features = NULL
      
    }
    
  } else {
    
    continuous.indel.features = NULL
    discrete.indel.features = NULL
    
  }
  
  return(list(continuous.snv.features, discrete.snv.features, continuous.indel.features, discrete.indel.features))
  
}
