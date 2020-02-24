#' Fit logistic regression prediction model for indel.
#'
#' @param mutCovariate.table.indel.file RDS file, indel covariates sparse matrix.
#' @param mutCovariate.count.indel.file RDS file, indel response matrix.
#' @param continuous.features.selected.indel.url.file Text file containing URLs of indel continuous features selected for model.
#' @param discrete.features.selected.indel.url.file Text file containing URLs of indel discrete features selected for model.
#' @param sample.specific.features.url.file Text file containing sample specific indel features, default = NULL.
#' @param fit.sparse To fit model using glmnet or glm, default = FALSE.
#' @param drop To drop insignificant features from fitted model or not, default = FALSE.
#' @param output.dir Save plot in given output directory.
#' @return Fitted indel prediction model.
#' @export

mutLRFit.indel = function(mutCovariate.table.indel.file, mutCovariate.count.indel.file, continuous.features.selected.indel.url.file, discrete.features.selected.indel.url.file, sample.specific.features.url.file = NULL, fit.sparse = FALSE, drop = FALSE, output.dir) {
  
  if (!fit.sparse) {
    
    print("Fit model using glm")
   
    # Define covariates matrix  
    mutfreq.aggregated = readRDS(file = mutCovariate.table.indel.file)
    # Define 2-column response matrix
    mutfreq.aggregated2 = readRDS(file = mutCovariate.count.indel.file)
    
    mutfreq.aggregated = as.data.frame(as.matrix(mutfreq.aggregated))
    mutfreq.aggregated = cbind(mutfreq.aggregated2, mutfreq.aggregated)
    
    if (!is.null(continuous.features.selected.indel.url.file)) {
      
      selected.continuous.urls <- read.delim(continuous.features.selected.indel.url.file, header = FALSE, stringsAsFactors = FALSE)
      selected.continuous.urls = selected.continuous.urls[ ,1]
      
    } else {
      
      selected.continuous.urls = NULL
      
    }
    
    if (!is.null(sample.specific.features.url.file)) {
      
      sample.specific.urls <- read.delim(sample.specific.features.url.file, stringsAsFactors = FALSE)
      continuous.sample.specific = discrete.sample.specific = NULL
      for (j in 1:ncol(sample.specific.urls)) {
  
        if (class(sample.specific.urls[ ,j]) != "character"){
          
          continuous.sample.specific = c(continuous.sample.specific, colnames(sample.specific.urls)[j])
          
        } else {
          
          discrete.sample.specific = c(discrete.sample.specific, colnames(sample.specific.urls)[j])
        }
        
      }

    } else {
      
      continuous.sample.specific = NULL
      discrete.sample.specific = NULL
      
    }
    
    for(i in colnames(mutfreq.aggregated)[which(!colnames(mutfreq.aggregated) %in% c("mut.count", "nonmut.count", "ind.mut.count", selected.continuous.urls, continuous.sample.specific))]) {
      
      mutfreq.aggregated[ ,i] = as.character(mutfreq.aggregated[ ,i])
      
    }
    
    # Fit model using glm with non-sparse matrix
    LRmodel <- stats::glm(cbind(mut.count, nonmut.count) ~ ., family = binomial(logit), data = mutfreq.aggregated, x = F, y = F)
    
    if (drop) {
      
      ignore = c("ind.mut.count")
      
      # Remove features that are not significant
      pval = summary(LRmodel)$coef[ ,4]
      pval = pval[-1]
      
      if (!is.null(sample.specific.features.url.file)) {
        
      special.feat = read.delim(sample.specific.features.url.file, stringsAsFactors = FALSE)[ ,-1]
      
      for (k in colnames(special.feat)) {
        
        if (class(special.feat[ ,k]) == "character") {
          
          if (sum(pval > 0.05 & grepl(k, names(pval))) != sum(grepl(k, names(pval)))) {
            
            ignore = c(ignore, names(pval)[(grepl(k, names(pval)))])
            
          }
          
        }
        
      }
      
      }
      
      if (sum(pval > 0.05 & !names(pval) %in% ignore) >= 1) {
        
        refit = TRUE
        
      } else {
        
        refit = FALSE
        
      }
      
      while(refit) {
        
        rm(LRmodel)
        gc(reset = TRUE)
        
        print(paste("Remove ", names(pval)[which(pval > 0.05 & !names(pval) %in% ignore)], sep = ""))
        pval = which(pval > 0.05 & !names(pval) %in% ignore)
        
        mutfreq.aggregated = mutfreq.aggregated[ ,-(pval + 2)]
        print("Refitting model...")
        LRmodel <- stats::glm(cbind(mut.count, nonmut.count) ~ ., family = binomial(logit), data = mutfreq.aggregated, x = F, y = F)
        
        pval = summary(LRmodel)$coef[ ,4]
        pval = pval[-1]
        
        if (!is.null(sample.specific.features.url.file)) {
          
          special.feat = read.delim(sample.specific.features.url.file, stringsAsFactors = FALSE)[ ,-1]
          
          for (k in colnames(special.feat)) {
            
            if (class(special.feat[ ,k]) == "character") {
              
              if (sum(pval > 0.05 & grepl(k, names(pval))) != sum(grepl(k, names(pval)))) {
                
                ignore = c(ignore, names(pval)[(grepl(k, names(pval)))])
                
              }
              
            }
            
          }
          
        }
        
        if (sum(pval > 0.05 & !names(pval) %in% ignore) >= 1) {
          
          refit = TRUE
          
        } else {
          
          refit = FALSE
          
        }
        
      }
      
      # Redefine selected features for prediction
      # features = rownames(summary(LRmodel)$coef)[-1]
      features = colnames(mutfreq.aggregated)[!colnames(mutfreq.aggregated) %in% c("mut.count", "nonmut.count")]
      
      # Continuous features
      if (!is.null(continuous.features.selected.indel.url.file)) {
        
        continuous.features = read.delim(continuous.features.selected.indel.url.file, stringsAsFactors = FALSE, header = FALSE)
        if (sum(continuous.features[ ,1] %in% features) != nrow(continuous.features)) {
          
          rem = which(!continuous.features[ ,1] %in% features)
          continuous.features = continuous.features[-rem, ]
          if (nrow(continuous.features) == 0) {
            
            continuous.features = NULL
            
          }
          
        } else {
          
          continuous.features = "unchanged"
          
        } 
        
        } else {
          
          continuous.features = "unchanged"
          
        }
      
      # Discrete features
      if (!is.null(discrete.features.selected.indel.url.file)) {
        
        discrete.features = read.delim(discrete.features.selected.indel.url.file, stringsAsFactors = FALSE, header = FALSE)
        if (sum(discrete.features[ ,1] %in% features) != nrow(discrete.features)) {
          
          rem = which(!discrete.features[ ,1] %in% features)
          discrete.features = discrete.features[-rem, ]
          if (nrow(discrete.features) == 0) {
            
            discrete.features = NULL
            
          }
          
        } else {
          
          discrete.features = "unchanged"
          
        } 
        
      } else {
        
        discrete.features = "unchanged"
        
      }
      
      # Sample specific features
      if (!is.null(sample.specific.features.url.file)) {
        
        sample.specific = read.delim(sample.specific.features.url.file, stringsAsFactors = FALSE)
        rem = NULL
        sample.feat = colnames(sample.specific)
        sample.feat = sample.feat[which(sample.feat != "SampleID")]
        for (p in sample.feat) {
          
          if (sum(grepl(p, features)) == 0) {
            rem = c(rem, p)
            
          }
          
        }
        if (length(rem) > 0) {
          
          sample.specific = sample.specific[ ,-which(colnames(sample.specific) %in% rem)]
          sample.specific = as.data.frame(sample.specific)
          
          if (ncol(sample.specific) == 1) {
            
            sample.specific = NULL
            
          }
          
        } else {
          
          sample.specific = "unchanged"
          
          } 
        
        } else {
          
          sample.specific = "unchanged"
          
        }
      
      # Plot feature importance barplot
      print("Plot feature importance")
      plot_feature_importance(LRmodel = LRmodel, mutCovariate.table.file = mutCovariate.table.indel.file, mutation.type = "indel", output.dir = output.dir, feature.names = c(discrete.sample.specific, colnames(mutfreq.aggregated)))
      
      return(list(stripGlmLR(LRmodel), continuous.features, discrete.features, sample.specific))
      
    } else {
      
      # Plot feature importance barplot
      print("Plot feature importance")
      plot_feature_importance(LRmodel = LRmodel, mutCovariate.table.file = mutCovariate.table.indel.file, mutation.type = "indel", output.dir = output.dir, feature.names = c(discrete.sample.specific, colnames(mutfreq.aggregated)))
      
      return(stripGlmLR(LRmodel))
    }
    
  } else {
    
    print("Fit model using glmnet")
    
  # Define covariates matrix  
  mutfreq.aggregated = readRDS(file = mutCovariate.table.indel.file)
  # Define 2-column response matrix
  mutfreq.aggregated2 = readRDS(file = mutCovariate.count.indel.file)
  
  myTryCatch <- function(expr) {
    
    warn <- err <- NULL
    value <- withCallingHandlers(
      
      tryCatch(expr, error = function(e) {
        
        err <<- e
        NULL
        
      } ), warning = function(w) {
        
        warn <<- w
        invokeRestart("muffleWarning")
        
      } )
    
    list(value = value, warning = warn, error = err)
    
  }
  
  print("Fitting model...")
  # Fit model using glmnet with sparse matrix and check for convergence issue
  LRmodel <- myTryCatch(glmnet::glmnet(x = mutfreq.aggregated, y = as.matrix(mutfreq.aggregated2[ ,c("nonmut.count", "mut.count")]), alpha = 1, lambda = 0, family = "binomial"))
  
  if(!is.null(LRmodel$warning)) {
    
    print("Convergence issue hence fitting model with glm")
    
    mutfreq.aggregated = as.data.frame(as.matrix(mutfreq.aggregated))
    mutfreq.aggregated = cbind(mutfreq.aggregated2, mutfreq.aggregated)
    
    if (!is.null(continuous.features.selected.indel.url.file)) {
      
      selected.continuous.urls <- read.delim(continuous.features.selected.indel.url.file, header = FALSE, stringsAsFactors = FALSE)
      selected.continuous.urls = selected.continuous.urls[ ,1]
      
    } else {
      
      selected.continuous.urls = NULL
      
    }
    
    if (!is.null(sample.specific.features.url.file)) {
      
      sample.specific.urls <- read.delim(sample.specific.features.url.file, stringsAsFactors = FALSE)
      continuous.sample.specific = discrete.sample.specific = NULL
      for (j in 1:ncol(sample.specific.urls)) {
        
        if (class(sample.specific.urls[ ,j]) != "character"){
          
          continuous.sample.specific = c(continuous.sample.specific, colnames(sample.specific.urls)[j])
          
        } else {
          
          discrete.sample.specific = c(discrete.sample.specific, paste(colnames(sample.specific.urls)[j], unique(sample.specific.urls[ ,j]), sep = ""))
          
        }
        
      }
      
    } else {
      
      continuous.sample.specific = NULL
      discrete.sample.specific = NULL
      
    }
    
    for(i in colnames(mutfreq.aggregated)[which(!colnames(mutfreq.aggregated) %in% c("mut.count", "nonmut.count", "ind.mut.count", selected.continuous.urls, continuous.sample.specific))]) {
      
      mutfreq.aggregated[ ,i] = as.character(mutfreq.aggregated[ ,i])
      
    }
    
    rm(LRmodel)
    gc(reset = TRUE)
    
    # Fit model using glm with non-sparse matrix
    LRmodel <- stats::glm(cbind(mut.count, nonmut.count) ~ ., family = binomial(logit), data = mutfreq.aggregated, x = F, y = F)
    
    if (drop) {
      
      ignore = c("ind.mut.count")
      
      # Remove features that are not significant
      pval = summary(LRmodel)$coef[ ,4]
      pval = pval[-1]
      
      if (!is.null(sample.specific.features.url.file)) {
        
      special.feat = read.delim(sample.specific.features.url.file, stringsAsFactors = FALSE)[ ,-1]
      
      for (k in colnames(special.feat)) {
        
        if (class(special.feat[ ,k]) == "character") {
          
          if (sum(pval > 0.05 & grepl(k, names(pval))) != sum(grepl(k, names(pval)))) {
            
            ignore = c(ignore, names(pval)[(grepl(k, names(pval)))])
            
          }
          
        }
        
      }
      
      }
      
      if (sum(pval > 0.05 & !names(pval) %in% ignore) >= 1) {
        
        refit = TRUE
        
      } else {
        
        refit = FALSE
        
      }
      
      while(refit) {
        
        rm(LRmodel)
        gc(reset = TRUE)
        
        print(paste("Remove ", names(pval)[which(pval > 0.05 & !names(pval) %in% ignore)], sep = ""))
        pval = which(pval > 0.05 & !names(pval) %in% ignore)
        
        mutfreq.aggregated = mutfreq.aggregated[ ,-(pval + 2)]
        print("Refitting model...")
        LRmodel <- stats::glm(cbind(mut.count, nonmut.count) ~ ., family = binomial(logit), data = mutfreq.aggregated, x = F, y = F)
        
        pval = summary(LRmodel)$coef[ ,4]
        pval = pval[-1]
  
        if (!is.null(sample.specific.features.url.file)) {
          
          special.feat = read.delim(sample.specific.features.url.file, stringsAsFactors = FALSE)[ ,-1]
          
          for (k in colnames(special.feat)) {
            
            if (class(special.feat[ ,k]) == "character") {
              
              if (sum(pval > 0.05 & grepl(k, names(pval))) != sum(grepl(k, names(pval)))) {
                
                ignore = c(ignore, names(pval)[(grepl(k, names(pval)))])
                
              }
              
            }
            
          }
          
        }
        
        if (sum(pval > 0.05 & !names(pval) %in% ignore) >= 1) {
          
          refit = TRUE
          
        } else {
          
          refit = FALSE
          
        }
        
      }
      
      # Redefine selected features for prediction
      # features = rownames(summary(LRmodel)$coef)[-1]
      features = colnames(mutfreq.aggregated)[!colnames(mutfreq.aggregated) %in% c("mut.count", "nonmut.count")]
      
      # Continuous features
      if (!is.null(continuous.features.selected.indel.url.file)) {
        
        continuous.features = read.delim(continuous.features.selected.indel.url.file, stringsAsFactors = FALSE, header = FALSE)
        if (sum(continuous.features[ ,1] %in% features) != nrow(continuous.features)) {
          
          rem = which(!continuous.features[ ,1] %in% features)
          continuous.features = continuous.features[-rem, ]
          if (nrow(continuous.features) == 0) {
            
            continuous.features = NULL
            
          }
          
        } else {
          
          continuous.features = "unchanged"
          
        }
        
        } else {
          
          continuous.features = "unchanged"
          
        }
      
      # Discrete features
      if (!is.null(discrete.features.selected.indel.url.file)) {
        
        discrete.features = read.delim(discrete.features.selected.indel.url.file, stringsAsFactors = FALSE, header = FALSE)
        if (sum(discrete.features[ ,1] %in% features) != nrow(discrete.features)) {
          
          rem = which(!discrete.features[ ,1] %in% features)
          discrete.features = discrete.features[-rem, ]
          if (nrow(discrete.features) == 0) {
            
            discrete.features = NULL
            
          }
          
        } else {
          
          discrete.features = "unchanged"
          
        }
        
      } else {
        
        discrete.features = "unchanged"
        
      }
      
      # Sample specific features
      if (!is.null(sample.specific.features.url.file)) {
        
        sample.specific = read.delim(sample.specific.features.url.file, stringsAsFactors = FALSE)
        rem = NULL
        sample.feat = colnames(sample.specific)
        sample.feat = sample.feat[which(sample.feat != "SampleID")]
        for (p in sample.feat) {
          
          if (sum(grepl(p, features)) == 0) {
            rem = c(rem, p)
            
          }
          
        }
        if (length(rem) > 0) {
          
          sample.specific = sample.specific[ ,-which(colnames(sample.specific) %in% rem)]
          sample.specific = as.data.frame(sample.specific)
          
          if (ncol(sample.specific) == 1) {
            
            sample.specific = NULL
            
          }
          
        } else {
          
          sample.specific = "unchanged"
          
          } 
        
        } else {
          
          sample.specific = "unchanged"
          
        }
      
      # Plot feature importance barplot
      print("Plot feature importance")
      plot_feature_importance(LRmodel = LRmodel, mutCovariate.table.file = mutCovariate.table.indel.file, mutation.type = "indel", output.dir = output.dir, feature.names = c(discrete.sample.specific, colnames(mutfreq.aggregated)))
      
      return(list(stripGlmLR(LRmodel), continuous.features, discrete.features, sample.specific))
      
    } else {
      
      # Plot feature importance barplot
      print("Plot feature importance")
      plot_feature_importance(LRmodel = LRmodel, mutCovariate.table.file = mutCovariate.table.indel.file, mutation.type = "indel", output.dir = output.dir, feature.names = c(discrete.sample.specific, colnames(mutfreq.aggregated)))
      
      return(stripGlmLR(LRmodel))
      
    }
    
  } else {
    
    if (!is.null(sample.specific.features.url.file)) {
      
      sample.specific.urls <- read.delim(sample.specific.features.url.file, stringsAsFactors = FALSE)
      continuous.sample.specific = discrete.sample.specific = NULL
      for (j in 1:ncol(sample.specific.urls)) {
        
        if (class(sample.specific.urls[ ,j]) != "character"){
          
          continuous.sample.specific = c(continuous.sample.specific, colnames(sample.specific.urls)[j])
          
        } else {
          
          discrete.sample.specific = c(discrete.sample.specific, colnames(sample.specific.urls)[j])
        }
        
      }
      
    } else {
      
      continuous.sample.specific = NULL
      discrete.sample.specific = NULL
      
    }
    
    # Plot feature importance barplot
    print("Plot feature importance")
    plot_feature_importance(LRmodel = LRmodel[[1]], mutCovariate.table.file = mutCovariate.table.indel.file, mutation.type = "indel", output.dir = output.dir, feature.names = c(discrete.sample.specific, colnames(mutfreq.aggregated)))
    
    return(LRmodel[[1]])
    
  }
  
}

}
