#' Select epigentic features through LASSO.
#'
#' @param sampled.sites.snv.file Sampled SNV mutated and non-mutated sites RDS file.
#' @param sampled.sites.indel.file Sampled indel mutated and non-mutated sites RDS file.
#' @param genomic.features.snv Text file containing URLs of potential continuous and discrete SNV epigenetic features to select from, default = NULL.
#' @param genomic.features.indel Text file containing URLs of potential continuous and discrete indel epigenetic features to select from, default = NULL.
#' @param genomic.features Text file containing URLs of potential continuous and discrete SNV and indel epigenetic features to select from, default = NULL.
#' @param genomic.features.fixed.snv Text file containing URLs of fixed continuous and discrete SNV epigenetic features, default = NULL.
#' @param genomic.features.fixed.indel Text file containing URLs of fixed continuous and discrete indel epigenetic features, default = NULL.
#' @param genomic.features.fixed Text file containing URLs of fixed continuous and discrete SNV and indel epigenetic features, default = NULL.
#' @param cores Number of cores, default = 1.
#' @param cutoff Frequency cutoff/threshold to determine epigenetic features used in prediction model, ranges from 0.5 to 1, default = 0.75.
#' @param feature.dir Directory containing binned feature bed files.
#' @param genome.build Reference genome build, default = Ch37.
#' @return A list containing frequency of SNV epigenetic features selected through lasso, SNV continuous and discrete features that passed the threshold, frequency of indel epigenetic features selected through lasso, indel continuous and discrete features that passed the threshold.
#' @export

epigenetic.selection = function(sampled.sites.snv.file, sampled.sites.indel.file, genomic.features.snv = NULL, genomic.features.indel = NULL, genomic.features = NULL, genomic.features.fixed.snv = NULL, genomic.features.fixed.indel = NULL, genomic.features.fixed = NULL, cores = 1, cutoff = 0.75, feature.dir, genome.build = "Ch37") {
  
  # Bin bigwig files
  if (!is.null(genomic.features)) {
    
    # Define all potential epigenetic features
    all.features = read.delim(genomic.features, stringsAsFactors = FALSE, header = TRUE)
    
    # Define continuous epigenetic features, if exists
    if (sum(all.features$feature_type == 1) > 0) {
      
      continuous.features = all.features[which(all.features$feature_type == 1), ]
      continuous.features.urls = continuous.features  
      
      # Bin features 
      if (sum(!is.na(continuous.features$nbins)) > 0) {
        
        to.bin = which(!is.na(continuous.features$nbins))
        
        # Bin continuous features based on number of bins provided
        for (x in to.bin) {
          
          print(paste("Binning ", continuous.features[x, "feature_name"], " into ", continuous.features[x, "nbins"], " bins", sep = ""))
          
          continuous.features.binned = bin.continuous(feature.name = continuous.features[x, "feature_name"], feature.url = continuous.features[x, "file_path"], nbins = continuous.features[x, "nbins"], genome.build = genome.build)
          
          write.table(continuous.features.binned, file = paste(feature.dir, continuous.features[x, "feature_name"], ".bed", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
          
          continuous.features.urls[x, "file_path"] <- paste(feature.dir, continuous.features[x, "feature_name"], ".bed", sep = "")
          all.features[which(all.features$feature_name == continuous.features.urls[x, "feature_name"]), "file_path"] <- paste(feature.dir,continuous.features[x, "feature_name"], ".bed", sep = "")
          
        }
        
      }
      
    }
      
    genomic.features.snv = genomic.features.indel = all.features
    
    } else if (!is.null(genomic.features.snv)) {
    
      # Define all potential epigenetic features
      all.features = read.delim(genomic.features.snv, stringsAsFactors = FALSE, header = TRUE)
      
      # Define continuous epigenetic features, if exists
      if (sum(all.features$feature_type == 1) > 0) {
        
        continuous.features = all.features[which(all.features$feature_type == 1), ]
        continuous.features.urls = continuous.features  
        
        # Bin features 
        if (sum(!is.na(continuous.features$nbins)) > 0) {
          
          to.bin = which(!is.na(continuous.features$nbins))
          
          # Bin continuous features based on number of bins provided
          for (x in to.bin) {
            
            print(paste("Binning ", continuous.features[x, "feature_name"], " into ", continuous.features[x, "nbins"], " bins", sep = ""))
            
            continuous.features.binned = bin.continuous(feature.name = continuous.features[x, "feature_name"], feature.url = continuous.features[x, "file_path"], nbins = continuous.features[x, "nbins"], genome.build = genome.build)
            
             write.table(continuous.features.binned, file = paste(feature.dir, continuous.features[x, "feature_name"], ".bed", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
            
            continuous.features.urls[x, "file_path"] <- paste(feature.dir, continuous.features[x, "feature_name"], ".bed", sep = "")
            all.features[which(all.features$feature_name == continuous.features.urls[x, "feature_name"]), "file_path"] <- paste(feature.dir, continuous.features[x, "feature_name"], ".bed", sep = "")
            
          }
          
        }
        
      }
      
        genomic.features.snv = all.features
      
    } else if (!is.null(genomic.features.indel)) {
    
      # Define all potential epigenetic features
      all.features = read.delim(genomic.features.indel, stringsAsFactors = FALSE, header = TRUE)
      
      # Define continuous epigenetic features, if exists
      if (sum(all.features$feature_type == 1) > 0) {
        
        continuous.features = all.features[which(all.features$feature_type == 1), ]
        continuous.features.urls = continuous.features  
        
        # Bin features 
        if (sum(!is.na(continuous.features$nbins)) > 0) {
          
          to.bin = which(!is.na(continuous.features$nbins))
          
          # Bin continuous features based on number of bins provided
          for (x in to.bin) {
            
            print(paste("Binning ", continuous.features[x, "feature_name"], " into ", continuous.features[x, "nbins"], " bins", sep = ""))
            
            continuous.features.binned = bin.continuous(feature.name = continuous.features[x, "feature_name"], feature.url = continuous.features[x, "file_path"], nbins = continuous.features[x, "nbins"], genome.build = genome.build)
            
            write.table(continuous.features.binned, file = paste(feature.dir, continuous.features[x, "feature_name"], ".bed", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
            
            continuous.features.urls[x, "file_path"] <- paste(feature.dir, continuous.features[x, "feature_name"], ".bed", sep = "")
            all.features[which(all.features$feature_name == continuous.features.urls[x, "feature_name"]), "file_path"] <- paste(feature.dir,continuous.features[x, "feature_name"],".bed",sep="")
            
          }
          
        }
        
      }
        
        genomic.features.indel = all.features
      
  }
  
  
  # If SNV mutations available, else skip this
  if (!is.null(sampled.sites.snv.file)) {
    
    # If other epigenetic files provided, run lasso selection, else select only local mutation rate
    if (!is.null(genomic.features.snv)) {
    
    create_lasso_input = create.lasso.input(sampled.sites.file = sampled.sites.snv.file, genomic.features.file = genomic.features.snv, mutation.type = "snv", cores = cores)
    
    # sd of each column for standardization later
    sds.snv = apply(as.matrix(create_lasso_input[[1]]), 2, sd)
    
    stabs.site <- stability.sel.epigenetic(x_data = create_lasso_input[[1]], y_data = create_lasso_input[[2]], threshold = cutoff, cores = cores) 
    
    feature = names(stabs.site$freq.1se)
    f = stabs.site$freq.1se
    freq = data.frame(feature, f)
    sel = freq
    sel = sel[order(sel$f, decreasing = TRUE), ]
    sel = as.character(sel[which(sel$f >= cutoff), "feature"])
    
    # Since there will always be continuous feature (local mutation rate)
    features.continuous.url.snv = create_lasso_input[[3]][which(create_lasso_input[[3]][ ,1] %in% sel), ]
    
    # Selected discrete feature urls
    if (!is.null(create_lasso_input[[4]])) {
      
    features.discrete.url.snv = create_lasso_input[[4]][which(create_lasso_input[[4]][ ,1] %in% sel), ]
    
    } else {
      
      features.discrete.url.snv <- data.frame(matrix(ncol = 2, nrow = 0))
      x <- c("V1", "V2")
      colnames(features.discrete.url.snv) <- x
      
    }
    freq.snv = freq
    rm(freq)
    
    feat.coef = colMeans(stabs.site$coef.1se)
    feat.coef = data.frame(feature, feat.coef)
    feat.coef.snv = feat.coef
    rm(feat.coef)

    } else {
      
      # If no epigenetic features provided, will automatically choose to keep local mutation rate
      freq.snv = NULL
      features.continuous.url.snv = data.frame(V1 = "local_mutrate", V2 = "./features/localmutrate_snv.bed")
      features.discrete.url.snv = NULL
      feat.coef.snv = NULL
      sds.snv = NULL
      
    }
    
  } else {
    
    freq.snv = NULL
    features.continuous.url.snv = NULL
    features.discrete.url.snv = NULL
    feat.coef.snv = NULL
    sds.snv = NULL
    
  }
  
  # If indel mutations available, else skip this
  if (!is.null(sampled.sites.indel.file)) {
    
    # If other epigenetic files provided, run lasso selection, else select only local mutation rate 
    if (!is.null(genomic.features.indel)) {
      
    create_lasso_input = create.lasso.input(sampled.sites.file = sampled.sites.indel.file, genomic.features.file = genomic.features.indel, mutation.type = "indel", cores = cores)
    
    # sd of each column for standardization later
    sds.indel = apply(as.matrix(create_lasso_input[[1]]), 2, sd)
    
    stabs.site <- stability.sel.epigenetic(x_data = create_lasso_input[[1]], y_data = create_lasso_input[[2]], threshold = 0.75, cores = cores) 
    
  feature = names(stabs.site$freq.1se)
  f = stabs.site$freq.1se
  freq = data.frame(feature, f)
  sel = freq
  sel = sel[order(sel$f, decreasing = TRUE), ]
  sel = as.character(sel[which(sel$f >= cutoff), "feature"])
  
  # Since there will always be continuous feature (local mutation rate)
  features.continuous.url.indel = create_lasso_input[[3]][which(create_lasso_input[[3]][ ,1] %in% sel), ]
  
  if (!is.null(create_lasso_input[[4]])) {
    
    features.discrete.url.indel = create_lasso_input[[4]][which(create_lasso_input[[4]][ ,1] %in% sel), ]
  
  } else {
    
    features.discrete.url.indel <- data.frame(matrix(ncol = 2, nrow = 0))
    x <- c("V1", "V2")
    colnames(features.discrete.url.indel) <- x
    
  }
  freq.indel = freq
  rm(freq)
  
  feat.coef = colMeans(stabs.site$coef.1se)
  feat.coef = data.frame(feature, feat.coef)
  feat.coef.indel = feat.coef
  rm(feat.coef)
  
    } else {
      
      # If no epigenetic features provided, will automatically choose to keep local mutation rate
      freq.indel = NULL
      features.continuous.url.indel = data.frame(V1 = "local_mutrate", V2 = "./features/localmutrate_indel.bed")
      features.discrete.url.indel = NULL
      feat.coef.indel = NULL
      sds.indel = NULL
      
    }
  
  } else {
    
    freq.indel = NULL
    features.continuous.url.indel = NULL
    features.discrete.url.indel = NULL
    feat.coef.indel = NULL
    sds.indel = NULL
    
  }
  
  # Consider fixed features that do not undergo selection
  if (!is.null(genomic.features.fixed)) {
    
    all.fixed = read.delim(genomic.features.fixed, stringsAsFactors = FALSE, header = TRUE)    
    if (sum(all.fixed$feature_type == 1) > 0) {
      
      fixed.continuous = all.fixed[which(all.fixed$feature_type == 1), ]
      fixed.continuous.urls = fixed.continuous  
      
      # Bin features 
      if (sum(!is.na(fixed.continuous$nbins)) > 0) {
        
        to.bin = which(!is.na(fixed.continuous$nbins))
        
        # Bin continuous features based on number of bins provided
        for (x in to.bin) {
          
          print(paste("Binning ", fixed.continuous[x, "feature_name"], " into ", fixed.continuous[x, "nbins"], " bins", sep = ""))
          
          fixed.continuous.binned = bin.continuous(feature.name = fixed.continuous[x, "feature_name"], feature.url = fixed.continuous[x, "file_path"], nbins = fixed.continuous[x, "nbins"], genome.build = genome.build)
          
          write.table(fixed.continuous.binned, file = paste(feature.dir, fixed.continuous[x, "feature_name"], ".bed", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
          
          fixed.continuous.urls[x, "file_path"] <- paste(feature.dir, fixed.continuous[x, "feature_name"], ".bed", sep = "")
          
        }
        
      }
      
      fixed.continuous.urls = fixed.continuous.urls[ ,c("feature_name", "file_path")]
      rownames(fixed.continuous.urls) = NULL
      if (!is.null(features.continuous.url.snv)) {
        
        rownames(features.continuous.url.snv) = NULL
        features.continuous.url.snv = rbind(features.continuous.url.snv, fixed.continuous.urls)
        
      } else {
        
        features.continuous.url.snv = fixed.continuous.urls
        
      }
      if (!is.null(features.continuous.url.indel)) {
        
        rownames(features.continuous.url.indel) = NULL
        features.continuous.url.indel = rbind(features.continuous.url.indel, fixed.continuous.urls)
        
      } else {
        
        features.continuous.url.indel = fixed.continuous.urls
        
      }
      
    } 
    
    if (sum(all.fixed$feature_type == 0) > 0) {
      
      fixed.discrete = all.fixed[which(all.fixed$feature_type == 0), 1:2]
      rownames(fixed.discrete) = NULL
      if (!is.null(features.discrete.url.snv)) {
        
        rownames(features.discrete.url.snv) = NULL
        features.discrete.url.snv = rbind(features.discrete.url.snv, fixed.discrete)
        
      } else {
        
        features.discrete.url.snv = fixed.discrete
        
      }
      if (!is.null(features.discrete.url.indel)) {
        
        rownames(features.discrete.url.indel) = NULL
        features.discrete.url.indel = rbind(features.discrete.url.indel, fixed.discrete)
        
      } else {
        
        features.discrete.url.indel = fixed.discrete
        
      }
      
    }
    
  } else if (!is.null(genomic.features.fixed.snv)) {
    
    all.fixed = read.delim(genomic.features.fixed.snv, stringsAsFactors = FALSE, header = TRUE)    
    if (sum(all.fixed$feature_type == 1) > 0) {
      
      fixed.continuous = all.fixed[which(all.fixed$feature_type == 1), ]
      fixed.continuous.urls = fixed.continuous  
      
      # Bin features 
      if (sum(!is.na(fixed.continuous$nbins)) > 0) {
        
        to.bin = which(!is.na(fixed.continuous$nbins))
        
        # Bin continuous features based on number of bins provided
        for (x in to.bin) {
          
          print(paste("Binning ", fixed.continuous[x, "feature_name"], " into ", fixed.continuous[x, "nbins"], " bins", sep = ""))
          
          fixed.continuous.binned = bin.continuous(feature.name = fixed.continuous[x, "feature_name"], feature.url = fixed.continuous[x, "file_path"], nbins = fixed.continuous[x, "nbins"], genome.build = genome.build)
          
          write.table(fixed.continuous.binned, file = paste(feature.dir, fixed.continuous[x, "feature_name"], ".bed", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
          
          fixed.continuous.urls[x, "file_path"] <- paste(feature.dir, fixed.continuous[x, "feature_name"], ".bed", sep = "")
          
        }
        
      }
      
      fixed.continuous.urls = fixed.continuous.urls[ ,c("feature_name", "file_path")]
      rownames(fixed.continuous.urls) = NULL
      if (!is.null(features.continuous.url.snv)) {
        
        rownames(features.continuous.url.snv) = NULL
        features.continuous.url.snv = rbind(features.continuous.url.snv, fixed.continuous.urls)
        
      } else {
        
        features.continuous.url.snv = fixed.continuous.urls
        
      }

    } 
    
    if (sum(all.fixed$feature_type == 0) > 0) {
      
      fixed.discrete = all.fixed[which(all.fixed$feature_type == 0), 1:2]
      rownames(fixed.discrete) = NULL
      if (!is.null(features.discrete.url.snv)) {
        
        rownames(features.discrete.url.snv) = NULL
        features.discrete.url.snv = rbind(features.discrete.url.snv, fixed.discrete)
        
      } else {
        
        features.discrete.url.snv = fixed.discrete
        
      }

    }
    
  } else if (!is.null(genomic.features.fixed.indel)) {
    
    all.fixed = read.delim(genomic.features.fixed.indel, stringsAsFactors = FALSE, header = TRUE)    
    if (sum(all.fixed$feature_type == 1) > 0) {
      
      fixed.continuous = all.fixed[which(all.fixed$feature_type == 1), ]
      fixed.continuous.urls = fixed.continuous  
      
      # Bin features 
      if (sum(!is.na(fixed.continuous$nbins)) > 0) {
        
        to.bin = which(!is.na(fixed.continuous$nbins))
        
        # Bin continuous features based on number of bins provided
        for (x in to.bin) {
          
          print(paste("Binning ", fixed.continuous[x, "feature_name"], " into ", fixed.continuous[x, "nbins"], " bins", sep = ""))
          
          fixed.continuous.binned = bin.continuous(feature.name = fixed.continuous[x, "feature_name"], feature.url = fixed.continuous[x, "file_path"], nbins = fixed.continuous[x, "nbins"], genome.build = genome.build)
          
          write.table(fixed.continuous.binned, file = paste(feature.dir, fixed.continuous[x, "feature_name"], ".bed", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
          
          fixed.continuous.urls[x, "file_path"] <- paste(feature.dir, fixed.continuous[x, "feature_name"], ".bed", sep = "")
          
        }
        
      }
      
      fixed.continuous.urls = fixed.continuous.urls[ ,c("feature_name", "file_path")]
      rownames(fixed.continuous.urls) = NULL
      if (!is.null(features.continuous.url.indel)) {
        
        rownames(features.continuous.url.indel) = NULL
        features.continuous.url.indel = rbind(features.continuous.url.indel, fixed.continuous.urls)
        
      } else {
        
        features.continuous.url.indel = fixed.continuous.urls
        
      }
      
    } 
    
    if (sum(all.fixed$feature_type == 0) > 0) {
      
      fixed.discrete = all.fixed[which(all.fixed$feature_type == 0), 1:2]
      rownames(fixed.discrete) = NULL
      if (!is.null(features.discrete.url.indel)) {
        
        rownames(features.discrete.url.indel) = NULL
        features.discrete.url.indel = rbind(features.discrete.url.indel, fixed.discrete)
        
      } else {
        
        features.discrete.url.indel = fixed.discrete
        
      }
      
    }
    
  }
  
return(list(freq.snv, feat.coef.snv, features.continuous.url.snv, features.discrete.url.snv, freq.indel, feat.coef.indel, features.continuous.url.indel, features.discrete.url.indel, sds.snv, sds.indel))
  
}
