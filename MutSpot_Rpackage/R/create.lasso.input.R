#' Create covariates and response matrices as input for LASSO.
#' 
#' @param sampled.sites.file Sampled mutated and non-mutated sites RDS file.
#' @param genomic.features.file Text file containing URLs of potential continuous and discrete epigenetic features to select from.
#' @param mutation.type Type of mutation, snv or indel.
#' @param cores Number of cores, default = 1.
#' @return List containing covariate matrix, response matrix, final continuous feature urls, final discrete feature urls.
#' @export

create.lasso.input = function(sampled.sites.file, genomic.features.file, mutation.type, cores = 1) {
  
  # Define sampled mutated and non-mutated sites 
  sites = readRDS(sampled.sites.file)
  
  # Define all potential epigenetic features
  all.features = read.delim(genomic.features.file, stringsAsFactors = FALSE, header = TRUE)
  
  # Define continuous epigenetic features, if exists
  if (sum(all.features$feature_type == 1) > 0) {
    
  continuous.features = all.features[which(all.features$feature_type == 1), ]
  continuous.features.urls=continuous.features  
  
  # Bin features 
  if (sum(!is.na(continuous.features$nbins)) > 0) {
  
    to.bin=which(!is.na(continuous.features$nbins))
    
  # Bin continuous features based on number of bins provided
  for (x in to.bin) {
    
    print(paste("Binning ", continuous.features[x, "feature_name"], " into ", continuous.features[x, "nbins"], " bins", sep=""))
    
    continuous.features.binned = bin.continuous(feature.name = continuous.features[x, "feature_name"], feature.url = continuous.features[x, "file_path"], nbins = continuous.features[x, "nbins"])

    feature.dir=gsub( gsub(".*./","",continuous.features[x,"file_path"]),"",continuous.features[x,"file_path"])
    write.table(continuous.features.binned, file = paste(feature.dir,continuous.features[x, "feature_name"],".bed",sep=""), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    
    continuous.features.urls[x,"file_path"] <- paste(feature.dir,continuous.features[x, "feature_name"],".bed",sep="")
    
  }

  }
  
  continuous.features.urls=continuous.features.urls[,c("feature_name","file_path")]
  
  # Add local mutation rate to continuous file
  continuous.features.urls = rbind(continuous.features.urls, c("local_mutrate", paste("./features/localmutrate_", mutation.type, ".bed", sep = "")))
  
  # Assign continuous epigenetic scores to each site
  selected.continuous.features.site1 = parallel::mclapply(1:nrow(continuous.features.urls), function(f) {
    
    print(continuous.features.urls[f, 2])
    df = read.delim(as.character(continuous.features.urls[f, 2]), stringsAsFactors = FALSE, header = FALSE)
    z = with(df, GenomicRanges::GRanges(V1, IRanges::IRanges(V2, V3), score = V4))
    sites.feat = sites
    sites.feat$score = 0
    sites.feat[S4Vectors::queryHits(suppressWarnings(IRanges::findOverlaps(sites.feat, z)))]$score = z[S4Vectors::subjectHits(suppressWarnings(IRanges::findOverlaps(sites.feat, z)))]$score
    if (length(unique(sites.feat$score)) == 1) {
      
      print(paste(continuous.features.urls[f, 1], "remove", sep = "  "))
      return(NULL)
      
    } else {
      
      mat = Matrix::Matrix(sites.feat$score, sparse = TRUE)
      colnames(mat) = continuous.features.urls[f, 1]
      return(mat)
      
    }
    
  }, mc.cores = cores)
  
  selected.continuous.features.site1 = do.call(cbind, selected.continuous.features.site1)
  
  } else {
    
    # Add local mutation rate as the only continuous feature
    continuous.features.urls = data.frame(V1 = "local_mutrate", V2 = paste("./features/localmutrate_", mutation.type, ".bed", sep = ""))
    continuous.features.urls$V1 = as.character(continuous.features.urls$V1)
    continuous.features.urls$V2 = as.character(continuous.features.urls$V2)
    
    # Assign continuous epigenetic scores to each site
    selected.continuous.features.site1 = parallel::mclapply(1:nrow(continuous.features.urls), function(f) {
      
      print(continuous.features.urls[f, 2])
      df = read.delim(as.character(continuous.features.urls[f, 2]), stringsAsFactors = FALSE, header = FALSE)
      z = with(df, GenomicRanges::GRanges(V1, IRanges::IRanges(V2, V3), score = V4))
      sites.feat = sites
      sites.feat$score = 0
      sites.feat[S4Vectors::queryHits(suppressWarnings(IRanges::findOverlaps(sites.feat, z)))]$score = z[S4Vectors::subjectHits(suppressWarnings(IRanges::findOverlaps(sites.feat, z)))]$score
      if (length(unique(sites.feat$score)) == 1) {
        
        print(paste(continuous.features.urls[f, 1], "remove", sep = "  "))
        return(NULL)
        
      } else {
        
        mat = Matrix::Matrix(sites.feat$score, sparse = TRUE)
        colnames(mat) = continuous.features.urls[f, 1]
        return(mat)
        
      }
      
    }, mc.cores = cores)
    
    selected.continuous.features.site1 = do.call(cbind, selected.continuous.features.site1)
    
  }
  
  # Define all potential discrete epigenetic features
  if (sum(all.features$feature_type == 0) > 0) {
    
    discrete.features = all.features[which(all.features$feature_type == 0), 1:2]

  # Assign discrete epigenetic scores to each site
  selected.discrete.features.site1 = parallel::mclapply(1:nrow(discrete.features), function(f) {
    
    print(discrete.features[f, 2])
    df = read.delim(as.character(discrete.features[f, 2]), stringsAsFactors = FALSE, header = FALSE)
    z = with(df, GenomicRanges::GRanges(V1, IRanges::IRanges(V2, V3)))
    sites.feat = sites
    sites.feat$score = "1"
    if (length(unique(S4Vectors::queryHits(suppressWarnings(IRanges::findOverlaps(sites.feat, z))))) == 0) {
      
      sites.feat$score = "0"
      
    } else if (length(unique(S4Vectors::queryHits(suppressWarnings(IRanges::findOverlaps(sites.feat, z))))) != length(sites.feat)) {
      
      sites.feat[-unique(S4Vectors::queryHits(suppressWarnings(IRanges::findOverlaps(sites.feat, z))))]$score = "0"
      
    }
    if (length(unique(sites.feat$score)) == 1) {
      
      print(paste(discrete.features[f, 1], "remove", sep = "  "))
      return(NULL)
      
    } else {
      
      mat = Matrix::Matrix(as.numeric(sites.feat$score), sparse = TRUE)
      colnames(mat) = discrete.features[f, 1]
      return(mat)
      
    }
    
  }, mc.cores = cores)
  
  selected.discrete.features.site1 = do.call(cbind, selected.discrete.features.site1)
  
  } else {
    
    discrete.features = NULL
    selected.discrete.features.site1 = NULL
    
  }
  
  # Covariate matrix
  t = cbind(selected.continuous.features.site1, selected.discrete.features.site1)
  # Response matrix
  y = factor(sites$mut)
  
  return(list(t, y, continuous.features.urls, discrete.features))
  
}
