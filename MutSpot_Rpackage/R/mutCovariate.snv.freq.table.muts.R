#' Prepare covariates matrix for mutated sites in a specfic chromosome.
#'
#' @param continuous.features Continuous epigenetic features selected for model fitting.
#' @param discrete.features Discrete epigenetic features selected for model fitting.
#' @param precompute.motif.pos All selected nucleotide contexts' positions in a specific chromosome.
#' @param nucleotide.selected Nucleotide contexts selected for model fitting.
#' @param sample.specific.features Sample-specific features.
#' @param sites Mutated sites in a specific chromosome.
#' @return Covariate matrix for mutated sites.
#' @export

mutCovariate.snv.freq.table.muts = function(continuous.features, discrete.features, precompute.motif.pos, nucleotide.selected, sample.specific.features, sites) {
  
  print (as.character(GenomeInfoDb::seqnames(sites[1])))
  
  # Replace sid by sample mutation count
  # sites$sid = ind.mut.count[as.character(sites$sid)]
  features1 = features = names(GenomicRanges::mcols(sites))
  # features1 = features1[which(features1 != "sid")]
  
  # Assign sample specific feature scores to each site
  for (i in colnames(sample.specific.features)) {
    
    df = data.frame(numeric(length(sites)))
    colnames(df) = i
    GenomicRanges::values(sites) = data.frame(GenomicRanges::values(sites), df, check.names = FALSE)
    names(GenomicRanges::mcols(sites))[length(names(GenomicRanges::mcols(sites)))] = i
    GenomicRanges::values(sites)[ ,i] = sample.specific.features[as.character(sites$sid), i]
  
  }
  features = c(features, colnames(sample.specific.features))
  
  # Assign continuous epigenetic scores to each site
  if (!is.null(continuous.features)) {
    
    for (i in names(continuous.features)) {
      
      df = data.frame(numeric(length(sites)))
      colnames(df) = i
      GenomicRanges::values(sites) = data.frame(GenomicRanges::values(sites), df, check.names = FALSE)
      names(GenomicRanges::mcols(sites))[length(names(GenomicRanges::mcols(sites)))] = i
      ovl = IRanges::findOverlaps(sites, continuous.features[[i]])
      if (length(ovl) != 0) {
        
        GenomicRanges::values(sites)[S4Vectors::queryHits(ovl), i] = continuous.features[[i]][S4Vectors::subjectHits(ovl)]$score
        
        }
      rm(ovl)
      rm(df)
      
    }
    features = c(features, names(continuous.features))
    
    }
  
  # Assign discrete epigenetic scores to each site
  if (!is.null(discrete.features)) {
    
    for (i in names(discrete.features)) {
      
      df = data.frame(numeric(length(sites)))
      colnames(df) = i
      GenomicRanges::values(sites) = data.frame(GenomicRanges::values(sites), df, check.names = FALSE)
      names(GenomicRanges::mcols(sites))[length(names(GenomicRanges::mcols(sites)))] = i
      ovl = IRanges::findOverlaps(sites, discrete.features[[i]])
      
      if (length(ovl) != 0) {
        
        GenomicRanges::values(sites)[S4Vectors::queryHits(ovl), i] = 1
        
        }
      rm(ovl)
      rm(df)
      
    }
    features = c(features, names(discrete.features))
    
    }
  
  # Assign nucleotide contexts scores to each site
  if (!is.null(nucleotide.selected)) {
    
    feat = NULL
    
    for (j in unique(nucleotide.selected$type)) {

      # print(j)
      for (k in nucleotide.selected[which(nucleotide.selected$type == j), "sequence"]) {
        
        # print(k)
        if (j == "oneMer") {
          
          df = data.frame(numeric(length(sites)))
          mot = paste(j, k, sep = "")
          colnames(df) = mot
          GenomicRanges::values(sites) = data.frame(GenomicRanges::values(sites), df, check.names = FALSE)
          GenomicRanges::values(sites)[which(GenomicRanges::start(sites) %in% precompute.motif.pos[[mot]]), mot] = 1
          rm(df)
          
        } else if (j == "threeMer") {
          
          df = data.frame(numeric(length(sites)))
          mot = paste(j, k, sep = "")
          colnames(df) = mot
          GenomicRanges::values(sites) = data.frame(GenomicRanges::values(sites), df, check.names = FALSE)
          GenomicRanges::values(sites)[which((GenomicRanges::start(sites) - 1) %in% precompute.motif.pos[[mot]]), mot] = 1
          rm(df)
          
        } else if (j == "threeRight") {
          
          df = data.frame(numeric(length(sites)))
          mot = paste(j, k, sep = "")
          colnames(df) = mot
          GenomicRanges::values(sites) = data.frame(GenomicRanges::values(sites), df, check.names = FALSE)
          GenomicRanges::values(sites)[which(GenomicRanges::start(sites) %in% precompute.motif.pos[[mot]]), mot] = 1
          rm(df)
          
        } else if (j == "threeLeft") {
          
          df = data.frame(numeric(length(sites)))
          mot = paste(j, k, sep = "")
          colnames(df) = mot
          GenomicRanges::values(sites) = data.frame(GenomicRanges::values(sites), df, check.names = FALSE)
          GenomicRanges::values(sites)[which(GenomicRanges::start(sites) %in% precompute.motif.pos[[mot]]), mot] = 1
          rm(df)
          
        } else if (j == "fiveMer") {
          
          df = data.frame(numeric(length(sites)))
          mot = paste(j, k, sep = "")
          colnames(df) = mot
          GenomicRanges::values(sites) = data.frame(GenomicRanges::values(sites), df, check.names = FALSE)
          GenomicRanges::values(sites)[which((GenomicRanges::start(sites) - 2) %in% precompute.motif.pos[[mot]]), mot] = 1
          rm(df)
          
        } else if (j == "fiveRight") {
          
          df = data.frame(numeric(length(sites)))
          mot = paste(j, k, sep = "")
          colnames(df) = mot
          GenomicRanges::values(sites) = data.frame(GenomicRanges::values(sites), df, check.names = FALSE)
          GenomicRanges::values(sites)[which(GenomicRanges::start(sites) %in% precompute.motif.pos[[mot]]), mot] = 1
          rm(df)
          
        } else {
          
          df = data.frame(numeric(length(sites)))
          mot = paste(j, k, sep="")
          colnames(df) = mot
          GenomicRanges::values(sites) = data.frame(GenomicRanges::values(sites), df, check.names = FALSE)
          GenomicRanges::values(sites)[which(GenomicRanges::start(sites) %in% precompute.motif.pos[[mot]]), mot] = 1
          rm(df)
          
        }
        feat = c(feat, paste(j, k, sep = ""))
        
      }
      
    }
    features = c(features, feat)
    
  }
  
  names(GenomicRanges::mcols(sites)) = features
  
  sites = data.frame(GenomicRanges::mcols(sites), check.names = FALSE)
  if (length(features1) != 0) {
    
    sites = sites[ ,-which(colnames(sites) %in% features1)]
    
  }
  sites = data.table::as.data.table(sites)
  sites2 = plyr::count(sites)
  colnames(sites2)[-ncol(sites2)] = colnames(sites)
  
  return(sites2)
  
}

