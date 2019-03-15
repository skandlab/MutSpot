#' Prepare covariates matrix for all sites in whole genome/specified region.
#'
#' @param continuous.features Continuous epigenetic features selected for model fitting.
#' @param discrete.features Discrete epigenetic features selected for model fitting.
#' @param precompute.motif.pos All selected nucleotide contexts' positions in a specific chromosome.
#' @param nucleotide.selected Nucleotide contexts selected for model fitting.
#' @param sites All sites in whole genome/specified region.
#' @return Covariates matrix for all sites in whole genome/specified region.
#' @export

mutCovariate.snv.freq.table.genome = function(continuous.features, discrete.features, precompute.motif.pos, nucleotide.selected, sites) {
  
  print(as.character(GenomicRanges::seqnames(sites[1])))
  
  sites = IRanges::tile(sites, width = 1)
  sites = BiocGenerics::unlist(sites)
  features = names(GenomicRanges::mcols(sites))
  
  # Assign continuous epigenetic scores to each site
  if (!is.null(continuous.features)) {
    
    for (i in names(continuous.features)) {
      
      df = data.frame(numeric(length(sites)))
      colnames(df) = i
      GenomicRanges::values(sites) = data.frame(GenomicRanges::values(sites), df, check.names = FALSE)
      names(GenomicRanges::mcols(sites))[length(names(GenomicRanges::mcols(sites)))] = i
      ovl = IRanges::findOverlaps(sites, continuous.features[[i]])
      if (length(ovl) != 0) {
        
        GenomicRanges::values(sites)[S4Vectors::queryHits(ovl), i] = continuous.features[[i]][S4Vectors::subjectHits(ovl)]$score}
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
      
    }
    features = c(features, names(discrete.features))
    
    }
  
  # Assign nucleotide contexts scores to each site
  if (!is.null(nucleotide.selected)) {
    
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
          mot = paste(j, k, sep = "")
          colnames(df) = mot
          GenomicRanges::values(sites) = data.frame(GenomicRanges::values(sites), df, check.names = FALSE)
          GenomicRanges::values(sites)[which(GenomicRanges::start(sites) %in% precompute.motif.pos[[mot]]), mot] = 1
          rm(df)
          
        }
        
      }
      
    }
    features = c(features, paste(nucleotide.selected$type, nucleotide.selected$sequence, sep = ""))
    
    }
  
  names(GenomicRanges::mcols(sites)) = features
  
  sites = data.frame(GenomicRanges::mcols(sites), check.names = FALSE)
  sites = data.table::as.data.table(sites)
  sites2 = plyr::count(sites)
  colnames(sites2)[-ncol(sites2)] = colnames(sites)
  
  return(sites2)

}
