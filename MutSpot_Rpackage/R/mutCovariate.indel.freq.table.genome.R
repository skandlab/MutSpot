#' Prepare covariates matrix for all sites in whole genome/specified region.
#'
#' @param continuous.features Continuous epigenetic features selected for model fitting.
#' @param discrete.features Discrete epigenetic features selected for model fitting.
#' @param polyAs All polyA positions in whole genome.
#' @param polyTs All polyT positions in whole genome.
#' @param polyCs All polyC positions in whole genome.
#' @param polyGs All polyG positions in whole genome.
#' @param sites All sites in whole genome/specified region.
#' @return Covariates matrix for all sites in whole genome/specified region.
#' @export

mutCovariate.indel.freq.table.genome = function(continuous.features, discrete.features, polyAs, polyTs, polyCs, polyGs, sites) {
  
  print (as.character(GenomicRanges::seqnames(sites[1])))
  
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
      
      GenomicRanges::values(sites)[IRanges::queryHits(ovl), i] = continuous.features[[i]][IRanges::subjectHits(ovl)]$score
      
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
      
      GenomicRanges::values(sites)[IRanges::queryHits(ovl), i] = 1
      
      }
    rm(ovl)
    rm(df)
    
  }
  features = c(features, names(discrete.features))
  
  }
  
  # Assign polyATCG scores to each site
  sites = sites + 5
  
  sites$polyT = 0
  sites$polyA = 0
  sites$polyG = 0
  sites$polyC = 0
  
  ovl = IRanges::findOverlaps(polyTs, sites, type = "within")
  if (length(ovl) != 0) {
    
    sites[IRanges::subjectHits(ovl)]$polyT = 1
    
    }
  rm(ovl)
  
  ovl = IRanges::findOverlaps(polyAs, sites, type = "within")
  if (length(ovl) != 0) {
    
    sites[IRanges::subjectHits(ovl)]$polyA = 1
    
    }
  rm(ovl)
  
  ovl = IRanges::findOverlaps(polyGs, sites, type = "within")
  if (length(ovl) != 0) {
    
    sites[IRanges::subjectHits(ovl)]$polyG = 1
    
    }
  rm(ovl)
  
  ovl = IRanges::findOverlaps(polyCs, sites, type = "within")
  if (length(ovl) != 0) {
    
    sites[IRanges::subjectHits(ovl)]$polyC = 1
    
    }
  rm(ovl)
  
  features = c(features, "polyT", "polyA", "polyG", "polyC")
  names(GenomicRanges::mcols(sites)) = features
  
  sites = data.frame(GenomicRanges::mcols(sites), check.names = FALSE)
  sites = data.table::as.data.table(sites)
  sites2 = plyr::count(sites)
  colnames(sites2)[-ncol(sites2)] = colnames(sites)
  return(sites2)
  
}
