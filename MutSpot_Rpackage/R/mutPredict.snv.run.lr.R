#' Run mutation recurrence for SNV.
#' 
#' @param roi GRanges containing all hotspots.
#' @param maf.snv GRanges containing SNV mutations in region of interest.
#' @param maf.snv2 GRanges containing all SNV mutations.
#' @param model.snv SNV model.
#' @param continuous.features.snv Full scores list of SNV continuous features selected for model.
#' @param discrete.features.snv Full scores list of SNV discrete features selected for model.
#' @param motifs Selected nucleotide contexts to be extracted.
#' @param nucleotide.selected Selected nucleotide contexts table.
#' @param sample.specific.features Sample-specific features.
#' @param continuous.sample.specific Continuous features in the final model.
#' @param min.count Minimum number of mutated samples in each hotspot, default = 2.
#' @param genome.size Total number of hotspots to run analysis on, default = 2533374732.
#' @param cores Number of cores, default = 1.
#' @param only.samples Restrict analysis to certain samples, default = NULL.
#' @param debug To print errors, default = FALSE.
#' @return Dataframe containing predicted hotspots significance.
#' @export

mutPredict.snv.run.lr <- function(roi, maf.snv, maf.snv2, model.snv, continuous.features.snv, discrete.features.snv, motifs, nucleotide.selected, sample.specific.features, continuous.sample.specific, min.count = 2, genome.size = 2533374732, cores = 1, only.samples = NULL, debug = FALSE) {
  
  # Check that roi names are unique
  if (any(duplicated(names(roi)))) {
    
    print("Error: regions with duplcated names")
    return()
    
  }
  
  # Print number of roi
  print(paste("Using ", length(roi), " ROIs from input", sep = ""))
  
  if (debug) {
    
    print(roi)
    
  }
  
  # Filter maf if sample list specified
  if (!is.null(only.samples)) {
    
    print("Filtering MAF based on supplied sample list")
    maf.snv = maf.snv[which(maf.snv$sid %in% only.samples)]
    
  }
  
  # Find overlaps between roi and maf to identify mutations in each roi
  print(">> Intersecting ROIs and MAF ...")
  maf.ovl <- IRanges::findOverlaps(roi, maf.snv, ignore.strand = TRUE)
  maf.ovl.m = IRanges::as.matrix(maf.ovl)
  
  if (debug) {
    
    print(maf.ovl.m)
    
  }
  
  # Select rois mutated in min.count.sample i.e. look at rois with at least min.count number of samples mutated
  nsamples = length(unique(maf.snv$sid))
  # Encode sample names as integers to increase the efficiency of the code
  roi.count.mutated = tapply(maf.snv$sid[maf.ovl.m[ ,2]], names(roi)[maf.ovl.m[ ,1]],function(s) length(unique(s)))
  roi.mutated = names(roi.count.mutated)[roi.count.mutated >= min.count]
  print(paste("Using ", length(roi.mutated), " ROIs mutated in >=", min.count, " samples", sep = ""))
  
  # Run analysis on all mutated roi with min.count
  rois.to.run = roi.mutated
  
  if (debug) {
    
    print(rois.to.run)
    
  }
  
  results = c()
  
  if (length(rois.to.run) == 0) {
    
    return(results)
    
  }
  
  process.single.roi <- function(x) {
    
    # Print progress for every 1% of data
    roi.progress <- which(rois.to.run == x)
    if ((roi.progress %% ceiling(length(rois.to.run) / 100)) == 0) {
      
      print(paste("Progress : ", roi.progress, "/", length(rois.to.run), sep = ""))
      
    }
    
    # Make sure we are only using first hit
    # Should not have multiple hits in new version, where we checked for duplicated roi names
    x.idx = which(names(roi) == x)[1]
    
    # Extract features for all sites in roi, note roi here is a GRange object, not GRangeList
    # roi.feat.snv = mutPredict.snv.get.features(BiocGenerics::unlist(GenomicRanges::reduce(roi[x.idx])), continuous.features.snv, discrete.features.snv, motifs)
    roi.feat.snv = mutPredict.snv.get.features(GenomicRanges::reduce(roi[x.idx]), continuous.features.snv, discrete.features.snv, motifs)
    
    x.len.snv = nrow(roi.feat.snv)
    
    # Vector of patient IDs
    sid = as.character(unique(maf.snv2$sid))
    
    # Compute probability of mutation in region for all samples
    p.bg = sapply (sid, function(s) {
      
      # Encode sid as the mutation count of each individual    
      sid = rep(s, x.len.snv)
      df = roi.feat.snv
      if (!is.null(nucleotide.selected)) {
        
        for (t in 1:nrow(nucleotide.selected)) {
          
          df[ ,paste(nucleotide.selected[t, 2], nucleotide.selected[t, 1], sep = "")] = ifelse(df[ ,nucleotide.selected[t, 2]] == nucleotide.selected[t, 1], "1", "0")
          
        }
        df = df[ ,!colnames(df) %in% unique(nucleotide.selected[ ,2])]
        
        }
      
      # Assign sample specific feature scores to each site
      for (i in colnames(sample.specific.features)) {
        
        df[ ,i] = sample.specific.features[as.character(s), i]
        
      }
      
      df = df[ ,c(colnames(sample.specific.features), colnames(df)[!colnames(df) %in% names(sample.specific.features)])]

      if (sum(!colnames(df) %in% c(continuous.sample.specific,"sample.count")) > 0 ) {
        
        for (i in colnames(df)[!colnames(df) %in% c(continuous.sample.specific, "sample.count")]) {
          
            df[ ,i] <- as.character(df[ ,i])
            
          }
      }
      
      # Compute background mutation rate foreach site in each individual
      if (!"glmnet" %in% class(model.snv)) {
        
      p = stats::predict(object = model.snv, newdata = df, type = "response") 
      
      } else {
        
        p = stats::predict(object = model.snv, newx = data.matrix(df), type = "response")
        
      }
      
      # Compute background mutation rate of the region in each individual
      p.roi.snv = 1 - prod(1 - p)   
      
      return(p.roi.snv)
      
      } )
    
    # Compute k, number of samples where x is mutated
    if (debug) {
      
      print(x)
      print(x.idx)
      print(maf.ovl.m[ ,1] == x.idx)
      
    }
    
    k = length(unique(maf.snv$sid[maf.ovl.m[maf.ovl.m[ ,1] == x.idx, 2]]))
    pval = poibin::ppoibin(length(p.bg) - k, 1 - p.bg, method = 'RF')
    c(x, pval, x.len.snv, mean(p.bg), k)
    
  }
  
  results = parallel::mclapply(rois.to.run, function(xr) { process.single.roi(xr) }, mc.cores = cores, mc.preschedule = FALSE, mc.silent = FALSE)
  
  # Error checking and remove entries with error
  results.error = sapply(results, function(x) { length(x) == 1 } )
  print(rois.to.run[which(results.error)])
  print(results[results.error])
  results = results[!results.error]
  
  results = do.call(rbind,results)
  rownames(results) = results[ ,1]
  # drop = FALSE --> forces matrix if single row / vector
  results = results[ ,-1, drop = FALSE]
  class(results) <- "numeric"
  # Compute BH fdr
  fdr = p.adjust(results[ ,1], method = 'BH') * genome.size / nrow(results)
  fdr = sapply(fdr, function(x) min(x, 1))  
  results = cbind(results, fdr)
  colnames(results) = c('pval', 'length', 'p.bg', 'k', 'fdr')
  results = results[order(results[ ,'pval']), , drop = FALSE]
  
  return(results)
  
}
