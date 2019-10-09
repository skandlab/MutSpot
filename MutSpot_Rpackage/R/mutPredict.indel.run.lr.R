#' Run mutation recurrence for indel.
#' 
#' @param output.dir Save plot in given output directory.
#' @param merge.hotspots To plot overlapping hotspots as 1 hotspot or individual hotspots, default = TRUE.
#' @param indel.mutations.file Indel mutations found in region of interest MAF file.
#' @param fdr.cutoff FDR cutoff, default = 0.1.
#' @param color.line Color given FDR cutoff, default = red.
#' @param color.dots Color hotspots that passed given FDR cutoff, default = maroon1.
#' @param color.muts Color points, default = orange.
#' @param top.no Number of top hotspots to plot, default = 3.
#' @param promoter.file Promoter regions bed file, default file = Ensembl75.promoters.coding.bed.
#' @param utr3.file 3'UTR regions bed file, default file = Ensembl75.3UTR.coding.bed.
#' @param utr5.file 5'UTR regions bed file, default file = Ensembl75.5UTR.coding.bed.
#' @param other.annotations Text file containing URLs of additional regions to be annotated, default = NULL.
#' @param debug To delete temporary files or not, default = FALSE.
#' @param cores Number of cores, default = 1.
#' @return Dataframe containing predicted hotspots significance.
#' @export

mutPredict.indel.run.lr <- function(output.dir, merge.hotspots = TRUE, indel.mutations.file, fdr.cutoff = 0.1, color.line = "red", color.dots = "maroon1", color.muts = "orange", top.no = 3,
                                    promoter.file = system.file("extdata", "Ensembl75.promoters.coding.bed", package = "MutSpot"),
                                    utr3.file = system.file("extdata", "Ensembl75.3UTR.coding.bed", package = "MutSpot"), 
                                    utr5.file = system.file("extdata", "Ensembl75.5UTR.coding.bed", package = "MutSpot"), 
                                    other.annotations = NULL, debug = FALSE, cores = 1) {
  
  if (!"temp_1.RDS" %in% list.files(output.dir)) {
    
    print("No hotspots found")
    mut.rec.hotspot = mut.rec.hotspot2 = ann.results.indel = NULL
    
  } else {
    
    mut.regions = readRDS(paste(output.dir, "temp_1.RDS", sep = ""))
    maf.indel = readRDS(paste(output.dir, "temp_2.RDS", sep = ""))
    maf.ovl.m = readRDS(paste(output.dir, "temp_3.RDS", sep = ""))
    continuous.selected.features.indel = readRDS(paste(output.dir, "temp_4.RDS", sep = ""))
    discrete.selected.features.indel = readRDS(paste(output.dir, "temp_5.RDS", sep = ""))
    # LRmodel.indel = readRDS(paste(output.dir,"temp_6.RDS",sep=""))
    load(paste(output.dir, "indel-LRmodel", sep = ""))
    rois.to.run = readRDS(paste(output.dir,"temp_7.RDS",sep=""))
    sample.specific.features = readRDS(paste(output.dir,"temp_8.RDS",sep=""))
    continuous.sample.specific = readRDS(paste(output.dir,"temp_9.RDS",sep=""))
    genome.size = readRDS(paste(output.dir,"temp_10.RDS",sep=""))
    
    # Remove temporary files if not debugging
    if (!debug) {
      
      delete.files = Sys.glob(paste(output.dir, "temp_*.RDS", sep = ""))
      for (i in delete.files) {
        
        unlink(i)
        
      }
      
    }
    
  process.single.roi <- function(x) {
    
    # Print progress for every 1% of data
    roi.progress <- which(rois.to.run == x)
    if ((roi.progress %% ceiling(length(rois.to.run) / 100)) == 0) {
      
      print(paste("Progress : ", roi.progress, "/", length(rois.to.run), sep = ""))
      
    }
    
    # Make sure we are only using first hit
    # Should not have multiple hits in new version, where we checked for duplicated roi names
    x.idx = which(names(mut.regions) == x)[1]
    
    # Extract features for all sites in roi, note roi here is a GRange object, not GRangeList
    roi.feat.indel = mutPredict.indel.get.features(GenomicRanges::reduce(mut.regions[x.idx]), continuous.selected.features.indel, discrete.selected.features.indel)
    
    # If poly A or poly G are no longer in model
    if (!"glmnet" %in% class(LRmodel)) {
      
    if (!"polyA1" %in% names(LRmodel$coefficients)) {
      
      roi.feat.indel = roi.feat.indel[ ,-which(colnames(roi.feat.indel) == "polyA")]
      
    }
    if (!"polyG1" %in% names(LRmodel$coefficients)) {
      
      roi.feat.indel = roi.feat.indel[ ,-which(colnames(roi.feat.indel) == "polyG")]
      
    }
      
    } else {
      
      if (!"polyA1" %in% rownames(LRmodel$beta)) {
        
        roi.feat.indel = roi.feat.indel[ ,-which(colnames(roi.feat.indel) == "polyA")]
        
      }
      if (!"polyG1" %in% rownames(LRmodel$beta)) {
        
        roi.feat.indel = roi.feat.indel[ ,-which(colnames(roi.feat.indel) == "polyG")]
        
      }
      
    }
    z=colnames(roi.feat.indel)[1]
    x.len.indel = nrow(roi.feat.indel)
    
    # Vector of patient IDs
    sid = as.character(unique(rownames(sample.specific.features)))
    nind = length(sid)
    
    roi.feat.indel = roi.feat.indel[rep(1:nrow(roi.feat.indel), nind), ]
    if (class(roi.feat.indel)=="numeric") {
      
      roi.feat.indel=as.data.frame(roi.feat.indel)
      colnames(roi.feat.indel)[1]=z
      
    }
    
    roi.feat.indel$sid = rep(sid, each = x.len.indel)
    rm(sid)
    
    # Assign sample specific feature scores to each site
    for (i in colnames(sample.specific.features)) {
      
      roi.feat.indel[ ,i] = sample.specific.features[as.character(roi.feat.indel$sid), i]
      
    }
    
    roi.feat.indel = roi.feat.indel[ ,c(colnames(sample.specific.features), colnames(roi.feat.indel)[!colnames(roi.feat.indel) %in% names(sample.specific.features)])]
    
    if (sum(!colnames(roi.feat.indel) %in% c(continuous.sample.specific, "ind.mut.count")) > 0 ) {
      
      for (i in colnames(roi.feat.indel)[!colnames(roi.feat.indel) %in% c(continuous.sample.specific, "ind.mut.count")]) {
        
        roi.feat.indel[ ,i] <- as.character(roi.feat.indel[ ,i])
        
      }
      
    }
    
    # Compute background mutation rate foreach site in each individual
    if (!"glmnet" %in% class(LRmodel)) {
      
      p = stats::predict(object = LRmodel, newdata = roi.feat.indel, type = "response") 
      
    } else {
      
      p = stats::predict(object = LRmodel, newx = data.matrix(roi.feat.indel), type = "response")
      
    }
    
    # Compute background mutation rate of the region in each individual
    p.bg = sapply(seq(1,nrow(roi.feat.indel),nrow(roi.feat.indel)/nind),FUN=function(pp) 1 - prod(1 - p[pp:(pp+(nrow(roi.feat.indel)/nind)-1)]))
    
    k = length(unique(maf.indel$sid[maf.ovl.m[maf.ovl.m[ ,1] == x.idx, 2]]))
    pval = poibin::ppoibin(length(p.bg) - k, 1 - p.bg, method = 'RF')
    
    rm(roi.feat.indel)
    gc()
    
    return(c(x, pval, x.len.indel, mean(p.bg), k))
    
  }
  
  if (length(rois.to.run) == 0) {
    
    print("No hotspots found")
    mut.rec.hotspot = mut.rec.hotspot2 = ann.results.indel = NULL
    
  } else {
    
  results = parallel::mclapply(rois.to.run, function(xr) { process.single.roi(xr) }, mc.cores = cores, mc.preschedule = FALSE, mc.silent = FALSE)
  
  # Error checking and remove entries with error
  results.error = sapply(results, function(x) { length(x) == 1 })
  print(rois.to.run[which(results.error)])
  print(results[results.error])
  results = results[!results.error]
  
  results = do.call(rbind, results)
  rownames(results) = results[ ,1]
  # drop = FALSE, forces matrix if single row / vector
  results = results[ ,-1, drop = FALSE]
  class(results) <- "numeric"
  # Compute BH fdr
  fdr = stats::p.adjust(results[ ,1], method = 'BH') * genome.size / nrow(results)
  fdr = sapply(fdr, function(x) min(x, 1))  
  results = cbind(results, fdr)
  colnames(results) = c('pval', 'length', 'p.bg', 'k', 'fdr')
  results = results[order(results[ ,'pval']), , drop = FALSE]
  
  mut.rec.hotspot = data.frame(chrom = as.character(GenomeInfoDb::seqnames(mut.regions[rownames(results)])), start = GenomicRanges::start(mut.regions[rownames(results)]), end = GenomicRanges::end(mut.regions[rownames(results)]), results)

  # merge overlapping hotspots in mut.rec.hotspot, assign smallest pvalue and recalculate k
  if (nrow(mut.rec.hotspot) != 0) {
    
    if (merge.hotspots) {
      
      print("Merge overlapping hotspots in final results...")
      mut.rec.hotspot2 = mut.rec.hotspot
      mut.rec.hotspot2$ID = as.character(rownames(mut.rec.hotspot2))
      mut.rec.hotspot2 = with(mut.rec.hotspot2, GenomicRanges::GRanges(chrom, IRanges::IRanges(start, end), pval = pval, length = length, p.bg = p.bg, k = k, fdr = fdr, ID = ID))
      hotspots = IRanges::reduce(mut.rec.hotspot2)
      hotspots$hs = paste("hs", 1:length(hotspots), sep = "")
      ovl.mut = IRanges::findOverlaps(maf.indel, hotspots)
      hotspots2 = hotspots[S4Vectors::subjectHits(ovl.mut)]
      hotspots2$sample = maf.indel[S4Vectors::queryHits(ovl.mut)]$sid
      hotspots2 = GenomicRanges::as.data.frame(hotspots2)
      hotspots2 = aggregate(sample ~ hs, hotspots2, FUN = function(k) length(unique(k)))
      colnames(hotspots2)[2] = "k"
      rownames(hotspots2) = hotspots2$hs
      hotspots$k = 0
      for (i in 1:length(hotspots)) {
        # print(i)
        hotspots$k[i] = hotspots2[which(hotspots2$hs == hotspots$hs[i]), "k"]
        
      }
      
      ovl = IRanges::findOverlaps(mut.rec.hotspot2, hotspots)
      mut.rec.hotspot2 = mut.rec.hotspot2[S4Vectors::queryHits(ovl)]
      mut.rec.hotspot2$hs = hotspots[S4Vectors::subjectHits(ovl)]$hs
      mut.rec.hotspot2$region.start = IRanges::start(hotspots[S4Vectors::subjectHits(ovl)])
      mut.rec.hotspot2$region.end = IRanges::end(hotspots[S4Vectors::subjectHits(ovl)])
      mut.rec.hotspot2$new.k = hotspots[S4Vectors::subjectHits(ovl)]$k
      mut.rec.hotspot2 = GenomicRanges::as.data.frame(mut.rec.hotspot2)
      mut.rec.hotspot2 = mut.rec.hotspot2[order(mut.rec.hotspot2$pval, decreasing = FALSE), ]
      mut.rec.hotspot2 = mut.rec.hotspot2[!duplicated(mut.rec.hotspot2$hs), ]
      mut.rec.hotspot2 = mut.rec.hotspot2[ ,c("seqnames", "region.start", "region.end", "pval", "length", "p.bg", "new.k", "fdr", "ID")]
      mut.rec.hotspot2$length = mut.rec.hotspot2$region.end - mut.rec.hotspot2$region.start + 1
      colnames(mut.rec.hotspot2) = c(colnames(mut.rec.hotspot), "ID")
      rownames(mut.rec.hotspot2) = mut.rec.hotspot2$ID
      mut.rec.hotspot2 = mut.rec.hotspot2[ ,-ncol(mut.rec.hotspot2)]
      
    } else {
      
      mut.rec.hotspot2 = NULL
      
    }
    
  } else {
    
    mut.rec.hotspot2 = NULL
    
  }
  
  # Plot hotspots manhattan
  manhattan.indel = paste(output.dir, "indel_manhattan.pdf", sep = "")
  
  if (!is.null(mut.rec.hotspot2)) {
    
    print("Plot manhattan figure for indel")
    pdf(manhattan.indel)
    plot_manhattan(hotspots.file = mut.rec.hotspot2, fdr.cutoff = fdr.cutoff, color.line = color.line, color.dots = color.dots)
    dev.off()
    
    print("Plot top hotspots for snv/indel")
    plot_top_hits(hotspots.file = mut.rec.hotspot2, fdr.cutoff = fdr.cutoff, color.muts = color.muts, mutations.file = indel.mutations.file, mutation.type = "indel", top.no = top.no, output.dir = output.dir)
    
    # Annotate hotspots
    ann.results.indel = mutAnnotate(hotspots.file = mut.rec.hotspot2, promoter.file = promoter.file,
                                    utr3.file = utr3.file, utr5.file = utr5.file,
                                    other.annotations = other.annotations)
    
  } else if (nrow(mut.rec.hotspot) != 0) {
    
    print("Plot manhattan figure for indel")
    pdf(manhattan.indel)
    plot_manhattan(hotspots.file = mut.rec.hotspot, fdr.cutoff = fdr.cutoff, color.line = color.line, color.dots = color.dots)
    dev.off()
    
    print("Plot top hotspots for indel")
    plot_top_hits(hotspots.file = mut.rec.hotspot, fdr.cutoff = fdr.cutoff, color.muts = color.muts, mutations.file = indel.mutations.file, mutation.type = "indel", top.no = top.no, output.dir = output.dir)
    
    # Annotate hotspots
    ann.results.indel = mutAnnotate(hotspots.file = mut.rec.hotspot, promoter.file = promoter.file,
                                    utr3.file = utr3.file, utr5.file = utr5.file,
                                    other.annotations = other.annotations)
    
  } else {
    
    ann.results.indel = NULL
    
  }
  
  }
  
  }
  
  return(list(mut.rec.hotspot, mut.rec.hotspot2, ann.results.indel))
  
}

