#' Run mutation recurrence for SNV.
#' 
#' @param output.dir Results directory.
#' @param merge.hotspots To plot overlapping hotspots as 1 hotspot or individual hotspots, default = TRUE.
#' @param snv.mutations.file SNV mutations found in region of interest MAF file.
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

mutPredict.snv.run.lr = function(output.dir, merge.hotspots = TRUE, snv.mutations.file, fdr.cutoff = 0.1, color.line = "red", color.dots = "maroon1", color.muts = "orange", top.no = 3,
                                 promoter.file = system.file("extdata", "Ensembl75.promoters.coding.bed", package = "MutSpot"),
                                 utr3.file = system.file("extdata", "Ensembl75.3UTR.coding.bed", package = "MutSpot"), 
                                 utr5.file = system.file("extdata", "Ensembl75.5UTR.coding.bed", package = "MutSpot"), 
                                 other.annotations = NULL, debug = FALSE, cores = 1) {
  
  if (!"temp-1.RDS" %in% list.files(output.dir)) {
    
    print("No hotspots found")
    mut.rec.hotspot = mut.rec.hotspot2 = ann.results.snv = NULL
    
  } else {
    
    mut.regions = readRDS(paste(output.dir, "temp-1.RDS", sep = ""))
    maf.snv = readRDS(paste(output.dir, "temp-2.RDS", sep = ""))
    maf.ovl.m = readRDS(paste(output.dir, "temp-3.RDS", sep = ""))
    continuous.selected.features.snv = readRDS(paste(output.dir, "temp-4.RDS", sep = ""))
    discrete.selected.features.snv = readRDS(paste(output.dir, "temp-5.RDS", sep = ""))
    sel.motif = readRDS(paste(output.dir, "temp-6.RDS",sep=""))
    # LRmodel.snv = readRDS(paste(output.dir,"temp-7.RDS",sep=""))
    load(paste(output.dir, "snv-LRmodel", sep = ""))
    nucleotide.selected = readRDS(paste(output.dir,"temp-8.RDS",sep=""))
    rois.to.run = readRDS(paste(output.dir,"temp-9.RDS",sep=""))
    sample.specific.features = readRDS(paste(output.dir,"temp-10.RDS",sep=""))
    continuous.sample.specific = readRDS(paste(output.dir,"temp-11.RDS",sep=""))
    genome.size = readRDS(paste(output.dir,"temp-12.RDS",sep=""))

    # Remove temporary files if not debugging
    if (!debug) {
      
      delete.files = Sys.glob(paste(output.dir, "temp-*.RDS", sep = ""))
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
      roi.feat.snv = mutPredict.snv.get.features(GenomicRanges::reduce(mut.regions[x.idx]), continuous.selected.features.snv, discrete.selected.features.snv, sel.motif)
      z=colnames(roi.feat.snv)[1]
      
      x.len.snv = nrow(roi.feat.snv)
      
      # Vector of patient IDs
      sid = as.character(unique(rownames(sample.specific.features)))
      nind = length(sid)
      
      if (!is.null(nucleotide.selected)) {
        
        for (t in 1:nrow(nucleotide.selected)) {
          
          roi.feat.snv[ ,paste(nucleotide.selected[t, 2], nucleotide.selected[t, 1], sep = "")] = ifelse(roi.feat.snv[ ,nucleotide.selected[t, 2]] == nucleotide.selected[t, 1], "1", "0")
          
        }
        roi.feat.snv = roi.feat.snv[ ,!colnames(roi.feat.snv) %in% unique(nucleotide.selected[ ,2])]
        
      }
      
      roi.feat.snv = roi.feat.snv[rep(1:nrow(roi.feat.snv), nind), ]
      if (class(roi.feat.snv)=="numeric") {
        
        roi.feat.snv=as.data.frame(roi.feat.snv)
        colnames(roi.feat.snv)[1]=z
        
      }
      
      roi.feat.snv$sid = rep(sid, each = x.len.snv)
      rm(sid)
      
      # Assign sample specific feature scores to each site
      for (i in colnames(sample.specific.features)) {
        
        roi.feat.snv[ ,i] = sample.specific.features[as.character(roi.feat.snv$sid), i]
        
      }
      
      roi.feat.snv = roi.feat.snv[ ,c(colnames(sample.specific.features), colnames(roi.feat.snv)[!colnames(roi.feat.snv) %in% names(sample.specific.features)])]
      
      if (sum(!colnames(roi.feat.snv) %in% c(continuous.sample.specific,"ind.mut.count")) > 0 ) {
        
        for (i in colnames(roi.feat.snv)[!colnames(roi.feat.snv) %in% c(continuous.sample.specific, "ind.mut.count")]) {
          
          roi.feat.snv[ ,i] <- as.character(roi.feat.snv[ ,i])
          
        }
      }
      
      # Compute background mutation rate foreach site in each individual
      if (!"glmnet" %in% class(LRmodel)) {
        
        p = stats::predict(object = LRmodel, newdata = roi.feat.snv, type = "response")
        
      } else {
        
        p = stats::predict(object = LRmodel, newx = data.matrix(roi.feat.snv), type = "response")
        
      }
      
      # Compute background mutation rate of the region in each individual
      p.bg = sapply(seq(1,nrow(roi.feat.snv),nrow(roi.feat.snv)/nind),FUN=function(pp) 1 - prod(1 - p[pp:(pp+(nrow(roi.feat.snv)/nind)-1)]))
      
      k = length(unique(maf.snv$sid[maf.ovl.m[maf.ovl.m[ ,1] == x.idx, 2]]))
      pval = poibin::ppoibin(length(p.bg) - k, 1 - p.bg, method = 'RF')
      
      rm(roi.feat.snv)
      gc()
      
      return(c(x, pval, x.len.snv, mean(p.bg), k))
      
    }
    
    if (length(rois.to.run) == 0) {
      
      print("No hotspots found")
      mut.rec.hotspot = mut.rec.hotspot2 = ann.results.snv = NULL
      
    } else {
      
      results = parallel::mclapply(rois.to.run, function(xr) { process.single.roi(xr) }, mc.cores = cores, mc.preschedule = FALSE, mc.silent = FALSE)
      
      # Error checking and remove entries with error
      results.error = sapply(results, function(x) { length(x) == 1 } )
      print(rois.to.run[which(results.error)])
      print(results[results.error])
      results = results[!results.error]
      
      results = do.call(rbind, results)
      rownames(results) = results[ ,1]
      # drop = FALSE --> forces matrix if single row / vector
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
          ovl.mut = IRanges::findOverlaps(maf.snv, hotspots)
          hotspots2 = hotspots[S4Vectors::subjectHits(ovl.mut)]
          hotspots2$sample = maf.snv[S4Vectors::queryHits(ovl.mut)]$sid
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
      
      # Plot hotspots manhattan and top hits
      manhattan.snv = paste(output.dir, "snv_manhattan.pdf", sep = "")
      
      if (!is.null(mut.rec.hotspot2)) {
        
        print("Plot manhattan figure for SNV")
        pdf(manhattan.snv)
        plot_manhattan(hotspots.file = mut.rec.hotspot2, fdr.cutoff = fdr.cutoff, color.line = color.line, color.dots = color.dots)
        dev.off()
        
        print("Plot top hotspots for SNV")
        plot_top_hits(hotspots.file = mut.rec.hotspot2, fdr.cutoff = fdr.cutoff, color.muts = color.muts, mutations.file  = snv.mutations.file, mutation.type = "SNV", top.no = top.no, output.dir = output.dir)
        
        # Annotate hotspots
        ann.results.snv = mutAnnotate(hotspots.file = mut.rec.hotspot2, promoter.file = promoter.file, 
                                      utr3.file = utr3.file, utr5.file = utr5.file,
                                      other.annotations = other.annotations)
        
      } else if (nrow(mut.rec.hotspot) != 0) {
        
        print("Plot manhattan figure for SNV")
        pdf(manhattan.snv)
        plot_manhattan(hotspots.file = mut.rec.hotspot, fdr.cutoff = fdr.cutoff, color.line = color.line, color.dots = color.dots)
        dev.off()
        
        print("Plot top hotspots for SNV")
        plot_top_hits(hotspots.file = mut.rec.hotspot, fdr.cutoff = fdr.cutoff, color.muts = color.muts, mutations.file  = snv.mutations.file, mutation.type = "SNV", top.no = top.no, output.dir = output.dir)
        
        # Annotate hotspots
        ann.results.snv = mutAnnotate(hotspots.file = mut.rec.hotspot, promoter.file = promoter.file, 
                                      utr3.file = utr3.file, utr5.file = utr5.file,
                                      other.annotations = other.annotations)
        
      } else {
        
        ann.results.snv = NULL
        
      }
      
    }
    
  }
  
  return(list(mut.rec.hotspot, mut.rec.hotspot2, ann.results.snv))
  
}

