#' Run mutation recurrence for SNV.
#' 
#' @param mask.regions.file Regions to mask in genome, for example, non-mappable regions/immunoglobin loci/CDS regions RDS file, depends on genome build, default file = mask_regions.RDS, Ch37.
#' @param nucleotide.selected.file Nucleotide context selected for model RDS file.
#' @param continuous.features.selected.snv.url.file Text file containing URLs of SNV continuous features selected for model.
#' @param discrete.features.selected.snv.url.file Text file containing URLs of SNV discrete features selected for model.
#' @param sample.specific.features.url.file Text file containing URLs of sample specific features, default = NULL.
#' @param snv.mutations.file SNV mutations found in region of interest MAF file.
#' @param snv.mutations.file2 SNV mutations MAF file.
#' @param collapse.regions To collapse region of interest or not, default = FALSE.
#' @param region.of.interest Region of interest bed file, default = NULL.
#' @param cores Number of cores, default = 1.
#' @param snv.model.file SNV model.
#' @param min.count Minimum number of mutated samples in each hotspot, default = 2.
#' @param hotspot.size Size of each hotspot, default = 21.
#' @param genome.size Genome size, depends on genome build, default = 2533374732, Ch37.
#' @param hotspots To run hotspot analysis or region-based analysis, default = TRUE.
#' @param output.dir Save temporary files in given output directory.
#' @param genome.build Reference genome build, default = Ch37.
#' @return Prepare intermediate files for hotspot prediction.
#' @export

mutPredict.snv.prepare = function(mask.regions.file = system.file("extdata", "mask_regions.RDS", package = "MutSpot"), nucleotide.selected.file, continuous.features.selected.snv.url.file, discrete.features.selected.snv.url.file, sample.specific.features.url.file = NULL, snv.mutations.file, snv.mutations.file2, collapse.regions = FALSE, region.of.interest, cores = 1, snv.model.file, min.count = 2, hotspot.size = 21, genome.size = 2533374732, hotspots = TRUE, output.dir, genome.build = "Ch37") {
  
  chrOrder <- c(paste("chr", 1:22, sep = ""), "chrX")
  if (genome.build == "Ch37") {
    
  seqi = GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[GenomeInfoDb::seqnames(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[1:23]]
  seqnames = GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens))[1:23]
  
  } else if (genome.build == "Ch38") {
    
    seqi = GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[GenomeInfoDb::seqnames(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[1:23]]
    seqnames = GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg38::Hsapiens))[1:23]
    
  }
  
  # Define masked region i.e. CDS, immunoglobulin loci and nonmappable
  mask.regions = readRDS(mask.regions.file)
  mask.regions = mask.regions[as.character(GenomeInfoDb::seqnames(mask.regions)) %in% seqnames]
  
  # Define SNV mutations in region of interest
  maf.snv <- maf.to.granges(snv.mutations.file)
  maf.snv = maf.snv[as.character(GenomeInfoDb::seqnames(maf.snv)) %in% seqnames]
  GenomeInfoDb::seqlevels(maf.snv) = as.character(unique(GenomeInfoDb::seqnames(maf.snv)))
  
  # Define SNV sample mutation count based on full SNV mutations file
  maf.snv2 <- maf.to.granges(snv.mutations.file2)
  maf.snv2 = maf.snv2[as.character(GenomeInfoDb::seqnames(maf.snv2)) %in% seqnames]
  maf.ind.snv = GenomicRanges::split(maf.snv2, maf.snv2$sid)
  ind.mut.count.snv = sapply(maf.ind.snv, length)
  rm(maf.snv2)
  rm(maf.ind.snv)
  
  # Define sample-specific features e.g. CIN index, COSMIC signatures
  if (!is.null(sample.specific.features.url.file)) {
    
    sample.specific.features = read.delim(sample.specific.features.url.file,stringsAsFactors = FALSE)
    if (nrow(sample.specific.features) != 0) {
      
      rownames(sample.specific.features) = as.character(sample.specific.features$SampleID)
      sample.specific.features=sample.specific.features[which(sample.specific.features$SampleID %in% names(ind.mut.count.snv)), ]
      sample.specific.features$ind.mut.count = ind.mut.count.snv[rownames(sample.specific.features)]
      sample.specific.features = sample.specific.features[ ,-which(colnames(sample.specific.features) == "SampleID")]
      
      sample.specific.features2 = parallel::mclapply(1:ncol(sample.specific.features), FUN = function(x) {
        
        print(colnames(sample.specific.features)[x])
        if (class(sample.specific.features[ ,x]) == "character") {
          
          t = factor(sample.specific.features[ ,x])
          t = model.matrix( ~ t)[ ,-1]
          if (class(t) == "matrix") {
            
            colnames(t) = substr(colnames(t), 2, nchar(colnames(t)))
            colnames(t) = paste(colnames(sample.specific.features)[x], colnames(t), sep = "")
            rownames(t) = rownames(sample.specific.features)
            
          } else {
            
            t = as.data.frame(t)
            colnames(t) = paste(colnames(sample.specific.features)[x], levels(factor(sample.specific.features[ ,x]))[2], sep = "")
            rownames(t) = rownames(sample.specific.features)
            
          }
          
        } else {
          
          t = as.data.frame(sample.specific.features[ ,x])
          colnames(t) = colnames(sample.specific.features)[x]
          rownames(t) = rownames(sample.specific.features)
          
        }
        return(t)
        
      }, mc.cores = cores)
      
      sample.specific.features = do.call(cbind, sample.specific.features2)
      
      sample.specific.urls <- read.delim(sample.specific.features.url.file, stringsAsFactors = FALSE)
      continuous.sample.specific = NULL
      for (j in 1:ncol(sample.specific.urls)) {
        
        if (class(sample.specific.urls[ ,j]) != "character"){
          
          continuous.sample.specific = c(continuous.sample.specific, colnames(sample.specific.urls)[j])
          
        }
        
      }
      
    } else {
      
      sample.specific.features = as.data.frame(ind.mut.count.snv)
      colnames(sample.specific.features) = "ind.mut.count"
      continuous.sample.specific = "ind.mut.count"
    }
    
  } else {
    
    sample.specific.features = as.data.frame(ind.mut.count.snv)
    colnames(sample.specific.features) = "ind.mut.count"
    continuous.sample.specific = "ind.mut.count"
    
  }
  
  # Remove masked regions from SNV mutations
  if (length(IRanges::findOverlaps(mask.regions, maf.snv)) != 0) {
    
    maf.masked.snv <- maf.snv[-S4Vectors::subjectHits(IRanges::findOverlaps(mask.regions, maf.snv))]
    
  } else {
    
    maf.masked.snv = maf.snv
    
  }
  dupl.snv = duplicated(maf.masked.snv)
  maf.uniq.snv = maf.masked.snv[!dupl.snv, ]
  
  if (!length(maf.uniq.snv) == 0) {
    
    # Extend each mutation with +/- 10 bases to define hotspots
    mut.regions = maf.uniq.snv + ceiling((hotspot.size - 1) / 2)
    names(mut.regions) = paste("Mutation", c(1:length(mut.regions)), sep = "")
    
    rm(list = c("maf.uniq.snv", "mask.regions"))
    gc()
    
    # load(file = snv.model.file)
    # LRmodel.snv = LRmodel
    # rm(LRmodel)
    
    # Define selected nucleotide contexts
    if (!is.null(nucleotide.selected.file)) {
      
      nucleotide.selected = readRDS(nucleotide.selected.file)
      if (!is.null(nucleotide.selected)) {
        
        nucleotide.selected = as.data.frame(nucleotide.selected)
        colnames(nucleotide.selected) = "sequence"
        nucleotide.selected$type = "type"
        nucleotide.selected$type = ifelse(grepl("one", nucleotide.selected$sequence), "oneMer", nucleotide.selected$type)
        nucleotide.selected$type = ifelse(grepl("three", nucleotide.selected$sequence), "threeMer", nucleotide.selected$type)
        nucleotide.selected$type = ifelse(grepl("three.right", nucleotide.selected$sequence), "threeRight", nucleotide.selected$type)
        nucleotide.selected$type = ifelse(grepl("three.left", nucleotide.selected$sequence), "threeLeft", nucleotide.selected$type)
        nucleotide.selected$type = ifelse(grepl("five", nucleotide.selected$sequence), "fiveMer", nucleotide.selected$type)
        nucleotide.selected$type = ifelse(grepl("five.right", nucleotide.selected$sequence), "fiveRight", nucleotide.selected$type)
        nucleotide.selected$type = ifelse(grepl("five.left", nucleotide.selected$sequence), "fiveLeft", nucleotide.selected$type)
        nucleotide.selected$sequence = gsub("one", "", nucleotide.selected$sequence)
        nucleotide.selected$sequence = gsub("three", "", nucleotide.selected$sequence)
        nucleotide.selected$sequence = gsub("five", "", nucleotide.selected$sequence)
        nucleotide.selected$sequence = gsub("*.*right", "", nucleotide.selected$sequence)
        nucleotide.selected$sequence = gsub("*.*left", "", nucleotide.selected$sequence)
        
        sel.motif = list()
        
        if ("oneMer" %in% nucleotide.selected$type){ sel.motif$oneMer = nucleotide.selected[which(nucleotide.selected$type == "oneMer"), "sequence"] }
        if ("threeMer" %in% nucleotide.selected$type){ sel.motif$threeMer = nucleotide.selected[which(nucleotide.selected$type == "threeMer"), "sequence"] }
        if ("threeRight" %in% nucleotide.selected$type){ sel.motif$threeRight = nucleotide.selected[which(nucleotide.selected$type == "threeRight"),"sequence"] }
        if ("threeLeft" %in% nucleotide.selected$type){ sel.motif$threeLeft = nucleotide.selected[which(nucleotide.selected$type == "threeLeft"), "sequence"] }
        if ("fiveMer" %in% nucleotide.selected$type){ sel.motif$fiveMer = nucleotide.selected[which(nucleotide.selected$type == "fiveMer"), "sequence"] }
        if ("fiveRight" %in% nucleotide.selected$type){ sel.motif$fiveRight = nucleotide.selected[which(nucleotide.selected$type == "fiveRight"), "sequence"] }
        if ("fiveLeft" %in% nucleotide.selected$type){ sel.motif$fiveLeft = nucleotide.selected[which(nucleotide.selected$type == "fiveLeft"), "sequence"] }
        
      } else {
        
        nucleotide.selected = sel.motif = NULL
        
      }
      
    } else {
      
      nucleotide.selected = sel.motif = NULL
      
    }
    
    # Read feature file paths
    if (!is.null(continuous.features.selected.snv.url.file)) {
      
      selected.continuous.urls.snv <- read.delim(continuous.features.selected.snv.url.file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
      
      if (nrow(selected.continuous.urls.snv) != 0) {
        
        continuous.selected.features.snv = parallel::mclapply(selected.continuous.urls.snv[ ,2], function(f) {
          
          print(f)
          df = bed.to.granges(as.character(f))
          
        }, mc.cores = cores)
        names(continuous.selected.features.snv) = as.character(selected.continuous.urls.snv[ ,1])
        
      } else {
        
        continuous.selected.features.snv = NULL
        
      }
      
    } else {
      
      continuous.selected.features.snv = NULL
      
    }
    
    if (!is.null(discrete.features.selected.snv.url.file)) {
      
      selected.discrete.urls.snv <- read.delim(discrete.features.selected.snv.url.file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
      if (nrow(selected.discrete.urls.snv) != 0) {
        
        discrete.selected.features.snv = parallel::mclapply(selected.discrete.urls.snv[ ,2], function(f) {
          
          print(f)
          df = bed.to.granges(as.character(f))
          
        }, mc.cores = cores)
        names(discrete.selected.features.snv) = as.character(selected.discrete.urls.snv[ ,1])
        
      } else {
        
        discrete.selected.features.snv = NULL
        
      }
    } else {
      
      discrete.selected.features.snv = NULL
      
    }
    
    continuous.sample.specific = c(continuous.sample.specific, names(continuous.selected.features.snv))
    
    gc()
    
    # Run hotspot recurrence analysis
    if (is.null(region.of.interest)) {
      
      if (any(duplicated(names(mut.regions)))) {
        
        print("Error: regions with duplcated names")
        return()
        
      }
      
      # Print number of roi
      print(paste("Using ", length(mut.regions), " ROIs from input", sep = ""))
      
      # Find overlaps between roi and maf to identify mutations in each roi
      print(">> Intersecting ROIs and MAF ...")
      maf.ovl <- IRanges::findOverlaps(mut.regions, maf.snv, ignore.strand = TRUE)
      maf.ovl.m = IRanges::as.matrix(maf.ovl)
      
      # Select rois mutated in min.count.sample i.e. look at rois with at least min.count number of samples mutated
      nsamples = length(unique(maf.snv$sid))
      # Encode sample names as integers to increase the efficiency of the code
      roi.count.mutated = tapply(maf.snv$sid[maf.ovl.m[ ,2]], names(mut.regions)[maf.ovl.m[ ,1]],function(s) length(unique(s)))
      rois.to.run = names(roi.count.mutated)[roi.count.mutated >= min.count]
      rm(roi.count.mutated)
      print(paste("Using ", length(rois.to.run), " ROIs mutated in >=", min.count, " samples", sep = ""))
      
      # Run analysis on all mutated roi with min.count
      # rois.to.run = roi.mutated
      mut.regions = mut.regions[which(names(mut.regions) %in% rois.to.run),]
      maf.ovl.m <- IRanges::findOverlaps(mut.regions, maf.snv, ignore.strand = TRUE)
      maf.ovl.m = IRanges::as.matrix(maf.ovl.m)
      
      saveRDS(mut.regions, file = paste(output.dir, "temp-1.RDS", sep = ""))
      saveRDS(maf.snv, file = paste(output.dir, "temp-2.RDS", sep = ""))
      saveRDS(maf.ovl.m, file = paste(output.dir, "temp-3.RDS", sep = ""))
      saveRDS(continuous.selected.features.snv, file = paste(output.dir, "temp-4.RDS", sep = ""))
      saveRDS(discrete.selected.features.snv, file = paste(output.dir, "temp-5.RDS", sep = ""))
      saveRDS(sel.motif, file = paste(output.dir, "temp-6.RDS", sep = ""))
      # saveRDS(LRmodel.snv,file = paste(output.dir, "temp-7.RDS", sep = ""))
      saveRDS(nucleotide.selected, file = paste(output.dir, "temp-8.RDS", sep = ""))
      saveRDS(rois.to.run,file = paste(output.dir, "temp-9.RDS", sep = ""))
      saveRDS(sample.specific.features, file = paste(output.dir, "temp-10.RDS", sep = ""))
      saveRDS(continuous.sample.specific, file = paste(output.dir, "temp-11.RDS", sep = ""))
      saveRDS(genome.size, file = paste(output.dir, "temp-12.RDS", sep = ""))
      
    } else {
      
      # Redefine hotspots if not whole genome analysis
      regions = bed.to.granges(region.of.interest)
      regions = regions[as.character(GenomeInfoDb::seqnames(regions)) %in% as.character(GenomeInfoDb::seqnames(seqi))]
      names(regions) = paste("Region", c(1:length(regions)), sep = "")
      
      # Define masked region i.e. CDS, immunoglobulin loci and nonmappable
      mask.regions = readRDS(mask.regions.file)
      
      # Remove masked regions from region of interest
      if (!collapse.regions) {
        
      maf.masked.regions <-regions[-S4Vectors::subjectHits(IRanges::findOverlaps(mask.regions, regions))]
      
      } else {
        
        maf.masked.regions <- subtract.regions.from.roi(regions, mask.regions, cores = cores)
        
      }
      
      if (hotspots) {
        
        if (any(duplicated(names(mut.regions)))) {
          
          print("Error: regions with duplcated names")
          return()
          
        }
        
        # Print number of roi
        print(paste("Using ", length(mut.regions), " ROIs from input", sep = ""))
        
        # Find overlaps between roi and maf to identify mutations in each roi
        print(">> Intersecting ROIs and MAF ...")
        maf.ovl <- IRanges::findOverlaps(mut.regions, maf.snv, ignore.strand = TRUE)
        maf.ovl.m = IRanges::as.matrix(maf.ovl)
        
        # Select rois mutated in min.count.sample i.e. look at rois with at least min.count number of samples mutated
        nsamples = length(unique(maf.snv$sid))
        # Encode sample names as integers to increase the efficiency of the code
        roi.count.mutated = tapply(maf.snv$sid[maf.ovl.m[ ,2]], names(mut.regions)[maf.ovl.m[ ,1]],function(s) length(unique(s)))
        rois.to.run = names(roi.count.mutated)[roi.count.mutated >= min.count]
        rm(roi.count.mutated)
        print(paste("Using ", length(rois.to.run), " ROIs mutated in >=", min.count, " samples", sep = ""))
        
        # Run analysis on all mutated roi with min.count
        # rois.to.run = roi.mutated
        mut.regions = mut.regions[which(names(mut.regions) %in% rois.to.run),]
        maf.ovl.m <- IRanges::findOverlaps(mut.regions, maf.snv, ignore.strand = TRUE)
        maf.ovl.m = IRanges::as.matrix(maf.ovl.m)
        
        genome.size = sum(GenomicRanges::width(GenomicRanges::reduce(maf.masked.regions)))
        
        saveRDS(mut.regions, file = paste(output.dir, "temp-1.RDS", sep = ""))
        saveRDS(maf.snv, file = paste(output.dir, "temp-2.RDS", sep = ""))
        saveRDS(maf.ovl.m, file = paste(output.dir, "temp-3.RDS", sep = ""))
        saveRDS(continuous.selected.features.snv, file = paste(output.dir, "temp-4.RDS", sep = ""))
        saveRDS(discrete.selected.features.snv, file = paste(output.dir, "temp-5.RDS", sep = ""))
        saveRDS(sel.motif, file = paste(output.dir, "temp-6.RDS", sep = ""))
        # saveRDS(LRmodel.snv,file = paste(output.dir, "temp-7.RDS", sep = ""))
        saveRDS(nucleotide.selected, file = paste(output.dir, "temp-8.RDS", sep = ""))
        saveRDS(rois.to.run,file = paste(output.dir, "temp-9.RDS", sep = ""))
        saveRDS(sample.specific.features, file = paste(output.dir, "temp-10.RDS", sep = ""))
        saveRDS(continuous.sample.specific, file = paste(output.dir, "temp-11.RDS", sep = ""))
        saveRDS(genome.size, file = paste(output.dir, "temp-12.RDS", sep = ""))
        
      } else {
        
        # mut.rec <- mutPredict.snv.run.lr(roi = maf.masked.regions, maf.snv = maf.snv, maf.snv2 = maf.snv2, model.snv = LRmodel.snv, continuous.features.snv = continuous.selected.features.snv, discrete.features.snv = discrete.selected.features.snv, motifs = sel.motif, nucleotide.selected = nucleotide.selected, sample.specific.features = sample.specific.features, continuous.sample.specific = continuous.sample.specific, min.count = min.count, genome.size = length(maf.masked.regions), cores = cores)
        # 
        # mut.regions2 = maf.masked.regions[names(maf.masked.regions) %in% rownames(mut.rec)]
        # mut.rec.hotspot = data.frame(chrom = as.character(GenomeInfoDb::seqnames(mut.regions2[rownames(mut.rec)])), start = GenomicRanges::start(mut.regions2[rownames(mut.rec)]), end = GenomicRanges::end(mut.regions2[rownames(mut.rec)]), mut.rec)
        
        if (any(duplicated(names(mut.regions)))) {
          
          print("Error: regions with duplcated names")
          return()
          
        }
        
        # Print number of roi
        print(paste("Using ", length(maf.masked.regions), " ROIs from input", sep = ""))
        
        # Find overlaps between roi and maf to identify mutations in each roi
        print(">> Intersecting ROIs and MAF ...")
        maf.ovl <- IRanges::findOverlaps(maf.masked.regions, maf.snv, ignore.strand = TRUE)
        maf.ovl.m = IRanges::as.matrix(maf.ovl)
        
        # Select rois mutated in min.count.sample i.e. look at rois with at least min.count number of samples mutated
        nsamples = length(unique(maf.snv$sid))
        # Encode sample names as integers to increase the efficiency of the code
        roi.count.mutated = tapply(maf.snv$sid[maf.ovl.m[ ,2]], names(maf.masked.regions)[maf.ovl.m[ ,1]],function(s) length(unique(s)))
        rois.to.run = names(roi.count.mutated)[roi.count.mutated >= min.count]
        rm(roi.count.mutated)
        print(paste("Using ", length(rois.to.run), " ROIs mutated in >=", min.count, " samples", sep = ""))
        
        # Run analysis on all mutated roi with min.count
        # rois.to.run = roi.mutated
        mut.regions = maf.masked.regions[which(names(maf.masked.regions) %in% rois.to.run),]
        maf.ovl.m <- IRanges::findOverlaps(mut.regions, maf.snv, ignore.strand = TRUE)
        maf.ovl.m = IRanges::as.matrix(maf.ovl.m)
        
        genome.size = length(maf.masked.regions)
        
        saveRDS(mut.regions, file = paste(output.dir, "temp-1.RDS", sep = ""))
        saveRDS(maf.snv, file = paste(output.dir, "temp-2.RDS", sep = ""))
        saveRDS(maf.ovl.m, file = paste(output.dir, "temp-3.RDS", sep = ""))
        saveRDS(continuous.selected.features.snv, file = paste(output.dir, "temp-4.RDS", sep = ""))
        saveRDS(discrete.selected.features.snv, file = paste(output.dir, "temp-5.RDS", sep = ""))
        saveRDS(sel.motif, file = paste(output.dir, "temp-6.RDS", sep = ""))
        # saveRDS(LRmodel.snv,file = paste(output.dir, "temp-7.RDS", sep = ""))
        saveRDS(nucleotide.selected, file = paste(output.dir, "temp-8.RDS", sep = ""))
        saveRDS(rois.to.run,file = paste(output.dir, "temp-9.RDS", sep = ""))
        saveRDS(sample.specific.features, file = paste(output.dir, "temp-10.RDS", sep = ""))
        saveRDS(continuous.sample.specific, file = paste(output.dir, "temp-11.RDS", sep = ""))
        saveRDS(genome.size, file = paste(output.dir, "temp-12.RDS", sep = ""))
        
      }
      
    }
    
  }
  
  rm(list = ls(all.names = TRUE))
  gc()
  
}