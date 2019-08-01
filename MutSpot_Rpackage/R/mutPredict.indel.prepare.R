#' Mutation hotspot recurrence prediction for indel.
#' 
#' @param mask.regions.file Regions to mask in genome, for example, non-mappable regions/immunoglobin loci/CDS regions RDS file, default file = mask_regions.RDS.
#' @param continuous.features.selected.indel.url.file Text file containing URLs of indel continuous features selected for model.
#' @param discrete.features.selected.indel.url.file Text file containing URLs of indel discrete features selected for model.
#' @param sample.specific.features.url.file Text file containing URLs of sample specific features, default = NULL.
#' @param indel.mutations.file Indel mutations found in region of interest MAF file.
#' @param indel.mutations.file2 Indel mutations MAF file.
#' @param indel.model.file Indel model.
#' @param region.of.interest Region of interest bed file, default = NULL.
#' @param cores Number of cores, default = 1.
#' @param min.count Minimum number of mutated samples in each hotspot, default = 2.
#' @param hotspot.size Size of each hotspot, default = 21.
#' @param genome.size Total number of hotspots to run analysis on, default = 2533374732.
#' @param hotspots To run hotspot analysis or region-based analysis, default = TRUE.
#' @param output.dir Save plot in given output directory.
#' @return Prepare intermediate files for hotspot prediction.
#' @export

mutPredict.indel.prepare = function(mask.regions.file = system.file("extdata", "mask_regions.RDS", package = "MutSpot"), continuous.features.selected.indel.url.file, discrete.features.selected.indel.url.file, sample.specific.features.url.file = NULL, indel.mutations.file, indel.mutations.file2, indel.model.file, region.of.interest, cores = 1, min.count = 2, hotspot.size = 21, genome.size = 2533374732, hotspots = TRUE, output.dir) {

# Chr1-X
chrOrder <- c(paste("chr", 1:22, sep = ""), "chrX")
seqi = GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[GenomeInfoDb::seqnames(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[1:23]]
seqnames = GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens))[1:23]

# Define masked region i.e. CDS, immunoglobulin loci and nonmappable
mask.regions = readRDS(mask.regions.file)
mask.regions = mask.regions[as.character(GenomeInfoDb::seqnames(mask.regions)) %in% seqnames]

# Define indel mutations in region of interest
maf.indel <- maf.to.granges(indel.mutations.file)
maf.indel = maf.indel[as.character(GenomeInfoDb::seqnames(maf.indel)) %in% seqnames]
GenomeInfoDb::seqlevels(maf.indel) = as.character(unique(GenomeInfoDb::seqnames(maf.indel)))

# Define indel sample mutation count based on full indel mutations file
maf.indel2 <- maf.to.granges(indel.mutations.file2)
maf.indel2 = maf.indel2[as.character(GenomeInfoDb::seqnames(maf.indel2)) %in% seqnames]
maf.ind.indel = GenomicRanges::split(maf.indel2, maf.indel2$sid)
ind.mut.count.indel = sapply(maf.ind.indel, length)
rm(maf.indel2)
rm(maf.ind.indel)

# Define sample-specific features e.g. CIN index, COSMIC signatures
if (!is.null(sample.specific.features.url.file)) {
  
  sample.specific.features = read.delim(sample.specific.features.url.file, stringsAsFactors = FALSE)
  if (nrow(sample.specific.features) != 0) {
    
    rownames(sample.specific.features) = as.character(sample.specific.features$SampleID)
    sample.specific.features = sample.specific.features[which(sample.specific.features$SampleID %in% names(ind.mut.count.indel)), ]
    sample.specific.features$sample.count = ind.mut.count.indel[rownames(sample.specific.features)]
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
    
    sample.specific.features = as.data.frame(ind.mut.count.indel)
    colnames(sample.specific.features) = "sample.count"
    continuous.sample.specific = "sample.count"
    
  }
  
} else {
  
  sample.specific.features = as.data.frame(ind.mut.count.indel)
  colnames(sample.specific.features) = "sample.count"
  continuous.sample.specific = "sample.count"
  
}

# Remove masked regions from indel mutations
if (length(IRanges::findOverlaps(mask.regions, maf.indel)) != 0) {
  
  maf.masked.indel <- maf.indel[-S4Vectors::subjectHits(IRanges::findOverlaps(mask.regions, maf.indel))]
  
} else {
  
  maf.masked.indel = maf.indel
  
}
dupl.indel = duplicated(maf.masked.indel)
maf.uniq.indel = maf.masked.indel[!dupl.indel, ]

if (length(maf.uniq.indel) != 0) {
  
  # Extract the middle nucleotide for each indel
  GenomicRanges::start(maf.uniq.indel) = GenomicRanges::start(maf.uniq.indel) + ceiling((GenomicRanges::width(maf.uniq.indel) - 1) / 2)
  GenomicRanges::end(maf.uniq.indel) = GenomicRanges::start(maf.uniq.indel)
  
  # Extend each mutation with +/- 10 bases to define hotspots
  mut.regions = maf.uniq.indel + ceiling((hotspot.size - 1) / 2)
  names(mut.regions) = paste("Mutation", c(1:length(mut.regions)), sep = "")
  
  rm(list = c("maf.uniq.indel", "mask.regions"))
  gc()
  
  # load(file = indel.model.file)
  # LRmodel.indel = LRmodel
  # rm(LRmodel)
  
  # Read feature file paths
  if (!is.null(continuous.features.selected.indel.url.file)) {
    
    selected.continuous.urls.indel <- read.delim(continuous.features.selected.indel.url.file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    if (nrow(selected.continuous.urls.indel) != 0) {
      
      continuous.selected.features.indel = parallel::mclapply(selected.continuous.urls.indel[ ,2], function(f) {
        
        print(f)
        df = bed.to.granges(as.character(f))
        
      }, mc.cores = cores)
      names(continuous.selected.features.indel) = as.character(selected.continuous.urls.indel[ ,1])
      
    } else {
      
      continuous.selected.features.indel = NULL
      
    }
    
  } else {
    
    continuous.selected.features.indel = NULL
    
  }
  
  if (!is.null(discrete.features.selected.indel.url.file)) {
    
    selected.discrete.urls.indel <- read.delim(discrete.features.selected.indel.url.file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    if (nrow(selected.discrete.urls.indel) != 0) {
      
      discrete.selected.features.indel = parallel::mclapply(selected.discrete.urls.indel[ ,2], function(f) {
        
        print(f)
        df = bed.to.granges(as.character(f))
        
      }, mc.cores = cores)
      names(discrete.selected.features.indel) = as.character(selected.discrete.urls.indel[ ,1])
      
    } else {
      
      discrete.selected.features.indel = NULL
      
    }
    
  } else {
    
    discrete.selected.features.indel = NULL
    
  }
  
  continuous.sample.specific = c(continuous.sample.specific, names(continuous.selected.features.indel))

  gc()
  
  # Run hotspot recurrence analysis
  if (is.null(region.of.interest)) {
    
    # Check that roi names are unique
    if (any(duplicated(names(mut.regions)))) {
      
      print("Error: regions with duplcated names")
      return()
      
    }
    
    # Print number of roi
    print(paste("Using ", length(mut.regions), " ROIs from input", sep = ""))
    
    # Find overlaps between roi and maf to identify mutations in each roi
    print(">> Intersecting ROIs and MAF ...")
    maf.ovl <- IRanges::findOverlaps(mut.regions, maf.indel, ignore.strand = TRUE)
    maf.ovl.m = IRanges::as.matrix(maf.ovl)
    
    # Select rois mutated in min.count.sample i.e. look at rois with at least min.count number of samples mutated
    nsamples = length(unique(maf.indel$sid))
    # Encode sample names as integers to increase the efficiency of the code
    roi.count.mutated = tapply(maf.indel$sid[maf.ovl.m[ ,2]], names(mut.regions)[maf.ovl.m[ ,1]], function(s) length(unique(s)))
    rois.to.run = names(roi.count.mutated)[roi.count.mutated >= min.count]
    rm(roi.count.mutated)
    print(paste("Using ", length(rois.to.run), " ROIs mutated in >=", min.count, " samples", sep = ""))
    
    # Run analysis on all mutated roi with min.count
    # rois.to.run = roi.mutated
    mut.regions = mut.regions[which(names(mut.regions) %in% rois.to.run),]
    maf.ovl.m <- IRanges::findOverlaps(mut.regions, maf.indel, ignore.strand = TRUE)
    maf.ovl.m = IRanges::as.matrix(maf.ovl.m)
    
    saveRDS(mut.regions, file = paste(output.dir, "temp_1.RDS", sep = ""))
    saveRDS(maf.indel, file = paste(output.dir, "temp_2.RDS", sep = ""))
    saveRDS(maf.ovl.m, file = paste(output.dir, "temp_3.RDS", sep = ""))
    saveRDS(continuous.selected.features.indel, file = paste(output.dir, "temp_4.RDS", sep = ""))
    saveRDS(discrete.selected.features.indel, file = paste(output.dir, "temp_5.RDS", sep = ""))
    # saveRDS(LRmodel.indel,file = paste(output.dir, "temp_6.RDS", sep = ""))
    saveRDS(rois.to.run,file = paste(output.dir, "temp_7.RDS", sep = ""))
    saveRDS(sample.specific.features, file = paste(output.dir, "temp_8.RDS", sep = ""))
    saveRDS(continuous.sample.specific, file = paste(output.dir, "temp_9.RDS", sep = ""))
    saveRDS(genome.size, file = paste(output.dir, "temp_10.RDS", sep = ""))

  } else {
    
    # Redefine hotspots if not whole genome analysis
    regions = bed.to.granges(region.of.interest)
    regions = regions[as.character(GenomeInfoDb::seqnames(regions)) %in% as.character(GenomeInfoDb::seqnames(seqi))]
    names(regions) = paste("Region", c(1:length(regions)), sep = "")
    
    # Define masked region i.e. CDS, immunoglobulin loci and nonmappable
    mask.regions = readRDS(mask.regions.file)
    
    # Remove masked regions from region of interest
    maf.masked.regions <- regions[-S4Vectors::subjectHits(IRanges::findOverlaps(mask.regions, regions))]
    
    if (hotspots) {
      
      # Check that roi names are unique
      if (any(duplicated(names(mut.regions)))) {
        
        print("Error: regions with duplcated names")
        return()
        
      }
      
      # Print number of roi
      print(paste("Using ", length(mut.regions), " ROIs from input", sep = ""))
      
      # Find overlaps between roi and maf to identify mutations in each roi
      print(">> Intersecting ROIs and MAF ...")
      maf.ovl <- IRanges::findOverlaps(mut.regions, maf.indel, ignore.strand = TRUE)
      maf.ovl.m = IRanges::as.matrix(maf.ovl)
      
      # Select rois mutated in min.count.sample i.e. look at rois with at least min.count number of samples mutated
      nsamples = length(unique(maf.indel$sid))
      # Encode sample names as integers to increase the efficiency of the code
      roi.count.mutated = tapply(maf.indel$sid[maf.ovl.m[ ,2]], names(mut.regions)[maf.ovl.m[ ,1]], function(s) length(unique(s)))
      rois.to.run = names(roi.count.mutated)[roi.count.mutated >= min.count]
      rm(roi.count.mutated)
      print(paste("Using ", length(rois.to.run), " ROIs mutated in >=", min.count, " samples", sep = ""))
      
      # Run analysis on all mutated roi with min.count
      # rois.to.run = roi.mutated
      mut.regions = mut.regions[which(names(mut.regions) %in% rois.to.run),]
      maf.ovl.m <- IRanges::findOverlaps(mut.regions, maf.indel, ignore.strand = TRUE)
      maf.ovl.m = IRanges::as.matrix(maf.ovl.m)
      
      genome.size = sum(GenomicRanges::width(maf.masked.regions))
      
      saveRDS(mut.regions, file = paste(output.dir, "temp_1.RDS", sep = ""))
      saveRDS(maf.indel, file = paste(output.dir, "temp_2.RDS", sep = ""))
      saveRDS(maf.ovl.m, file = paste(output.dir, "temp_3.RDS", sep = ""))
      saveRDS(continuous.selected.features.indel, file = paste(output.dir, "temp_4.RDS", sep = ""))
      saveRDS(discrete.selected.features.indel, file = paste(output.dir, "temp_5.RDS", sep = ""))
      # saveRDS(LRmodel.indel,file = paste(output.dir, "temp_6.RDS", sep = ""))
      saveRDS(rois.to.run,file = paste(output.dir, "temp_7.RDS", sep = ""))
      saveRDS(sample.specific.features, file = paste(output.dir, "temp_8.RDS", sep = ""))
      saveRDS(continuous.sample.specific, file = paste(output.dir, "temp_9.RDS", sep = ""))
      saveRDS(genome.size, file = paste(output.dir, "temp_10.RDS", sep = ""))
      
    } else {
      
      # Check that roi names are unique
      if (any(duplicated(names(maf.masked.regions)))) {
        
        print("Error: regions with duplcated names")
        return()
        
      }
      
      # Print number of roi
      print(paste("Using ", length(maf.masked.regions), " ROIs from input", sep = ""))
      
      # Find overlaps between roi and maf to identify mutations in each roi
      print(">> Intersecting ROIs and MAF ...")
      maf.ovl <- IRanges::findOverlaps(maf.masked.regions, maf.indel, ignore.strand = TRUE)
      maf.ovl.m = IRanges::as.matrix(maf.ovl)
      
      # Select rois mutated in min.count.sample i.e. look at rois with at least min.count number of samples mutated
      nsamples = length(unique(maf.indel$sid))
      # Encode sample names as integers to increase the efficiency of the code
      roi.count.mutated = tapply(maf.indel$sid[maf.ovl.m[ ,2]], names(maf.masked.regions)[maf.ovl.m[ ,1]], function(s) length(unique(s)))
      rois.to.run = names(roi.count.mutated)[roi.count.mutated >= min.count]
      rm(roi.count.mutated)
      print(paste("Using ", length(rois.to.run), " ROIs mutated in >=", min.count, " samples", sep = ""))
      
      # Run analysis on all mutated roi with min.count
      # rois.to.run = roi.mutated
      mut.regions = maf.masked.regions[which(names(maf.masked.regions) %in% rois.to.run),]
      maf.ovl.m <- IRanges::findOverlaps(mut.regions, maf.indel, ignore.strand = TRUE)
      maf.ovl.m = IRanges::as.matrix(maf.ovl.m)
      
      genome.size = sum(GenomicRanges::width(maf.masked.regions))
      
      saveRDS(mut.regions, file = paste(output.dir, "temp_1.RDS", sep = ""))
      saveRDS(maf.indel, file = paste(output.dir, "temp_2.RDS", sep = ""))
      saveRDS(maf.ovl.m, file = paste(output.dir, "temp_3.RDS", sep = ""))
      saveRDS(continuous.selected.features.indel, file = paste(output.dir, "temp_4.RDS", sep = ""))
      saveRDS(discrete.selected.features.indel, file = paste(output.dir, "temp_5.RDS", sep = ""))
      # saveRDS(LRmodel.indel,file = paste(output.dir, "temp_6.RDS", sep = ""))
      saveRDS(rois.to.run,file = paste(output.dir, "temp_7.RDS", sep = ""))
      saveRDS(sample.specific.features, file = paste(output.dir, "temp_8.RDS", sep = ""))
      saveRDS(continuous.sample.specific, file = paste(output.dir, "temp_9.RDS", sep = ""))
      saveRDS(genome.size, file = paste(output.dir, "temp_10.RDS", sep = ""))
      
    }
    
  }
  
}

rm(list = ls(all.names = TRUE))
gc()

}

