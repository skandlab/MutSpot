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
#' @param genome.size Total number of hotspots to run analysis on, default = 2533374732.
#' @param hotspots To run hotspot analysis or region-based analysis, default = TRUE.
#' @param merge.hotspots To plot overlapping hotspots as 1 hotspot or individual hotspots, default = TRUE.
#' @return Dataframe containing predicted hotspots significance with hotspots information for indel.
#' @export

mutPredict.indel = function(mask.regions.file = system.file("extdata", "mask_regions.RDS", package = "MutSpot"), continuous.features.selected.indel.url.file, discrete.features.selected.indel.url.file, sample.specific.features.url.file = NULL, indel.mutations.file, indel.mutations.file2, indel.model.file, region.of.interest, cores = 1, min.count = 2, genome.size = 2533374732, hotspots = TRUE, merge.hotspots = TRUE) {

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
maf.masked.indel <- maf.indel[-S4Vectors::subjectHits(IRanges::findOverlaps(mask.regions, maf.indel))]
dupl.indel = duplicated(maf.masked.indel)
maf.uniq.indel = maf.masked.indel[!dupl.indel, ]

if (length(maf.uniq.indel) == 0) {
  
  mut.rec.hotspot <- data.frame(chrom=character(),
                                start=integer(),
                                end=integer(),
                                pval=numeric(),
                                length=integer(),
                                p.bg=numeric(),
                                k=integer(),
                                fdr=numeric(),
                                stringsAsFactors=FALSE)
  
  mut.rec.hotspot2=NULL
  
} else {

# Extract the middle nucleotide for each indel
GenomicRanges::start(maf.uniq.indel) = GenomicRanges::start(maf.uniq.indel) + ceiling((GenomicRanges::width(maf.uniq.indel) - 1) / 2)
GenomicRanges::end(maf.uniq.indel) = GenomicRanges::start(maf.uniq.indel)

# Extend each mutation with +/- 10 bases to define hotspots
mut.regions = maf.uniq.indel + 10
names(mut.regions) = paste("Mutation", c(1:length(mut.regions)), sep = "")

rm(list = c("maf.uniq.indel", "mask.regions"))
gc()

load(file = indel.model.file)
LRmodel.indel = LRmodel

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

# Run hotspot recurrence analysis
if (is.null(region.of.interest)) {
  
  mut.rec <- mutPredict.indel.run.lr(roi = mut.regions, maf.indel = maf.indel, maf.indel2 = maf.indel2, model.indel = LRmodel.indel, continuous.features.indel = continuous.selected.features.indel, discrete.features.indel = discrete.selected.features.indel, sample.specific.features = sample.specific.features, continuous.sample.specific = continuous.sample.specific, min.count = min.count, genome.size = genome.size, cores = cores)
  
  mut.regions2 = mut.regions[names(mut.regions) %in% rownames(mut.rec)]
  mut.rec.hotspot = data.frame(chrom = as.character(GenomeInfoDb::seqnames(mut.regions2[rownames(mut.rec)])), start = GenomicRanges::start(mut.regions2[rownames(mut.rec)]), end = GenomicRanges::end(mut.regions2[rownames(mut.rec)]), mut.rec)
  
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
  
  mut.rec <- mutPredict.indel.run.lr(roi = mut.regions, maf.indel = maf.indel, maf.indel2 = maf.indel2, model.indel = LRmodel.indel, continuous.features.indel = continuous.selected.features.indel, discrete.features.indel = discrete.selected.features.indel, sample.specific.features = sample.specific.features, continuous.sample.specific = continuous.sample.specific, min.count = min.count, genome.size = sum(GenomicRanges::width(maf.masked.regions)), cores = cores)
  
  mut.regions2 = mut.regions[names(mut.regions) %in% rownames(mut.rec)]
  mut.rec.hotspot = data.frame(chrom = as.character(GenomeInfoDb::seqnames(mut.regions2[rownames(mut.rec)])), start = GenomicRanges::start(mut.regions2[rownames(mut.rec)]), end = GenomicRanges::end(mut.regions2[rownames(mut.rec)]), mut.rec)

  } else {
    
    mut.rec <- mutPredict.indel.run.lr(roi = maf.masked.regions, maf.indel = maf.indel, maf.indel2 = maf.indel2, model.indel = LRmodel.indel, continuous.features.indel = continuous.selected.features.indel, discrete.features.indel = discrete.selected.features.indel, sample.specific.features = sample.specific.features, continuous.sample.specific = continuous.sample.specific, min.count = min.count, genome.size = length(maf.masked.regions), cores = cores)
    
    mut.regions2 = maf.masked.regions[names(maf.masked.regions) %in% rownames(mut.rec)]
    mut.rec.hotspot = data.frame(chrom = as.character(GenomeInfoDb::seqnames(mut.regions2[rownames(mut.rec)])), start = GenomicRanges::start(mut.regions2[rownames(mut.rec)]), end = GenomicRanges::end(mut.regions2[rownames(mut.rec)]), mut.rec)
    
  }
  
  }

# merge overlapping hotspots in mut.rec.hotspot, assign smallest pvalue and recalculate k
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

}

return(list(mut.rec.hotspot, mut.rec.hotspot2))

}
