#' Calculate local mutation rate.
#' 
#' @param snv.mutations.file SNV mutations MAF file.
#' @param indel.mutations.file Indel mutations MAF file.
#' @return A list containing binned SNV local mutation rates and binned indel local mutation rates.
#' @export

local.mutrate = function(snv.mutations.file, indel.mutations.file){
  # local.mutrate = function(snv.mutations.file, indel.mutations.file, mask.regions.file = system.file("extdata", "mask_regions.RDS", package = "mutrec2"), cores = 1){
    
  # If SNV mutations available, else skip this
  if (!is.null(snv.mutations.file)) {
    
    print("Calculate local mutation rate for SNV...")
    
    # Define SNV mutations
    maf.snv.mutations <- maf.to.granges(snv.mutations.file)
    
    # Chr1-ChrX
    seqi = GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[intersect(GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens))[1:23], as.character(GenomeInfoDb::seqnames(maf.snv.mutations)))]
    nind.snv = length(unique(maf.snv.mutations$sid))
    
    # Tile genome into equally sized bins (100Kb)
    binsize = 100*1e3
    genome.bins <- GenomicRanges::tileGenome(seqi, tilewidth = binsize, cut.last.tile.in.chrom = TRUE)
    nbin = length(genome.bins)
    names(genome.bins) = paste("n", seq(1:nbin), sep = "")
    genome.bins.grl <- GenomicRanges::split(genome.bins, names(genome.bins))
    
    # # Define masked region i.e. CDS, immunoglobulin loci and nonmappable
    # mask.regions = readRDS(mask.regions.file)
    
    # # Remove bins which are in masked region
    # genome.bins.grl.masked = subtract.regions.from.roi(genome.bins.grl, mask.regions, cores = cores)
    
    genome.bins.grl.masked = genome.bins.grl
    
    # Calculate mutation rate for each bin
    genome.bins.length = sum(GenomicRanges::width(genome.bins.grl.masked))
    genome.bins.length = data.frame(name = names(genome.bins.length), length = genome.bins.length)
    maf.ovl <- IRanges::findOverlaps(genome.bins.grl.masked, maf.snv.mutations)
    maf.ovl.m = IRanges::as.matrix(maf.ovl)
    genome.bins.mutcount = tapply(maf.ovl.m[ ,2], names(genome.bins.grl.masked)[maf.ovl.m[ ,1]], function(s) length(s))
    genome.bins.mutcount = data.frame(name = names(genome.bins.mutcount), mutcount = genome.bins.mutcount)
    genome.bins.mutrate = merge(genome.bins.length, genome.bins.mutcount, all = T)
    
    
    genome.bins.mutrate[is.na(genome.bins.mutrate)] = 0
    genome.bins.mutrate$mut.rate = genome.bins.mutrate$mutcount / genome.bins.mutrate$length / nind.snv
    
    mean.mutrate = sum(genome.bins.mutrate$mutcount) / sum(as.numeric(genome.bins.mutrate$length)) / nind.snv
    genome.bins.mutrate2 = rep(mean.mutrate, length(genome.bins.grl))
    names(genome.bins.mutrate2) = names(genome.bins)
    genome.bins.mutrate2[as.character(genome.bins.mutrate$name)] = genome.bins.mutrate$mut.rate
    GenomicRanges::values(genome.bins) = data.frame(mutrate = genome.bins.mutrate2[names(genome.bins)])
    
    # Adjust start coordinate to be 0-based
    df.snv <- data.frame(seqnames = GenomeInfoDb::seqnames(genome.bins), starts = as.integer(GenomicRanges::start(genome.bins) - 1), ends = GenomicRanges::end(genome.bins), mutrate = genome.bins$mutrate)

    # Divide local mutation rate into n bins and calculate the mean rate in each bin
    nbins = 10
    mutrateBins = quantile(df.snv$mutrate, seq(0, 1, 1 / nbins))
    mutrateBins[1] = mutrateBins[1] - 1e-5
    mutrate.bin = .bincode(df.snv$mutrate, mutrateBins)
    bin.mean = lapply(1:nbins, function(x) mean(df.snv$mutrate[mutrate.bin == x]))

    feat.bins = .bincode(df.snv$mutrate, mutrateBins)
    out.snv = sapply(feat.bins, function(x) bin.mean[x])
    
    df.snv$mutrate = unlist(out.snv)
    
  } else {
    
    df.snv = NULL
    
    }
  
  if (!is.null(indel.mutations.file)){
    
    ## Define indel mutations
    maf.indel.mutations <- maf.to.granges(indel.mutations.file)
    
    # Chr1-ChrX
    seqi = GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[intersect(GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens))[1:23], as.character(GenomeInfoDb::seqnames(maf.indel.mutations)))]
    nind.indel = length(unique(maf.indel.mutations$sid))
    
    # Tile genome into equally sized bins (100kb)
    binsize = 100*1e3
    genome.bins <- GenomicRanges::tileGenome(seqi, tilewidth = binsize, cut.last.tile.in.chrom = TRUE)
    nbin = length(genome.bins)
    names(genome.bins) = paste("n", seq(1:nbin), sep = "")
    genome.bins.grl <- GenomicRanges::split(genome.bins, names(genome.bins))
    
    # # Define masked region i.e. CDS, immunoglobulin loci and nonmappable
    # mask.regions = readRDS(mask.regions.file)
    
    # # Remove bins which are in masked region
    # genome.bins.grl.masked = subtract.regions.from.roi(genome.bins.grl, mask.regions, cores = cores)
    
    genome.bins.grl.masked = genome.bins.grl
    
    # Calculate mutation rate for each bin
    genome.bins.length = sum(GenomicRanges::width(genome.bins.grl.masked))
    genome.bins.length = data.frame(name = names(genome.bins.length), length = genome.bins.length)
    maf.ovl <- IRanges::findOverlaps(genome.bins.grl.masked, maf.indel.mutations)
    maf.ovl.m = IRanges::as.matrix(maf.ovl)
    genome.bins.mutcount = tapply(maf.ovl.m[ ,2], names(genome.bins.grl.masked)[maf.ovl.m[ ,1]], function(s) length(s))
    genome.bins.mutcount = data.frame(name = names(genome.bins.mutcount), mutcount = genome.bins.mutcount)
    genome.bins.mutrate = merge(genome.bins.length, genome.bins.mutcount, all=T)
    genome.bins.mutrate[is.na(genome.bins.mutrate)] = 0
    genome.bins.mutrate$mut.rate = genome.bins.mutrate$mutcount / genome.bins.mutrate$length / nind.indel
    
    mean.mutrate = sum(genome.bins.mutrate$mutcount) / sum(as.numeric(genome.bins.mutrate$length)) / nind.indel
    genome.bins.mutrate2 = rep(mean.mutrate, length(genome.bins.grl))
    names(genome.bins.mutrate2) = names(genome.bins)
    genome.bins.mutrate2[as.character(genome.bins.mutrate$name)] = genome.bins.mutrate$mut.rate
    GenomicRanges::values(genome.bins) = data.frame(mutrate = genome.bins.mutrate2[names(genome.bins)])
    
    # Adjust start coordinate to be 0-based
    df.indel <- data.frame(seqnames = GenomeInfoDb::seqnames(genome.bins), starts = as.integer(GenomicRanges::start(genome.bins) - 1), ends = GenomicRanges::end(genome.bins), mutrate = genome.bins$mutrate)

    # Divide local mutation rate in to n bins and calculate the mean rate in each bin
    nbins = 10
    mutrateBins = quantile(df.indel$mutrate, seq(0, 1, 1 / nbins))
    mutrateBins[1] = mutrateBins[1] - 1e-5
    mutrate.bin = .bincode(df.indel$mutrate, mutrateBins)
    bin.mean = lapply(1:nbins, function(x) mean(df.indel$mutrate[mutrate.bin == x]))
    
    feat.bins = .bincode(df.indel$mutrate, mutrateBins)
    out.indel = sapply(feat.bins, function(x) bin.mean[x])
    
    df.indel$mutrate = unlist(out.indel)
    
  } else {
    
    df.indel = NULL
    
    }
  
  return(list(df.snv,df.indel))
  
}

