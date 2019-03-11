#' Prepare covariates matrix by chromosome for indel model fitting.
#' 
#' @param mask.regions.file Regions to mask in genome, for example, non-mappable regions/immunoglobin loci/CDS regions RDS file, default file = mask_regions.RDS.
#' @param all.sites.file All sites in whole genome RDS file, default file = all_sites.RDS.
#' @param continuous.features.selected.indel.url.file Text file containing URLs of indel continuous features selected for model.
#' @param discrete.features.selected.indel.url.file Text file containing URLs of indel discrete features selected for model.
#' @param sample.specific.features.url.file Text file containing sample specific indel features, default = NULL.
#' @param indel.mutations.file Indel mutations found in region of interest MAF file.
#' @param region.of.interest Region of interest bed file, default = NULL.
#' @param indel.mutations.file2 Indel mutations MAF file.
#' @param cores Number of cores, default = 1.
#' @param chr.interest Chromosome of interest, default = "chr1".
#' @return Saves covariates matrices for mutated sites and whole genome/specified region by chromosome.
#' @export

mutCovariate.indel = function(mask.regions.file = system.file("extdata", "mask_regions.RDS", package = "MutSpot"), all.sites.file = system.file("extdata", "all_sites.RDS", package = "MutSpot"), continuous.features.selected.indel.url.file, discrete.features.selected.indel.url.file, sample.specific.features.url.file = NULL, indel.mutations.file, region.of.interest = NULL, indel.mutations.file2, cores = 1, chr.interest = "chr1") {
  
  # Chr1-X
  chrOrder <- c(paste("chr", 1:22, sep = ""), "chrX")
seqi = GenomicRanges::seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[GenomicRanges::seqnames(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[1:23]]
seqnames = GenomicRanges::seqnames(GenomicRanges::seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens))[1:23]

# Define masked region i.e. CDS, immunoglobulin loci and nonmappable
mask.regions = readRDS(mask.regions.file)
mask.regions = mask.regions[as.character(GenomicRanges::seqnames(mask.regions)) %in% seqnames]

# Define all sites in whole genome
all.sites = readRDS(all.sites.file)
all.sites = all.sites[as.character(GenomicRanges::seqnames(all.sites)) %in% seqnames]
all.sites.masked = subtract.regions.from.roi(all.sites, mask.regions, cores = cores)
sum(as.numeric(GenomicRanges::width(all.sites.masked)))

# If specified region, redefine all sites to be in specified region
if (!is.null(region.of.interest)) {
  
  print("specified region")
  
  all.sites = read.delim(region.of.interest, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  all.sites = with(all.sites, GenomicRanges::GRanges(V1, IRanges::IRanges(V2, V3)))
  all.sites.masked = subtract.regions.from.roi(all.sites, mask.regions, cores = cores)
  all.sites.masked = all.sites.masked[as.character(GenomicRanges::seqnames(all.sites.masked)) %in% seqnames]
  
}

# Extract all polyA, poly C, poly G and poly T positions in a specific chromosome
chrs <- names(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[1:23]
polyA <- lapply(chr.interest, function(x) BSgenome::start(Biostrings::matchPattern("AAAAA",
                                                                           BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[x]])))
names(polyA) <- chr.interest
polyAs = parallel::mclapply(chr.interest, FUN = function(x) {
  
  print(x)
  gr = GenomicRanges::GRanges(x, IRanges::IRanges(polyA[[x]], polyA[[x]] + 5 - 1))
  
} )
polyAs = BiocGenerics::unlist(GenomicRanges::GRangesList(polyAs))

polyC <- lapply(chr.interest, function(x) BSgenome::start(Biostrings::matchPattern("CCCCC",
                                                                           BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[x]])))
names(polyC) <- chr.interest
polyCs = parallel::mclapply(chr.interest, FUN = function(x) {
  
  print(x)
  gr = GenomicRanges::GRanges(x, IRanges::IRanges(polyC[[x]], polyC[[x]] + 5 - 1))
  
} )
polyCs = BiocGenerics::unlist(GenomicRanges::GRangesList(polyCs))

polyG <- lapply(chr.interest, function(x) BSgenome::start(Biostrings::matchPattern("GGGGG",
                                                                           BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[x]])))
names(polyG) <- chr.interest
polyGs = parallel::mclapply(chr.interest, FUN = function(x) {
  
  print(x)
  gr = GenomicRanges::GRanges(x, IRanges::IRanges(polyG[[x]], polyG[[x]] + 5 - 1))
  
} )
polyGs = BiocGenerics::unlist(GenomicRanges::GRangesList(polyGs))

polyT <- lapply(chr.interest, function(x) BSgenome::start(Biostrings::matchPattern("TTTTT",
                                                                           BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[x]])))
names(polyT) <- chr.interest
polyTs = parallel::mclapply(chr.interest, FUN = function(x) {
  
  print(x)
  gr = GenomicRanges::GRanges(x, IRanges::IRanges(polyT[[x]], polyT[[x]] + 5 - 1))
  
} )
polyTs = BiocGenerics::unlist(GenomicRanges::GRangesList(polyTs))

# Read feature file paths
if (!is.null(continuous.features.selected.indel.url.file)) {
  
selected.continuous.urls <- read.delim(continuous.features.selected.indel.url.file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
continuous.selected.features = parallel::mclapply(selected.continuous.urls[ ,2], function(f) {
  
  print(f)
  df = read.delim(as.character(f), stringsAsFactors = FALSE, header = FALSE)
  with(df, GenomicRanges::GRanges(V1, IRanges::IRanges(V2, V3), score = V4)) }, mc.cores = cores)
names(continuous.selected.features) = as.character(selected.continuous.urls[ ,1])

} else {
  
  continuous.selected.features = NULL
  
}

if (!is.null(discrete.features.selected.indel.url.file)) {
  
selected.discrete.urls <- read.delim(discrete.features.selected.indel.url.file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
discrete.selected.features = parallel::mclapply(selected.discrete.urls[ ,2], function(f) {
  
  print(f)
  df = read.delim(as.character(f), stringsAsFactors = FALSE, header = FALSE)
  with(df, GenomicRanges::GRanges(V1, IRanges::IRanges(V2, V3))) }, mc.cores = cores)
names(discrete.selected.features) = as.character(selected.discrete.urls[ ,1])

} else {
  
  discrete.selected.features = NULL
  
}

# Define indel mutations
maf.mutations <- maf.to.granges(indel.mutations.file)
maf.mutations = maf.mutations[as.character(GenomicRanges::seqnames(maf.mutations)) %in% seqnames]
mut.masked <- maf.mutations[S4Vectors::subjectHits(IRanges::findOverlaps(all.sites.masked, maf.mutations))]
mut.masked.sites = mut.masked
GenomicRanges::start(mut.masked.sites) = GenomicRanges::start(mut.masked.sites) + ceiling((GenomicRanges::width(mut.masked.sites) - 1) / 2)
GenomicRanges::end(mut.masked.sites) = GenomicRanges::start(mut.masked.sites)

# Define indel sample mutation count based on full indel mutations file
maf.mutations2 = maf.to.granges(indel.mutations.file2)
maf.mutations2 = maf.mutations2[as.character(GenomicRanges::seqnames(maf.mutations2)) %in% seqnames]
maf.ind = GenomicRanges::split(maf.mutations2, maf.mutations2$sid)
ind.mut.count = sapply(maf.ind, length)
nind = length(ind.mut.count) 

# Define sample-specific features e.g. CIN index, COSMIC signatures
if (!is.null(sample.specific.features.url.file)) {

  sample.specific.features = read.delim(sample.specific.features.url.file,stringsAsFactors = FALSE)
  rownames(sample.specific.features)=as.character(sample.specific.features$SampleID)
  sample.specific.features=sample.specific.features[which(sample.specific.features$SampleID %in% names(ind.mut.count)),]
  sample.specific.features$sample.count=ind.mut.count[rownames(sample.specific.features)]
  sample.specific.features=sample.specific.features[,-which(colnames(sample.specific.features)=="SampleID")]
  
  sample.specific.features2=parallel::mclapply(1:ncol(sample.specific.features), FUN=function(x) {
    
    print(colnames(sample.specific.features)[x])
    if(class(sample.specific.features[,x])=="character") {
      
      t=factor(sample.specific.features[,x])
      t=model.matrix(~t)[,-1]
      colnames(t)=substr(colnames(t),2,nchar(colnames(t)))
      colnames(t)=paste(colnames(sample.specific.features)[x],colnames(t),sep="")
      rownames(t)=rownames(sample.specific.features)
      
    } else {
      
      t=as.data.frame(sample.specific.features[,x])
      colnames(t)=colnames(sample.specific.features)[x]
      rownames(t)=rownames(sample.specific.features)
      
    }
    return(t)
    
  },mc.cores=cores)
  
  sample.specific.features=do.call(cbind,sample.specific.features2)
  
} else {
  
  sample.specific.features=as.data.frame(ind.mut.count)
  colnames(sample.specific.features)="sample.count"
  
}

# Remove larger objects before tabulating
sort(sapply(ls(), function(x) { object.size(get(x)) / 10 ^ 6 } ))
rm(list = c("maf.mutations", "maf.ind", "mask.regions", "all.sites", "maf.mutations2"))
gc(reset = T)

# Tabulate covariates for mutations
GenomeInfoDb::seqlevels(mut.masked.sites) = as.character(unique(GenomicRanges::seqnames(mut.masked.sites)))
mut.chr = GenomicRanges::split(mut.masked.sites, GenomicRanges::seqnames(mut.masked.sites))
mut.chr = mut.chr[names(mut.chr) %in% chrs]

if (chr.interest %in% names(mut.chr)) {
  
mut.freq <- mutCovariate.indel.freq.table.muts(continuous.features = continuous.selected.features, discrete.features = discrete.selected.features, sample.specific.features = sample.specific.features, polyAs = polyAs, polyTs = polyTs, polyCs = polyCs, polyGs = polyGs, sites = mut.chr[[chr.interest]])
rm(list = c("mut.chr"))

} else {
  
  mut.freq = NULL
  
}

# Tabulate covariates for all positions in indels
GenomeInfoDb::seqlevels(mut.masked) = as.character(unique(GenomicRanges::seqnames(mut.masked)))
mut.indel = mut.masked
mut.indel.chr = GenomicRanges::split(mut.indel, GenomicRanges::seqnames(mut.indel))
mut.indel.chr = mut.indel.chr[names(mut.indel.chr) %in% chrs]

if (chr.interest %in% names(mut.indel.chr)) {
  
indel.freq <- mutCovariate.indel.freq.table.muts(continuous.features = continuous.selected.features, discrete.features = discrete.selected.features, sample.specific.features = sample.specific.features, polyAs = polyAs, polyTs = polyTs, polyCs = polyCs, polyGs = polyGs, sites = mut.indel.chr[[chr.interest]])
rm(list = c("mut.indel.chr"))

} else {
  
  indel.freq = NULL
  
}

if (chr.interest %in% as.character(GenomicRanges::seqnames(all.sites.masked))) {
  
all.sites.masked = all.sites.masked[as.character(GenomicRanges::seqnames(all.sites.masked)) == chr.interest]
len = sapply(GenomicRanges::split(all.sites.masked, GenomicRanges::seqnames(all.sites.masked)), length)
len = len[len != 0]
len2 = sapply(1:length(len), function(i) { sum(len[1:i]) })
len2 = c(0, len2)
chunk.size = 5000
chunks <- lapply(1:length(len), function(j) { lapply(1:ceiling(len[j] / chunk.size), function(i) ((len2[j] + (i-1) * chunk.size + 1):(len2[j] + min((i) * chunk.size, len[j])))) })
chunks = do.call(c, chunks)

genome.freq <- parallel::mclapply(chunks, function(x) mutCovariate.indel.freq.table.genome(continuous.features = continuous.selected.features, discrete.features = discrete.selected.features, polyAs = polyAs, polyTs = polyTs, polyCs = polyCs, polyGs = polyGs, sites = all.sites.masked[x]), mc.cores = cores, mc.preschedule = FALSE, mc.silent = FALSE)
genome.freq = data.table::rbindlist(genome.freq)
genome.freq = data.frame(genome.freq, check.names = F)
# Sum up number of sites with same covariates combination
if (ncol(genome.freq) > 2) {
  
  genome.freq.aggregated = aggregate(genome.freq$freq, by = genome.freq[ ,1:(ncol(genome.freq) - 1)], FUN = sum)
  
} else { # Special case: only 1 feature
  
  genome.freq.aggregated = aggregate(genome.freq$freq ~ genome.freq[,colnames(genome.freq)[1]], FUN = sum)
  colnames(genome.freq.aggregated) = c(colnames(genome.freq)[1], "x")
  
}
rm(list = c("genome.freq"))

} else {
  
  genome.freq.aggregated = NULL
  
}

saveRDS(list(mut.freq, indel.freq, genome.freq.aggregated), file = paste("mutCovariate_indel_", chr.interest, ".RDS", sep = ""))

}
