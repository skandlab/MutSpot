#' Prepare covariates matrix by chromosome for SNV model fitting.
#' 
#' @param mask.regions.file Regions to mask in genome, for example, non-mappable regions/immunoglobin loci/CDS regions RDS file, default file = mask_regions.RDS.
#' @param all.sites.file All sites in whole genome RDS file, default file = all_sites.RDS.
#' @param nucleotide.selected.file Nucleotide context selected for model RDS file.
#' @param continuous.features.selected.snv.url.file Text file containing URLs of SNV continuous features selected for model.
#' @param discrete.features.selected.snv.url.file Text file containing URLs of SNV discrete features selected for model.
#' @param sample.specific.features.url.file Text file sample specific SNV features, default = NULL.
#' @param snv.mutations.file SNV mutations found in region of interest MAF file.
#' @param region.of.interest Region of interest bed file, default = NULL.
#' @param snv.mutations.file2 SNV mutations MAF file.
#' @param cores Number of cores, default = 1.
#' @param chr.interest Chromosome of interest, default = "chr1".
#' @return Saves covariates matrices for mutated sites and whole genome/specified region by chromosome.
#' @export

mutCovariate.snv = function(mask.regions.file = system.file("extdata", "mask_regions.RDS", package = "MutSpot"), all.sites.file = system.file("extdata", "all_sites.RDS", package = "MutSpot"), nucleotide.selected.file, continuous.features.selected.snv.url.file, discrete.features.selected.snv.url.file, sample.specific.features.url.file = NULL, snv.mutations.file, region.of.interest = NULL, snv.mutations.file2, cores = 1, chr.interest = "chr1") {
  
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
# sum(as.numeric(GenomicRanges::width(all.sites.masked)))

# If specified region, redefine all sites to be in specified region
if (!is.null(region.of.interest)) {
  
  print("specified region")
  
  all.sites = read.delim(region.of.interest, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  all.sites = with(all.sites, GenomicRanges::GRanges(V1, IRanges::IRanges(V2, V3)))
  all.sites.masked = subtract.regions.from.roi(all.sites,mask.regions, cores = cores)
  all.sites.masked = all.sites.masked[as.character(GenomicRanges::seqnames(all.sites.masked)) %in% seqnames]
  
}

# If nucleotide contexts selected for model fitting, else skip this
if (!is.null(nucleotide.selected.file)) {
  
nucleotide.selected = readRDS(nucleotide.selected.file)
nucleotide.selected = as.data.frame(nucleotide.selected)
colnames(nucleotide.selected) = "sequence"
nucleotide.selected$type = "type"
nucleotide.selected$type = ifelse(grepl("one",nucleotide.selected$sequence), "oneMer", nucleotide.selected$type)
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

precompute.motif.pos = as.list(numeric(nrow(nucleotide.selected)))
names(precompute.motif.pos) = paste(nucleotide.selected$type, nucleotide.selected$sequence, sep = "")
chrs <- names(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[1:23]

# Extract all nucleotide contexts' positions in specific chromosome
for (i in 1:nrow(nucleotide.selected)) {
  
  print(nucleotide.selected[i, ])
  precompute.motif.pos[[paste(nucleotide.selected[i, "type"], nucleotide.selected[i, "sequence"], sep = "")]] <- unique(c(BSgenome::start(Biostrings::matchPattern(nucleotide.selected[i, "sequence"], BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[chr.interest]])), BSgenome::start(Biostrings::matchPattern(as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(nucleotide.selected[i, "sequence"]))), BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[chr.interest]]))))
  
}
} else {
  
  precompute.motif.pos = NULL
  nucleotide.selected = NULL
  
}

# Read feature file paths
if (!is.null(continuous.features.selected.snv.url.file)) {
  
selected.continuous.urls <- read.delim(continuous.features.selected.snv.url.file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
continuous.selected.features = parallel::mclapply(selected.continuous.urls[ ,2], function(f) {
  print(f)
  df = read.delim(as.character(f), stringsAsFactors = FALSE, header = FALSE)
  with(df, GenomicRanges::GRanges(V1, IRanges::IRanges(V2, V3), score = V4))
  
}, mc.cores = cores)

names(continuous.selected.features) = as.character(selected.continuous.urls[ ,1])

} else {
  
  continuous.selected.features = NULL
  
}

if (!is.null(discrete.features.selected.snv.url.file)) {
  
selected.discrete.urls <- read.delim(discrete.features.selected.snv.url.file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
discrete.selected.features = parallel::mclapply(selected.discrete.urls[ ,2], function(f) {
  
  print(f)
  df = read.delim(as.character(f), stringsAsFactors = FALSE, header = FALSE)
  with(df, GenomicRanges::GRanges(V1, IRanges::IRanges(V2, V3)))
  
}, mc.cores = cores)

names(discrete.selected.features) = as.character(selected.discrete.urls[ ,1])

} else {
  
  discrete.selected.features = NULL
  
}

# Define SNV mutations 
maf.mutations <- maf.to.granges(snv.mutations.file)
maf.mutations = maf.mutations[as.character(GenomicRanges::seqnames(maf.mutations)) %in% seqnames]
mut.masked <- maf.mutations[S4Vectors::subjectHits(IRanges::findOverlaps(all.sites.masked, maf.mutations))]

# Define SNV sample mutation count based on full SNV mutations file
maf.mutations2 <- maf.to.granges(snv.mutations.file2)
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
sort(sapply(ls(), function(x) { object.size(get(x)) / 10 ^ 6 }))
rm(list = c("maf.mutations", "maf.ind", "mask.regions", "all.sites", "maf.mutations2"))
gc(reset = T)

# Tabulate covariates for mutations
GenomeInfoDb::seqlevels(mut.masked) = as.character(unique(GenomicRanges::seqnames(mut.masked)))
mut.chr = GenomicRanges::split(mut.masked, GenomicRanges::seqnames(mut.masked))
chrs <- names(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[1:23]
mut.chr = mut.chr[names(mut.chr) %in% chrs]

if (chr.interest %in% names(mut.chr)) {
  
mut.freq <- mutCovariate.snv.freq.table.muts(continuous.features = continuous.selected.features, discrete.features = discrete.selected.features, precompute.motif.pos = precompute.motif.pos, nucleotide.selected = nucleotide.selected, sample.specific.features = sample.specific.features, sites = mut.chr[[chr.interest]])
rm(list = c("mut.chr"))

} else {
  
  mut.freq = NULL
  
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

genome.freq <- parallel::mclapply(chunks, function(x) mutCovariate.snv.freq.table.genome(continuous.features = continuous.selected.features, discrete.features = discrete.selected.features, precompute.motif.pos = precompute.motif.pos, nucleotide.selected = nucleotide.selected, sites = all.sites.masked[x]), mc.cores = cores, mc.preschedule = FALSE, mc.silent = FALSE)
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
saveRDS(list(mut.freq, genome.freq.aggregated), file = paste("mutCovariate_", chr.interest, ".RDS", sep = ""))

}

