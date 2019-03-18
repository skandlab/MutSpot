#' Runs all or selected steps in MutSpot analysis.
#'
#' @param run.to Numeric vector defining which steps to run, default = 1, 2, 3.1, 4.1, 5.1, 5.2, 5.3, 5.4, 5.5, 6, 7, 8, 9.1, 9.2, 9.3.
#' @param chromosomes Character vector defining which chromosomes to compute feature matrix on, default = chr1-chrX.
#' @param snv.mutations SNV mutations MAF file, default = NULL.
#' @param indel.mutations Indel mutations MAF file, default = NULL.
#' @param mask.regions.file Regions to mask in genome, for example, non-mappable regions/immunoglobin loci/CDS regions RDS file, default file = mask_regions.RDS.
#' @param all.sites.file All sites in whole genome RDS file, default file = all_sites.RDS.
#' @param region.of.interest Region of interest bed file, default = NULL.
#' @param ratio Sampling ratio, default = 1.
#' @param sample To sample for non-mutated sites or to use all sites in region of interest, default = TRUE.
#' @param cores Number of cores, default = 1.
#' @param cutoff.nucleotide Frequency cutoff/threshold to determine nucleotide contexts used in prediction model, ranges from 0.5 to 1, default = 0.90.
#' @param cutoff.nucleotide.new Updated frequency cutoff/threshold to determine nucleotide contexts used in prediction model, ranges from 0.5 to 1.
#' @param genomic.features.snv Text file containing URLs of potential continuous and discrete SNV epigenetic features to select from, default = NULL.
#' @param genomic.features.indel Text file containing URLs of potential continuous and discrete indel epigenetic features to select from, default = NULL.
#' @param genomic.features Text file containing URLs of potential continuous and discrete SNV and indel epigenetic features to select from, default = NULL.
#' @param genomic.features.fixed.snv Text file containing URLs of fixed continuous and discrete SNV epigenetic features, default = NULL.
#' @param genomic.features.fixed.indel Text file containing URLs of fixed continuous and discrete indel epigenetic features, default = NULL.
#' @param genomic.features.fixed Text file containing URLs of fixed continuous and discrete SNV and indel epigenetic features, default = NULL.
#' @param sample.snv.features Text file containing sample specific SNV features, default = NULL.
#' @param sample.indel.features Text file containing sample specific indel features, default = NULL.
#' @param cutoff.features Frequency cutoff/threshold to determine epigenetic features used in prediction model, ranges from 0.5 to 1, default = 0.75.
#' @param cutoff.features.new.snv Updated frequency cutoff/threshold to determine SNV epigenetic features used in prediction model, ranges from 0.5 to 1.
#' @param cutoff.features.new.indel Updated frequency cutoff/threshold to determine indel epigenetic features used in prediction model, ranges from 0.5 to 1.
#' @param fit.sparse To fit model using glmnet or glm, default = FALSE.
#' @param drop To drop insignificant features from fitted model or not, default = FALSE.
#' @param min.count.snv Minimum number of mutated samples in each SNV hotspot, default = 2.
#' @param min.count.indel Minimum number of mutated samples in each indel hotspot, default = 2.
#' @param genome.size Total number of hotspots to run analysis on, default = 2533374732.
#' @param hotspots To run hotspot analysis or region-based analysis, default = TRUE.
#' @param promoter.file Promoter regions bed file, default file = Ensembl75.promoters.coding.bed
#' @param utr3.file 3'UTR regions bed file, default file = Ensembl75.3UTR.coding.bed
#' @param utr5.file 5'UTR regions bed file, default file = Ensembl75.5UTR.coding.bed
#' @param other.annotations Text file containing URLs of additional regions to be annotated, default = NULL.
#' @param fdr.cutoff FDR cutoff, default = 0.1.
#' @param color.line Color given FDR cutoff, default = red.
#' @param color.dots Color hotspots that passed given FDR cutoff, default = maroon1.
#' @param merge.hotspots To plot overlapping hotspots as 1 hotspot or individual hotspots, default = TRUE.
#' @param color.muts Color points, default = orange.
#' @param z.value To use z-value for plot or coefficients, default = FALSE.
#' @param top.no Number of top hotspots to plot, default = 3.
#' @param debug To keep intermediate output files or not, default = FALSE.
#' @return Corresponding output from each step in MutSpot analysis.
#' @export

MutSpot = function(run.to = c(1:2, 3.1, 4.1, 5.1, 5.2, 5.3, 5.4, 5.5, 6:8, 9.1, 9.2, 9.3), chromosomes = c(1:22,"X"), snv.mutations = NULL, indel.mutations = NULL, mask.regions.file = system.file("extdata", "mask_regions.RDS", package = "MutSpot"), all.sites.file = system.file("extdata", "all_sites.RDS", package = "MutSpot"), region.of.interest = NULL, ratio = 1, sample = T, cores = 1, cutoff.nucleotide = 0.90, cutoff.nucleotide.new = NULL, genomic.features.snv = NULL, genomic.features.indel = NULL, genomic.features = NULL, genomic.features.fixed.snv = NULL, genomic.features.fixed.indel = NULL, genomic.features.fixed = NULL, sample.snv.features = NULL, sample.indel.features = NULL, cutoff.features = 0.75, cutoff.features.new.snv = NULL, cutoff.features.new.indel = NULL, fit.sparse = FALSE, drop = FALSE, min.count.snv = 2, min.count.indel = 2, genome.size = 2533374732, hotspots = TRUE,
                  promoter.file = system.file("extdata", "Ensembl75.promoters.coding.bed", package = "MutSpot"), utr3.file = system.file("extdata", "Ensembl75.3UTR.coding.bed", package = "MutSpot"), utr5.file = system.file("extdata", "Ensembl75.5UTR.coding.bed", package = "MutSpot"), other.annotations = NULL, fdr.cutoff = 0.1, color.line = "red", color.dots = "maroon1", merge.hotspots = TRUE, color.muts = "orange", z.value = FALSE, top.no = 3, debug = FALSE) {
  
## check format of working directory ##
if (substr(working.dir, nchar(working.dir), nchar(working.dir)) != "/") {
  
  working.dir = paste(working.dir, "/", sep = "")
  
}  
  setwd(working.dir)
  output.dir = paste(working.dir, "results/", sep = "")
  dir.create(output.dir)
  
  ## Load all dependencies ##
  if("data.table" %in% rownames(installed.packages()) == FALSE) {
    
    print("install [data.table]")
    install.packages("data.table")
    
  }
  if ("BSgenome.Hsapiens.UCSC.hg19" %in% rownames(installed.packages()) == FALSE) {
    
    print("install [BSgenome.Hsapiens.UCSC.hg19]")
    source("https://bioconductor.org/biocLite.R")
    biocLite("BSgenome.Hsapiens.UCSC.hg19")
    
  }
  if("Biostrings" %in% rownames(installed.packages()) == FALSE) {
    
    print("install [Biostrings]")
    source("https://bioconductor.org/biocLite.R")
    biocLite("Biostrings")
    
  }
  if("GenomicRanges" %in% rownames(installed.packages()) == FALSE) {
    
    print("install [GenomicRanges]")
    source("https://bioconductor.org/biocLite.R")
    biocLite("GenomicRanges")
    
  }
  if("ggplot2" %in% rownames(installed.packages()) == FALSE) {
    
    print("install [ggplot2]")
    install.packages("ggplot2")
    
  }
  if("glmnet" %in% rownames(installed.packages()) == FALSE) {
    
    print("install [glmnet]")
    install.packages("glmnet")
    
  }
  if("Matrix" %in% rownames(installed.packages()) == FALSE) {
    
    print("install [Matrix]")
    install.packages("Matrix")
    
  }
  if("poibin" %in% rownames(installed.packages()) == FALSE) {
    
    print("install [poibin]")
    install.packages("poibin")
    
  }
  if("stringi" %in% rownames(installed.packages()) == FALSE) {
    
    print("install [stringi]")
    install.packages("stringi")
    
  }
  if("stringr" %in% rownames(installed.packages()) == FALSE) {
    
    print("install [stringr]")
    install.packages("stringr")
    
  }
  if("plyr" %in% rownames(installed.packages()) == FALSE) {
    
    print("install [plyr]")
    install.packages("plyr")
    
  }
  if("parallel" %in% rownames(installed.packages()) == FALSE) {
    
    print("install [parallel]")
    install.packages("parallel")
    
  }
  
  ## Load all required packages ##
  suppressWarnings(suppressMessages(library(Biostrings)))
  suppressWarnings(suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19)))
  suppressWarnings(suppressMessages(library(GenomicRanges)))
  suppressWarnings(suppressMessages(library(data.table)))
  suppressWarnings(suppressMessages(library(ggplot2)))
  suppressWarnings(suppressMessages(library(glmnet)))
  suppressWarnings(suppressMessages(library(Matrix)))
  suppressWarnings(suppressMessages(library(parallel)))
  suppressWarnings(suppressMessages(library(plyr)))
  suppressWarnings(suppressMessages(library(stringr)))
  suppressWarnings(suppressMessages(library(stringi)))
  suppressWarnings(suppressMessages(library(poibin)))
  
## Step 1 ##
sampled.sites.snv.file = paste(output.dir, "sampled.sites.snv.RDS", sep = "")
snv.mutations.region.file = paste(output.dir, "SNV_region.MAF", sep = "")
sampled.sites.indel.file = paste(output.dir, "sampled.sites.indel.RDS", sep = "")
indel.mutations.region.file = paste(output.dir, "indel_region.MAF", sep = "")
  
if (1 %in% run.to) {
    
  print("Sample sites from SNV/indel")
  
sample_sites = sample.sites(snv.mutations.file = snv.mutations, indel.mutations.file = indel.mutations, mask.regions.file = mask.regions.file, all.sites.file = all.sites.file, region.of.interest = region.of.interest, ratio = ratio, sample = sample, cores = cores)

# Save SNV mutations in specified region as MAF file
if (!is.null(sample_sites[[1]])) {
  
  if (nrow(sample_sites[[1]]) != 0) {
  
  write.table(sample_sites[[1]], file = snv.mutations.region.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
  }
  
}

# Save SNV mutated and non-mutated sites as RDS file
if (!is.null(sample_sites[[2]])) {
  
  saveRDS(sample_sites[[2]], sampled.sites.snv.file)
  
}

# Save indel mutations in specified region as MAF file
if (!is.null(sample_sites[[3]])) {
  
  if (nrow(sample_sites[[3]]) != 0) {
    
  write.table(sample_sites[[3]], file = indel.mutations.region.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
  }
  
}


# Save indel mutated and non-mutated sites as RDS file
if (!is.null(sample_sites[[4]])) {
  
  saveRDS(sample_sites[[4]], sampled.sites.indel.file)
  
}

}

if (!is.null(region.of.interest) & !is.null(snv.mutations)) {
  
  snv.mutations.int = snv.mutations.region.file
  
} else {
  
  snv.mutations.int = snv.mutations
  
}

if (!is.null(region.of.interest) & !is.null(indel.mutations)) {
  
  indel.mutations.int = indel.mutations.region.file
  
} else {
  
  indel.mutations.int = indel.mutations
  
}

if (is.null(snv.mutations) | !file.exists(sampled.sites.snv.file)) {
  
  sampled.sites.snv.file = NULL
  
}

if (is.null(indel.mutations) | !file.exists(sampled.sites.indel.file)) {
  
  sampled.sites.indel.file = NULL
  
}


## Step 2 ##
local.mutrate.snv.file = paste(working.dir, "features/", "localmutrate_snv.bed", sep = "")
local.mutrate.indel.file = paste(working.dir, "features/", "localmutrate_indel.bed", sep = "")

  if (2 %in% run.to) {
    
    print("Calculate local mutation rates for SNV/indel")
    
local_mutrate = local.mutrate(snv.mutations.file = snv.mutations, indel.mutations.file = indel.mutations)

# Save binned SNV local mutation rate as bed file
if (!is.null(local_mutrate[[1]])) {
  
  write.table(local_mutrate[[1]], file = local.mutrate.snv.file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
}

# Save binned indel local mutation rate as bed file
if (!is.null(local_mutrate[[2]])) {
  
  write.table(local_mutrate[[2]], file = local.mutrate.indel.file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
}

}

## Step 3 ##
nucleotide.stabs.freq.file = paste(output.dir, "nucleotide_stabs_freq.RDS", sep = "")
nucleotide.selected.file = paste(output.dir, "nucleotide_selected.RDS", sep = "")
indels.polyAT.file = paste(output.dir, "indel_polyAT.MAF", sep = "")

## Step 3I ##
if (3.1 %in% run.to) {
  
  print("Select nucleotide context for SNV and filter polyAT from indel")
    
nucleotide_selection = nucleotide.selection(sampled.sites.snv.file = sampled.sites.snv.file, indel.mutations.file = indel.mutations.int, cutoff = cutoff.nucleotide, cores = cores)

# Check for need to change threshold
print(head(nucleotide_selection[[1]][order(nucleotide_selection[[1]]$f, decreasing=TRUE), ]))

# Save lasso selection frequency table as RDS file
if (!is.null(nucleotide_selection[[1]])) {
  
  saveRDS(nucleotide_selection[[1]], file = nucleotide.stabs.freq.file)
  
}

# Save selected nucleotide contexts as RDS file
if (!is.null(nucleotide_selection[[2]]) & length(nucleotide_selection[[2]] != 0)) {
  
  saveRDS(nucleotide_selection[[2]], file = nucleotide.selected.file)
  
}

# Save filtered indel mutations as MAF file
if (!is.null(nucleotide_selection[[3]])) {
  
  if (nrow(nucleotide_selection[[3]]) != 0) {
  
  write.table(nucleotide_selection[[3]], file = indels.polyAT.file, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
  indel.mutations = indels.polyAT.file
  
  }
  
}

}

## Step 3II ##
if (3.2 %in% run.to) {
  
  print("Reselect nucleotide context for SNV using altered threshold")
  
if (!is.null(cutoff.nucleotide.new)) {
  
nucleotide_selection_adjust = nucleotide.selection.adjust(nucleotide.stabs.file = nucleotide.stabs.freq.file, new.cutoff = cutoff.nucleotide.new)

saveRDS(nucleotide_selection_adjust, file = nucleotide.selected.file)

}

}

## Step 4 ##
features.stabs.snv.file = paste(output.dir, "features_stabs_snv.RDS", sep = "")
continuous.features.selected.snv.url.file = paste(output.dir, "continuous_features_selected_snv_url.txt", sep = "")
discrete.features.selected.snv.url.file = paste(output.dir, "discrete_features_selected_snv_url.txt", sep = "")
features.stabs.indel.file = paste(output.dir, "features_stabs_indel.RDS", sep = "")
continuous.features.selected.indel.url.file = paste(output.dir, "continuous_features_selected_indel_url.txt", sep = "")
discrete.features.selected.indel.url.file = paste(output.dir, "discrete_features_selected_indel_url.txt", sep = "")

if (4.1 %in% run.to) {
  
  print("Select epigenetic features for SNV/indel")
  
    ## Step 4I ##
epigenetic_selection = epigenetic.selection(sampled.sites.snv.file = sampled.sites.snv.file, sampled.sites.indel.file = sampled.sites.indel.file, genomic.features.snv = genomic.features.snv, genomic.features.indel = genomic.features.indel, genomic.features = genomic.features, genomic.features.fixed.snv = genomic.features.fixed.snv, genomic.features.fixed.indel = genomic.features.fixed.indel, genomic.features.fixed = genomic.features.fixed, cores = cores, cutoff = cutoff.features)

# Check for need to change threshold for SNVs
if (!is.null(epigenetic_selection[[1]])) {
  
print(head(epigenetic_selection[[1]][order(epigenetic_selection[[1]]$f, decreasing = TRUE), ]))

}

# Save lasso selection frequency table as RDS file 
if (!is.null(epigenetic_selection[[1]])) {
  
  saveRDS(epigenetic_selection[[1]], file = features.stabs.snv.file)
  
}

# Save URLs of selected SNV continuous features as text file
if (!is.null(epigenetic_selection[[2]])) {
  
  if (nrow(epigenetic_selection[[2]]) != 0) {
  
  write.table(epigenetic_selection[[2]], file = continuous.features.selected.snv.url.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  }
  
}

# Save URLs of selected SNV discrete features as text file
if (!is.null(epigenetic_selection[[3]])) {
  
  if (nrow(epigenetic_selection[[3]]) != 0) {
  
  write.table(epigenetic_selection[[3]], file = discrete.features.selected.snv.url.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  }

}

# Check for need to change threshold for indels
if (!is.null(epigenetic_selection[[4]])) {
  
print(head(epigenetic_selection[[4]][order(epigenetic_selection[[4]]$f, decreasing = TRUE), ]))

}

# Save lasso selection frequency table as RDS file
if (!is.null(epigenetic_selection[[4]])) {
  
  saveRDS(epigenetic_selection[[4]], file = features.stabs.indel.file)
  
}

# Save URLs of selected indel continuous features as text file
if (!is.null(epigenetic_selection[[5]])) {
  
  if (nrow(epigenetic_selection[[5]]) != 0) {
  
  write.table(epigenetic_selection[[5]], file = continuous.features.selected.indel.url.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  }
  
}

# Save URLs of selected indel discrete features as text file
if (!is.null(epigenetic_selection[[6]])) {
  
  if (nrow(epigenetic_selection[[6]]) != 0) {
  
  write.table(epigenetic_selection[[6]], file = discrete.features.selected.indel.url.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  }
  
}

}

## Step 4II ##
if (4.2 %in% run.to) {
  
  print("Reselect epigenetic features for SNV/indel using altered threshold")
  
if (!is.null(cutoff.features.new.snv) | !is.null(cutoff.features.new.indel)) {

  epigenetic_selection_adjust = epigenetic.selection.adjust(feature.stabs.snv.file = features.stabs.snv.file , feature.stabs.indel.file = features.stabs.indel.file, continuous.features.selected.snv.url.file = continuous.features.selected.snv.url.file, discrete.features.selected.snv.url.file = discrete.features.selected.snv.url.file, continuous.features.selected.indel.url.file = continuous.features.selected.indel.url.file, discrete.features.selected.indel.url.file = discrete.features.selected.indel.url.file, new.cutoff.snv = cutoff.features.new.snv, new.cutoff.indel = cutoff.features.new.indel)
  
  # Save URLs of selected SNV continuous features as text file
  if (!is.null(epigenetic_selection_adjust[[1]])) {
    
    write.table(epigenetic_selection_adjust[[1]], file = continuous.features.selected.snv.url.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
  }
  
  # Save URLs of selected SNV discrete features as text file
  if (!is.null(epigenetic_selection_adjust[[2]])) {
    
    write.table(epigenetic_selection_adjust[[2]], file = discrete.features.selected.snv.url.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
  }
  
  # Save URLs of selected indel continuous features as text file
  if (!is.null(epigenetic_selection_adjust[[3]])) {
    
    write.table(epigenetic_selection_adjust[[3]], file = continuous.features.selected.indel.url.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
  }
  
  # Save URLs of selected indel discrete features as text file
  if (!is.null(epigenetic_selection_adjust[[4]])) {
    
    write.table(epigenetic_selection_adjust[[4]], file = discrete.features.selected.indel.url.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
  }
  
}

  }
  
## Step 5 ##
mutCovariate.snv.output.file = paste(output.dir, "mutCovariate-compile-part1.RDS", sep = "")
mutCovariate.snv.output.p1 = paste(output.dir, "mutCovariate-sparse-p1.RDS", sep = "")
mutCovariate.snv.output.p2 = paste(output.dir, "mutCovariate-sparse-p2.RDS", sep = "")
mutCovariate.indel.output.file = paste(output.dir, "mutCovariate-indel-compile-part1.RDS", sep = "")
mutCovariate.indel.output.p1 = paste(output.dir, "indel-mutCovariate-sparse-p1.RDS", sep = "")
mutCovariate.indel.output.p2 = paste(output.dir, "indel-mutCovariate-sparse-p2.RDS", sep = "")

if (file.exists(indels.polyAT.file)) {
  
  indel.mutations.int = indels.polyAT.file
  
}

if (!file.exists(continuous.features.selected.snv.url.file)) {
  
  continuous.features.selected.snv.url.file = NULL
  
}
if (!file.exists(discrete.features.selected.snv.url.file)) {
  
  discrete.features.selected.snv.url.file = NULL
  
}

if (!file.exists(continuous.features.selected.indel.url.file)) {
  
  continuous.features.selected.indel.url.file = NULL
  
}
if (!file.exists(discrete.features.selected.indel.url.file)) {
  
  discrete.features.selected.indel.url.file = NULL
  
}
if (!file.exists(nucleotide.selected.file)) {
  
  nucleotide.selected.file = NULL
  
}

if (5.1 %in% run.to) {
  
  print("Prepare feature matrix for SNV/indel for each chromosome")
  
    ## Step 5a ##
  if ((!is.null(continuous.features.selected.snv.url.file) | !is.null(discrete.features.selected.snv.url.file) | !is.null(sample.snv.features) | !is.null(nucleotide.selected.file)) & !is.null(sampled.sites.snv.file)) {
    
    sample.specific.features.url.file = sample.snv.features
    snv.mutations.file = snv.mutations.int
    snv.mutations.file2 = snv.mutations
    
    # Chr1-X
    chrOrder <- c(paste("chr", 1:22, sep = ""), "chrX")
    seqi = GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[GenomeInfoDb::seqnames(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[1:23]]
    seqnames = GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens))[1:23]
    
    # Define masked region i.e. CDS, immunoglobulin loci and nonmappable
    mask.regions = readRDS(mask.regions.file)
    mask.regions = mask.regions[as.character(GenomeInfoDb::seqnames(mask.regions)) %in% seqnames]
    
    # Define all sites in whole genome
    all.sites = readRDS(all.sites.file)
    all.sites = all.sites[as.character(GenomeInfoDb::seqnames(all.sites)) %in% seqnames]
    all.sites.masked = subtract.regions.from.roi(all.sites, mask.regions, cores = cores)
    # sum(as.numeric(GenomicRanges::width(all.sites.masked)))
    
    # If specified region, redefine all sites to be in specified region
    if (!is.null(region.of.interest)) {
      
      print("specified region")
      
      all.sites = read.delim(region.of.interest, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      all.sites = with(all.sites, GenomicRanges::GRanges(V1, IRanges::IRanges(V2, V3)))
      all.sites.masked = subtract.regions.from.roi(all.sites,mask.regions, cores = cores)
      all.sites.masked = all.sites.masked[as.character(GenomeInfoDb::seqnames(all.sites.masked)) %in% seqnames]
      
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
    maf.mutations = maf.mutations[as.character(GenomeInfoDb::seqnames(maf.mutations)) %in% seqnames]
    mut.masked <- maf.mutations[S4Vectors::subjectHits(IRanges::findOverlaps(all.sites.masked, maf.mutations))]
    
    # Define SNV sample mutation count based on full SNV mutations file
    maf.mutations2 <- maf.to.granges(snv.mutations.file2)
    maf.mutations2 = maf.mutations2[as.character(GenomeInfoDb::seqnames(maf.mutations2)) %in% seqnames]
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
          if (class(t) == "matrix") {
            
          colnames(t)=substr(colnames(t),2,nchar(colnames(t)))
          colnames(t)=paste(colnames(sample.specific.features)[x], colnames(t), sep = "")
          rownames(t)=rownames(sample.specific.features)
          
          } else {
            
            t = as.data.frame(t)
            colnames(t) = paste(colnames(sample.specific.features)[x], levels(factor(sample.specific.features[,x]))[2], sep = "")
            rownames(t) = rownames(sample.specific.features)
            
          }
          
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
    GenomeInfoDb::seqlevels(mut.masked) = as.character(unique(GenomeInfoDb::seqnames(mut.masked)))
    mut.chr = GenomicRanges::split(mut.masked, GenomeInfoDb::seqnames(mut.masked))
    chrs <- names(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[1:23]
    mut.chr = mut.chr[names(mut.chr) %in% chrs]
    
    
    for (chr.interest in paste("chr", chromosomes, sep="")) {
      
      
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
          
          print(paste(chr.interest,nucleotide.selected[i, "sequence"], nucleotide.selected[i,"type"],sep=":"))
          if (nucleotide.selected$type[i] %in% c("oneMer","threeMer","fiveMer")){
            precompute.motif.pos[[paste(nucleotide.selected[i, "type"], nucleotide.selected[i, "sequence"], sep = "")]] <- unique(c(BSgenome::start(Biostrings::matchPattern(nucleotide.selected[i, "sequence"], BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[chr.interest]])), BSgenome::start(Biostrings::matchPattern(as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(nucleotide.selected[i, "sequence"]))), BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[chr.interest]]))))
          } else if (nucleotide.selected$type[i] == "threeRight") {
            a <- unique(BSgenome::start(Biostrings::matchPattern(nucleotide.selected[i, "sequence"], BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[chr.interest]])))
            b <- unique(BSgenome::start(Biostrings::matchPattern(as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(nucleotide.selected[i, "sequence"]))), BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[chr.interest]])))+1
            precompute.motif.pos[[paste(nucleotide.selected[i, "type"], nucleotide.selected[i, "sequence"], sep = "")]]=unique(c(a,b))
          } else if (nucleotide.selected$type[i]=="threeLeft") {
            a <- unique(BSgenome::start(Biostrings::matchPattern(nucleotide.selected[i, "sequence"], BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[chr.interest]])))+1
            b <- unique(BSgenome::start(Biostrings::matchPattern(as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(nucleotide.selected[i, "sequence"]))), BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[chr.interest]])))
            precompute.motif.pos[[paste(nucleotide.selected[i, "type"], nucleotide.selected[i, "sequence"], sep = "")]]=unique(c(a,b))
          }  else if (nucleotide.selected$type[i] =="fiveRight") {
            a <- unique(BSgenome::start(Biostrings::matchPattern(nucleotide.selected[i, "sequence"], BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[chr.interest]])))
            b <- unique(BSgenome::start(Biostrings::matchPattern(as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(nucleotide.selected[i, "sequence"]))), BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[chr.interest]])))+2
            precompute.motif.pos[[paste(nucleotide.selected[i, "type"], nucleotide.selected[i, "sequence"], sep = "")]]=unique(c(a,b))
          } else if (nucleotide.selected$type[i]=="fiveLeft") {
            a <- unique(BSgenome::start(Biostrings::matchPattern(nucleotide.selected[i, "sequence"], BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[chr.interest]])))+2
            b <- unique(BSgenome::start(Biostrings::matchPattern(as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(nucleotide.selected[i, "sequence"]))), BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[chr.interest]])))
            precompute.motif.pos[[paste(nucleotide.selected[i, "type"], nucleotide.selected[i, "sequence"], sep = "")]]=unique(c(a,b))
          }
          
        }
      } else {
        
        precompute.motif.pos = NULL
        nucleotide.selected = NULL
        
      }
      
      if (chr.interest %in% names(mut.chr)) {
        
        mut.freq <- mutCovariate.snv.freq.table.muts(continuous.features = continuous.selected.features, discrete.features = discrete.selected.features, precompute.motif.pos = precompute.motif.pos, nucleotide.selected = nucleotide.selected, sample.specific.features = sample.specific.features, sites = mut.chr[[chr.interest]])

      } else {
        
        mut.freq = NULL
        
      }
      
      if (chr.interest %in% as.character(GenomeInfoDb::seqnames(all.sites.masked))) {
        
        all.sites.masked2 = all.sites.masked[as.character(GenomeInfoDb::seqnames(all.sites.masked)) == chr.interest]
        len = sapply(GenomicRanges::split(all.sites.masked2, GenomeInfoDb::seqnames(all.sites.masked2)), length)
        len = len[len != 0]
        len2 = sapply(1:length(len), function(i) { sum(len[1:i]) })
        len2 = c(0, len2)
        chunk.size = 5000
        chunks <- lapply(1:length(len), function(j) { lapply(1:ceiling(len[j] / chunk.size), function(i) ((len2[j] + (i-1) * chunk.size + 1):(len2[j] + min((i) * chunk.size, len[j])))) })
        chunks = do.call(c, chunks)
        
        genome.freq <- parallel::mclapply(chunks, function(x) mutCovariate.snv.freq.table.genome(continuous.features = continuous.selected.features, discrete.features = discrete.selected.features, precompute.motif.pos = precompute.motif.pos, nucleotide.selected = nucleotide.selected, sites = all.sites.masked2[x]), mc.cores = cores, mc.preschedule = FALSE, mc.silent = FALSE)
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
      
      saveRDS(list(mut.freq, genome.freq.aggregated), file = paste(output.dir,"mutCovariate_", chr.interest, ".RDS", sep = ""))
      
    }
    
  }
  
  ## Step 5b ##
  if ((!is.null(continuous.features.selected.indel.url.file) | !is.null(discrete.features.selected.indel.url.file) | !is.null(sample.indel.features)) & !is.null(sampled.sites.indel.file)) {
    
    sample.specific.features.url.file = sample.indel.features
    indel.mutations.file = indel.mutations.int
    indel.mutations.file2 = indel.mutations
    
    chrOrder <- c(paste("chr", 1:22, sep = ""), "chrX")
    seqi = GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[GenomeInfoDb::seqnames(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[1:23]]
    seqnames = GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens))[1:23]
    
    # Define masked region i.e. CDS, immunoglobulin loci and nonmappable
    mask.regions = readRDS(mask.regions.file)
    mask.regions = mask.regions[as.character(GenomeInfoDb::seqnames(mask.regions)) %in% seqnames]
    
    # Define all sites in whole genome
    all.sites = readRDS(all.sites.file)
    all.sites = all.sites[as.character(GenomeInfoDb::seqnames(all.sites)) %in% seqnames]
    all.sites.masked = subtract.regions.from.roi(all.sites, mask.regions, cores = cores)
    sum(as.numeric(GenomicRanges::width(all.sites.masked)))
    
    # If specified region, redefine all sites to be in specified region
    if (!is.null(region.of.interest)) {
      
      print("specified region")
      
      all.sites = read.delim(region.of.interest, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      all.sites = with(all.sites, GenomicRanges::GRanges(V1, IRanges::IRanges(V2, V3)))
      all.sites.masked = subtract.regions.from.roi(all.sites, mask.regions, cores = cores)
      all.sites.masked = all.sites.masked[as.character(GenomeInfoDb::seqnames(all.sites.masked)) %in% seqnames]
      
    }
    
    # Extract all polyA, poly C, poly G and poly T positions in a specific chromosome
    chrs <- names(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[1:23]
    
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
    maf.mutations = maf.mutations[as.character(GenomeInfoDb::seqnames(maf.mutations)) %in% seqnames]
    mut.masked <- maf.mutations[S4Vectors::subjectHits(IRanges::findOverlaps(all.sites.masked, maf.mutations))]
    mut.masked.sites = mut.masked
    GenomicRanges::start(mut.masked.sites) = GenomicRanges::start(mut.masked.sites) + ceiling((GenomicRanges::width(mut.masked.sites) - 1) / 2)
    GenomicRanges::end(mut.masked.sites) = GenomicRanges::start(mut.masked.sites)
    
    # Define indel sample mutation count based on full indel mutations file
    maf.mutations2 = maf.to.granges(indel.mutations.file2)
    maf.mutations2 = maf.mutations2[as.character(GenomeInfoDb::seqnames(maf.mutations2)) %in% seqnames]
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
          if (class(t) == "matrix") {
            
          colnames(t)=substr(colnames(t),2,nchar(colnames(t)))
          colnames(t)=paste(colnames(sample.specific.features)[x],colnames(t),sep="")
          rownames(t)=rownames(sample.specific.features)
          
          } else {
            
            t = as.data.frame(t)
            colnames(t) = paste(colnames(sample.specific.features)[x], levels(factor(sample.specific.features[,x]))[2], sep = "")
            rownames(t) = rownames(sample.specific.features)
            
          }
          
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
    GenomeInfoDb::seqlevels(mut.masked.sites) = as.character(unique(GenomeInfoDb::seqnames(mut.masked.sites)))
    mut.chr = GenomicRanges::split(mut.masked.sites, GenomeInfoDb::seqnames(mut.masked.sites))
    mut.chr = mut.chr[names(mut.chr) %in% chrs]

    # Tabulate covariates for all positions in indels
    GenomeInfoDb::seqlevels(mut.masked) = as.character(unique(GenomeInfoDb::seqnames(mut.masked)))
    mut.indel = mut.masked
    mut.indel.chr = GenomicRanges::split(mut.indel, GenomeInfoDb::seqnames(mut.indel))
    mut.indel.chr = mut.indel.chr[names(mut.indel.chr) %in% chrs]
    
    for (chr.interest in paste("chr", chromosomes, sep="")) {
      
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
      
      if (chr.interest %in% names(mut.chr)) {
        
        mut.freq <- mutCovariate.indel.freq.table.muts(continuous.features = continuous.selected.features, discrete.features = discrete.selected.features, sample.specific.features = sample.specific.features, polyAs = polyAs, polyTs = polyTs, polyCs = polyCs, polyGs = polyGs, sites = mut.chr[[chr.interest]])
        # rm(list = c("mut.chr"))
        
      } else {
        
        mut.freq = NULL
        
      }
      
      if (chr.interest %in% names(mut.indel.chr)) {
        
        indel.freq <- mutCovariate.indel.freq.table.muts(continuous.features = continuous.selected.features, discrete.features = discrete.selected.features, sample.specific.features = sample.specific.features, polyAs = polyAs, polyTs = polyTs, polyCs = polyCs, polyGs = polyGs, sites = mut.indel.chr[[chr.interest]])
        # rm(list = c("mut.indel.chr"))
        
      } else {
        
        indel.freq = NULL
        
      }
      
      if (chr.interest %in% as.character(GenomeInfoDb::seqnames(all.sites.masked))) {
        
        all.sites.masked2 = all.sites.masked[as.character(GenomeInfoDb::seqnames(all.sites.masked)) == chr.interest]
        len = sapply(GenomicRanges::split(all.sites.masked2, GenomeInfoDb::seqnames(all.sites.masked2)), length)
        len = len[len != 0]
        len2 = sapply(1:length(len), function(i) { sum(len[1:i]) })
        len2 = c(0, len2)
        chunk.size = 5000
        chunks <- lapply(1:length(len), function(j) { lapply(1:ceiling(len[j] / chunk.size), function(i) ((len2[j] + (i-1) * chunk.size + 1):(len2[j] + min((i) * chunk.size, len[j])))) })
        chunks = do.call(c, chunks)
        
        genome.freq <- parallel::mclapply(chunks, function(x) mutCovariate.indel.freq.table.genome(continuous.features = continuous.selected.features, discrete.features = discrete.selected.features, polyAs = polyAs, polyTs = polyTs, polyCs = polyCs, polyGs = polyGs, sites = all.sites.masked2[x]), mc.cores = cores, mc.preschedule = FALSE, mc.silent = FALSE)
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
      
      saveRDS(list(mut.freq, indel.freq, genome.freq.aggregated), file = paste(output.dir, "mutCovariate_indel_", chr.interest, ".RDS", sep = ""))
      
    }
    
  }
  
}
  
if (5.2 %in% run.to) {
  
  print("Compile feature matrix for snv for all chromosomes")
  
  if ((!is.null(continuous.features.selected.snv.url.file) | !is.null(discrete.features.selected.snv.url.file) | !is.null(sample.snv.features)) & !is.null(sampled.sites.snv.file)) {
  
mutCovariate_snv_compile = mutCovariate.snv.compile(mask.regions.file = mask.regions.file, all.sites.file = all.sites.file, snv.mutations.file = snv.mutations.int, sample.specific.features.url.file = sample.snv.features, region.of.interest = region.of.interest, cores = cores, snv.mutations.file2 = snv.mutations, chrom.dir = output.dir)

saveRDS(mutCovariate_snv_compile, file = mutCovariate.snv.output.file)

delete.files = Sys.glob(paste(output.dir, "mutCovariate_chr*.RDS", sep = ""))
for (i in delete.files) {
  
  unlink(i)
  
}

  }
  
}

if (5.3 %in% run.to) {
  
  print("Convert feature matrix for snv to sparse matrix")
  
  if (file.exists(mutCovariate.snv.output.file)) {
    
mutCovariate_snv_sparse = mutCovariate.snv.sparse(compiled = mutCovariate.snv.output.file)

saveRDS(mutCovariate_snv_sparse[[1]], file = mutCovariate.snv.output.p1)
saveRDS(mutCovariate_snv_sparse[[2]], file = mutCovariate.snv.output.p2)
  
unlink(mutCovariate.snv.output.file)

  }
  
}

if (5.4 %in% run.to) {
  
  print("Compile feature matrix for indel for all chromosomes")
  
  if ((!is.null(continuous.features.selected.indel.url.file) | !is.null(discrete.features.selected.indel.url.file) | !is.null(sample.indel.features)) & !is.null(sampled.sites.indel.file)) {
    
mutCovariate_indel_compile = mutCovariate.indel.compile(mask.regions.file = mask.regions.file, all.sites.file = all.sites.file, indel.mutations.file = indel.mutations.int, sample.specific.features.url.file = sample.indel.features, region.of.interest = region.of.interest, cores = cores, indel.mutations.file2 = indel.mutations)

saveRDS(mutCovariate_indel_compile, file = mutCovariate.indel.output.file)

delete.files = Sys.glob(paste(output.dir, "mutCovariate_indel_chr*.RDS", sep = ""))
for (i in delete.files) {
  
  unlink(i)
  
}

}

}

if (5.5 %in% run.to) {
  
  print("Convert feature matrix for indel to sparse matrix")
  
  if (file.exists(mutCovariate.indel.output.file)) {
  
mutCovariate_indel_sparse = mutCovariate.indel.sparse(compiled = mutCovariate.indel.output.file)

saveRDS(mutCovariate_indel_sparse[[1]], file = mutCovariate.indel.output.p1)
saveRDS(mutCovariate_indel_sparse[[2]], file = mutCovariate.indel.output.p2)
  
unlink(mutCovariate.indel.output.file)

  }
  
}

## Step 6 ##
LRmodel.snv.file = paste(output.dir, "snv-LRmodel", sep = "")
LRmodel.indel.file = paste(output.dir, "indel-LRmodel", sep = "")

if (6 %in% run.to) {
  
  print("Fit prediction model for SNV/indel")
  
  ## Step 6a ##
  if (file.exists(mutCovariate.snv.output.p1)) {
    
  LRmodel = mutLRFit.snv(mutCovariate.table.snv.file = mutCovariate.snv.output.p1, mutCovariate.count.snv.file = mutCovariate.snv.output.p2, continuous.features.selected.snv.url.file = continuous.features.selected.snv.url.file, discrete.features.selected.snv.url.file = discrete.features.selected.snv.url.file, nucleotide.selected.file = nucleotide.selected.file, sample.specific.features.url.file = sample.snv.features, fit.sparse = fit.sparse, drop = drop)

  if (class(LRmodel)[1]!="list") {
    
save(LRmodel, file = LRmodel.snv.file)
  } else {
    if (!"unchanged" %in% LRmodel[[2]]) {
      saveRDS(LRmodel[[2]],file=nucleotide.selected.file)
    }
    if (!"unchanged" %in% LRmodel[[3]]){
      if(!is.null(LRmodel[[3]])) {
      write.table(LRmodel[[3]],file=continuous.features.selected.snv.url.file,sep = "\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
      } else {
        z <- data.frame(V1=character(),V2=character(),stringsAsFactors=FALSE)
        write.table(z,file=continuous.features.selected.snv.url.file,sep = "\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
      }
    }
    if (!"unchanged" %in% LRmodel[[4]]){
      if(!is.null(LRmodel[[4]])) {
      write.table(LRmodel[[4]],file=discrete.features.selected.snv.url.file,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
      } else {
        z <- data.frame(V1=character(),V2=character(),stringsAsFactors=FALSE)
        write.table(z,file=discrete.features.selected.snv.url.file,sep = "\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
      }
    }
    if (!"unchanged" %in% LRmodel[[5]]){
      if(!is.null(LRmodel[[5]])) {
      write.table(LRmodel[[5]], file=sample.snv.features,sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
      } else {
        z <- data.frame(V1=character(),V2=character(),stringsAsFactors=FALSE)
        write.table(z,file=sample.snv.features,sep = "\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
      }
    }
    
    LRmodel=LRmodel[[1]]
    save(LRmodel,file=LRmodel.snv.file)
    
  }

  }
  
## Step 6b ##
  if (file.exists(mutCovariate.indel.output.p1)) {
    
LRmodel = mutLRFit.indel(mutCovariate.table.indel.file = mutCovariate.indel.output.p1, mutCovariate.count.indel.file = mutCovariate.indel.output.p2, continuous.features.selected.indel.url.file = continuous.features.selected.indel.url.file, discrete.features.selected.indel.url.file = discrete.features.selected.indel.url.file, sample.specific.features.url.file = sample.indel.features, fit.sparse = fit.sparse, drop = drop)

if (class(LRmodel)[1]!="list") {
  
  save(LRmodel, file = LRmodel.indel.file)
} else {
  if (!"unchanged" %in% LRmodel[[2]]){
    if (!is.null(LRmodel[[2]])) {
    write.table(LRmodel[[2]],file=continuous.features.selected.indel.url.file,sep = "\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
    } else {
      z <- data.frame(V1=character(),V2=character(),stringsAsFactors=FALSE)
      write.table(z,file=continuous.features.selected.indel.url.file,sep = "\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
    }
  }
  if (!"unchanged" %in% LRmodel[[3]]){
    if (!is.null(LRmodel[[3]])) {
    write.table(LRmodel[[3]],file=discrete.features.selected.indel.url.file,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
    } else {
      z <- data.frame(V1=character(),V2=character(),stringsAsFactors=FALSE)
      write.table(z,file=discrete.features.selected.indel.url.file,sep = "\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
    }
  }
  if (!"unchanged" %in% LRmodel[[4]]){
    if (!is.null(LRmodel[[4]])) {
    write.table(LRmodel[[4]], file=sample.indel.features,sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
    } else {
      z <- data.frame(V1=character(),V2=character(),stringsAsFactors=FALSE)
      write.table(z,file=sample.indel.features,sep = "\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
    }
  }
  
  LRmodel=LRmodel[[1]]
  save(LRmodel,file=LRmodel.indel.file)
  
}

  }
  
}
  
## Step 7 ##
snv.hotspots = paste(output.dir, "snv_hotspots.tsv", sep = "")
snv.hotspots.merged = paste(output.dir,"snv_hotspots_merged.tsv",sep="")
indel.hotspots = paste(output.dir, "indel_hotspots.tsv", sep = "")
indel.hotspots.merged=paste(output.dir,"indel_hotspots_merged.tsv",sep="")

if (7 %in% run.to) {
  
  print("Predict hotspot mutations for snv/indel")
  
  ## Step 7a ##
  if (file.exists(LRmodel.snv.file)) {
    
  results.snv = mutPredict.snv(mask.regions.file = mask.regions.file, nucleotide.selected.file = nucleotide.selected.file, continuous.features.selected.snv.url.file = continuous.features.selected.snv.url.file, discrete.features.selected.snv.url.file = discrete.features.selected.snv.url.file, sample.specific.features.url.file = sample.snv.features,
                         snv.mutations.file = snv.mutations.int, snv.mutations.file2 = snv.mutations, region.of.interest = region.of.interest, cores = cores, snv.model.file = LRmodel.snv.file, min.count = min.count.snv, genome.size = genome.size, hotspots = hotspots, merge.hotspots = merge.hotspots)

  if(nrow(results.snv[[1]]) != 0) {
    
write.table(results.snv[[1]], file = snv.hotspots, col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")

  }
  
  if (!is.null(results.snv[[2]])) {
    
    write.table(results.snv[[2]], file=snv.hotspots.merged, col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
    
  }
  
  }
  
## Step 7b ##
  if (file.exists(LRmodel.indel.file)) {
    
results.indel = mutPredict.indel(mask.regions.file = mask.regions.file, continuous.features.selected.indel.url.file = continuous.features.selected.indel.url.file, discrete.features.selected.indel.url.file = discrete.features.selected.indel.url.file, sample.specific.features.url.file = sample.indel.features,
                             indel.mutations.file = indel.mutations.int, indel.mutations.file2 = indel.mutations, indel.model.file = LRmodel.indel.file, region.of.interest = region.of.interest, cores = cores, min.count = min.count.indel, genome.size = genome.size, hotspots = hotspots, merge.hotspots = merge.hotspots)

if (nrow(results.indel[[1]]) != 0) {
  
write.table(results.indel[[1]], file = indel.hotspots, col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
  
}

if (!is.null(results.indel[[2]])) {
  
  write.table(results.indel[[2]],file=indel.hotspots.merged,col.names=TRUE,row.names=TRUE,quote=FALSE,sep="\t")
  
}

  }
  
}

## Step 8 ##
ann.snv.hotspots = paste(output.dir, "snv_hotspots_annotated.tsv", sep = "")
ann.indel.hotspots = paste(output.dir, "indel_hotspots_annotated.tsv", sep = "")

if (8 %in% run.to) {
  
  print("Annotate hotspots for snv/indel")
  
  ## Step 8a ##
  if (file.exists(snv.hotspots.merged)) {
    
    ann.results.snv = mutAnnotate(hotspots.file = snv.hotspots.merged, promoter.file = promoter.file, 
                                      utr3.file = utr3.file, utr5.file = utr5.file,
                                      other.annotations = other.annotations)
    
    write.table(ann.results.snv, file = ann.snv.hotspots, col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
    
  } else if (file.exists(snv.hotspots)) {
    
    ann.results.snv = mutAnnotate(hotspots.file = snv.hotspots, promoter.file = promoter.file, 
                                  utr3.file = utr3.file, utr5.file = utr5.file,
                                  other.annotations = other.annotations)
    
    write.table(ann.results.snv, file = ann.snv.hotspots, col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
    
    
  }
  
  ## Step 8b ##
  if (file.exists(indel.hotspots.merged)) {
    
    ann.results.indel = mutAnnotate(hotspots.file = indel.hotspots.merged, promoter.file = promoter.file,
                                         utr3.file = utr3.file, utr5.file = utr5.file,
                                    other.annotations = other.annotations)
    
    write.table(ann.results.indel, file = ann.indel.hotspots, col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
    
  } else if (file.exists(indel.hotspots)) {
    
    ann.results.indel = mutAnnotate(hotspots.file = indel.hotspots, promoter.file = promoter.file,
                                    utr3.file = utr3.file, utr5.file = utr5.file, 
                                    other.annotations = other.annotations)
    
    write.table(ann.results.indel, file = ann.indel.hotspots, col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
    
    
  }
  
}

## Step 9 ##
if (9.1 %in% run.to) {
  
  print("Plot feature importance")
  
  ## Step 9a ##
  if (file.exists(LRmodel.snv.file)) {
    
    plot_feature_importance(model = LRmodel.snv.file, mutCovariate.table.file = mutCovariate.snv.output.p1, mutCovariate.count.file = mutCovariate.snv.output.p2, continuous.features.selected.url.file = continuous.features.selected.snv.url.file, z.value = z.value, mutation.type = "SNV", output.dir = output.dir)
    
  }
  
  ## Step 9b ##
  if (file.exists(LRmodel.indel.file)) {
    
    plot_feature_importance(model = LRmodel.indel.file, mutCovariate.table.file = mutCovariate.indel.output.p1, mutCovariate.count.file = mutCovariate.indel.output.p2, continuous.features.selected.url.file = continuous.features.selected.indel.url.file, z.value = z.value, mutation.type = "indel", output.dir = output.dir)
    
  }
  
}

manhattan.snv = paste(output.dir, "snv_manhattan.pdf", sep = "")
manhattan.indel = paste(output.dir, "indel_manhattan.pdf", sep = "")

if (9.2 %in% run.to) {
  
  print("Plot manhattan figure for snv/indel")
  
  ## Step 9a ##
  if (file.exists(snv.hotspots.merged)) {
    
    pdf(manhattan.snv)
    plot_manhattan(hotspots.file = snv.hotspots.merged, fdr.cutoff = fdr.cutoff, color.line = color.line, color.dots = color.dots)
    dev.off()
    
  } else if (file.exists(snv.hotspots)) {
    
    pdf(manhattan.snv)
    plot_manhattan(hotspots.file = snv.hotspots, fdr.cutoff = fdr.cutoff, color.line = color.line, color.dots = color.dots)
    dev.off()
    
    
  }
  
  ## Step 9b ##
  if (file.exists(indel.hotspots.merged)) {
    
    pdf(manhattan.indel)
    plot_manhattan(hotspots.file = indel.hotspots.merged, fdr.cutoff = fdr.cutoff, color.line = color.line, color.dots = color.dots)
    dev.off()
    
  } else if (file.exists(indel.hotspots)) {
    
    pdf(manhattan.indel)
    plot_manhattan(hotspots.file = indel.hotspots, fdr.cutoff = fdr.cutoff, color.line = color.line, color.dots = color.dots)
    dev.off()
    
    
  }
  
}

if (9.3 %in% run.to) {
  
  print("Plot top hotspots for snv/indel")
  
  ## Step 9a
  if (file.exists(snv.hotspots.merged)) {
    
    plot_top_hits(hotspots.file = snv.hotspots.merged, fdr.cutoff = fdr.cutoff, color.muts = color.muts, mutations.file  = snv.mutations, mutation.type = "SNV", top.no = top.no, output.dir = output.dir)
    
  } else if (file.exists(snv.hotspots)) {
    
    plot_top_hits(hotspots.file = snv.hotspots, fdr.cutoff = fdr.cutoff, color.muts = color.muts, mutations.file  = snv.mutations, mutation.type = "SNV", top.no = top.no, output.dir = output.dir)
    
    
  }
  
  ## Step 9b
  if (file.exists(indel.hotspots.merged)) {
    
    plot_top_hits(hotspots.file = indel.hotspots.merged, fdr.cutoff = fdr.cutoff, color.muts = color.muts, mutations.file = indel.mutations, mutation.type = "indel", top.no = top.no, output.dir = output.dir)
    
  } else if (file.exists(indel.hotspots)) {
    
    plot_top_hits(hotspots.file = indel.hotspots, fdr.cutoff = fdr.cutoff, color.muts = color.muts, mutations.file = indel.mutations, mutation.type = "indel", top.no = top.no, output.dir = output.dir)
    
  }
  
}

}
