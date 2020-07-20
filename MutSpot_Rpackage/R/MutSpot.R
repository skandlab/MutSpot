#' Runs all or selected steps in MutSpot analysis.
#'
#' @param run.to Numeric vector defining which steps to run, default = 1, 2, 3.1, 3.2, 4.1, 4.2, 5.1, 5.2, 5.3, 5.4, 5.5, 6, 7.
#' @param working.dir Working directory, default = NULL will use current working directory.
#' @param chromosomes Character vector defining which chromosomes to compute feature matrix on, default = chr1-chrX.
#' @param snv.mutations SNV mutations MAF file, default = NULL.
#' @param indel.mutations Indel mutations MAF file, default = NULL.
#' @param mask.regions.file Regions to mask in genome, for example, non-mappable regions/immunoglobin loci/CDS regions RDS file, depends on genome build, default = mask_regions.RDS, Ch37.
#' @param all.sites.file All sites in whole genome RDS file, depends on genome build, default file = all_sites.RDS, Ch37.
#' @param region.of.interest Region of interest bed file, default = NULL.
#' @param sample To sample for non-mutated sites or to use all sites in region of interest, default = TRUE.
#' @param cores Number of cores, default = 1.
#' @param cutoff.nucleotide Frequency cutoff/threshold to determine nucleotide contexts used in prediction model, ranges from 0.5 to 1, default = 0.90.
#' @param cutoff.nucleotide.new Updated frequency cutoff/threshold to determine nucleotide contexts used in prediction model, ranges from 0.5 to 1.
#' @param top.nucleotide Number of top nucleotide contexts to select, default = NULL.
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
#' @param top.features Number of top genomic features to select, default = NULL.
#' @param fit.sparse To fit model using glmnet or glm, default = FALSE.
#' @param drop To drop insignificant features from fitted model or not, default = FALSE.
#' @param min.count.snv Minimum number of mutated samples in each SNV hotspot, default = 2.
#' @param min.count.indel Minimum number of mutated samples in each indel hotspot, default = 2.
#' @param hotspot.size Size of each hotspot, default = 21.
#' @param genome.size Genome size, depends on genome build, default = 2533374732, Ch37.
#' @param hotspots To run hotspot analysis or region-based analysis, default = TRUE.
#' @param collapse.regions To collapse region of interest or not, default = FALSE.
#' @param promoter.file Promoter regions bed file, depends on genome build, default file = Ensembl75.promoters.coding.bed, Ch37.
#' @param utr3.file 3'UTR regions bed file, depends on genome.build, default file = Ensembl75.3UTR.coding.bed, Ch37.
#' @param utr5.file 5'UTR regions bed file, depends on genome build, default file = Ensembl75.5UTR.coding.bed, Ch37.
#' @param other.annotations Text file containing URLs of additional regions to be annotated, default = NULL.
#' @param fdr.cutoff FDR cutoff, default = 0.05.
#' @param color.line Color given FDR cutoff, default = red.
#' @param color.dots Color hotspots that passed given FDR cutoff, default = maroon1.
#' @param merge.hotspots To plot overlapping hotspots as 1 hotspot or individual hotspots, default = TRUE.
#' @param color.muts Color points, default = orange.
#' @param top.no Number of top hotspots to plot, default = 3.
#' @param debug To keep intermediate output files or not, default = FALSE.
#' @param dilution.analysis To run dilution test or not, default = FALSE.
#' @param genome.build Reference genome build, default = Ch37.
#' @return Corresponding output from each step in MutSpot analysis.
#' @export

MutSpot = function(run.to = c(1:2, 3.1, 3.2, 4.1, 4.2, 5.1, 5.2, 5.3, 5.4, 5.5, 6, 7), working.dir = NULL, chromosomes = c(1:22, "X"), snv.mutations = NULL, indel.mutations = NULL, mask.regions.file = system.file("extdata", "mask_regions.RDS", package = "MutSpot"), all.sites.file = system.file("extdata", "all_sites.RDS", package = "MutSpot"), region.of.interest = NULL, sample = T, cores = 1, cutoff.nucleotide = 0.90, cutoff.nucleotide.new = NULL, top.nucleotide = NULL, genomic.features.snv = NULL, genomic.features.indel = NULL, genomic.features = NULL, genomic.features.fixed.snv = NULL, genomic.features.fixed.indel = NULL, genomic.features.fixed = NULL, sample.snv.features = NULL, sample.indel.features = NULL, cutoff.features = 0.75, cutoff.features.new.snv = NULL, cutoff.features.new.indel = NULL, top.features = NULL, fit.sparse = FALSE, drop = FALSE, min.count.snv = 2, min.count.indel = 2, hotspot.size = 21, genome.size = 2533374732, hotspots = TRUE, collapse.regions = FALSE,
                  promoter.file = system.file("extdata", "Ensembl75.promoters.coding.bed", package = "MutSpot"), utr3.file = system.file("extdata", "Ensembl75.3UTR.coding.bed", package = "MutSpot"), utr5.file = system.file("extdata", "Ensembl75.5UTR.coding.bed", package = "MutSpot"), other.annotations = NULL, fdr.cutoff = 0.05, color.line = "red", color.dots = "maroon1", merge.hotspots = TRUE, color.muts = "orange", top.no = 3, debug = FALSE, dilution.analysis = FALSE, genome.build = "Ch37") {
  
  # Set default files according to genome.build
  if (genome.build == "Ch38") {
    
  mask.regions.file = system.file("extdata", "mask_regions_hg38.RDS", package = "MutSpot")
  all.sites.file = system.file("extdata", "all_sites_hg38.RDS", package = "MutSpot")
  genome.size = 2797357708
  promoter.file = system.file("extdata", "promoter_hg38.bed", package = "MutSpot")
  utr3.file = system.file("extdata", "utr3_hg38.bed", package = "MutSpot")
  utr5.file = system.file("extdata", "utr5_hg38.bed", package = "MutSpot")
  
  }
  
  # Set working directory
  if (is.null(working.dir)) {
    
    working.dir = getwd()
    
  }
  
## check format of working directory ##
if (substr(working.dir, nchar(working.dir), nchar(working.dir)) != "/") {
  
  working.dir = paste(working.dir, "/", sep = "")
  
}  
  setwd(working.dir)
  
  output.dir = paste(working.dir, "results/", sep = "")
  if (!dir.exists(output.dir)) {
    
  dir.create(output.dir)
    
  }
  
  features.dir = paste(working.dir, "features/", sep = "")
  if (!dir.exists(features.dir)) {
    
    dir.create(features.dir)
    
  }
  
  ## Load all dependencies ##
  if("data.table" %in% rownames(installed.packages()) == FALSE) {
    
    print("install [data.table]")
    install.packages("data.table")
    
  }
  if("pscl" %in% rownames(installed.packages()) == FALSE) {
    
    print("install [pscl]")
    install.packages("pscl")
    
  }
  if (grepl("R version 3.6", R.Version()$version.string)) {
    
    if("BiocManager" %in% rownames(installed.packages()) == FALSE) {
      
      print("install [BiocManager]")
      install.packages("BiocManager")
      
    }
    if ("BSgenome.Hsapiens.UCSC.hg19" %in% rownames(installed.packages()) == FALSE & genome.build == "Ch37") {
      
      print("install [BSgenome.Hsapiens.UCSC.hg19]")
      BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")

    }
    if ("BSgenome.Hsapiens.UCSC.hg38" %in% rownames(installed.packages()) == FALSE & genome.build == "Ch38") {
      
      print("install [BSgenome.Hsapiens.UCSC.hg38]")
      BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
      
    }
    if("Biostrings" %in% rownames(installed.packages()) == FALSE) {
      
      print("install [Biostrings]")
      BiocManager::install("Biostrings")

    }
    if("GenomicRanges" %in% rownames(installed.packages()) == FALSE) {
      
      print("install [GenomicRanges]")
      BiocManager::install("GenomicRanges")

    }
    
      } else if(grepl("R version 3.2", R.Version()$version.string)) {
        
        if ("BSgenome.Hsapiens.UCSC.hg19" %in% rownames(installed.packages()) == FALSE & genome.build == "Ch37") {
          
          print("install [BSgenome.Hsapiens.UCSC.hg19]")
          source("http://bioconductor.org/biocLite.R")
          biocLite("BSgenome.Hsapiens.UCSC.hg19")
          
        }
        if ("BSgenome.Hsapiens.UCSC.hg38" %in% rownames(installed.packages()) == FALSE & genome.build == "Ch38") {
          
          print("install [BSgenome.Hsapiens.UCSC.hg38]")
          source("http://bioconductor.org/biocLite.R")
          biocLite("BSgenome.Hsapiens.UCSC.hg38")
          
        }
        if("Biostrings" %in% rownames(installed.packages()) == FALSE) {
          
          print("install [Biostrings]")
          source("http://bioconductor.org/biocLite.R")
          biocLite("Biostrings")
          
        }
        if("GenomicRanges" %in% rownames(installed.packages()) == FALSE) {
          
          print("install [GenomicRanges]")
          source("http://bioconductor.org/biocLite.R")
          biocLite("GenomicRanges")
          
        }
        
      } else {
    
  if ("BSgenome.Hsapiens.UCSC.hg19" %in% rownames(installed.packages()) == FALSE & genome.build == "Ch37") {
    
    print("install [BSgenome.Hsapiens.UCSC.hg19]")
    source("https://bioconductor.org/biocLite.R")
    biocLite("BSgenome.Hsapiens.UCSC.hg19")
    
  }
  if ("BSgenome.Hsapiens.UCSC.hg38" %in% rownames(installed.packages()) == FALSE & genome.build == "Ch38") {
          
          print("install [BSgenome.Hsapiens.UCSC.hg38]")
          source("https://bioconductor.org/biocLite.R")
          biocLite("BSgenome.Hsapiens.UCSC.hg38")
          
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
  if (genome.build == "Ch37") {
    
  suppressWarnings(suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19)))
    
  } else if (genome.build == "Ch38") {
    
    suppressWarnings(suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38)))
    
  }
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
downsampled.sites.snv.file = paste(output.dir, "downsampled.sites.snv.RDS", sep = "")
downsampled.sites.indel.file = paste(output.dir, "downsampled.sites.indel.RDS", sep = "")
  
if (1 %in% run.to) {
    
  print("Sample sites from SNV/indel")
  
sample_sites = sample.sites(snv.mutations.file = snv.mutations, indel.mutations.file = indel.mutations, mask.regions.file = mask.regions.file, all.sites.file = all.sites.file, region.of.interest = region.of.interest, sample = sample, cores = cores, genome.build = genome.build)

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

# Save snv downsampled mutated sites as RDS file
if (!is.null(sample_sites[[5]])) {
  
  saveRDS(sample_sites[[5]], downsampled.sites.snv.file)
  
}

# Save indel downsampled mutated sites as RDS file
if (!is.null(sample_sites[[6]])) {
  
  saveRDS(sample_sites[[6]], downsampled.sites.indel.file)
  
}

rm(sample_sites)
gc()

}

if (!is.null(region.of.interest) & !is.null(snv.mutations)) {
  
  if (file.exists(snv.mutations.region.file)) {
    
  snv.mutations.int = snv.mutations.region.file
  
  } else {
    
    snv.mutations.int = NULL
    
  }
  
} else {
  
  snv.mutations.int = snv.mutations
  
}

if (!is.null(region.of.interest) & !is.null(indel.mutations)) {
  
  if (file.exists(indel.mutations.region.file)) {
  
  indel.mutations.int = indel.mutations.region.file
  
  } else {
    
    indel.mutations.int = NULL
    
  }
  
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
local.mutrate.snv.file = paste(features.dir, "localmutrate_snv.bed", sep = "")
local.mutrate.indel.file = paste(features.dir, "localmutrate_indel.bed", sep = "")

  if (2 %in% run.to) {
    
    print("Calculate local mutation rates for SNV/indel")
    
local_mutrate = local.mutrate(snv.mutations.file = snv.mutations, indel.mutations.file = indel.mutations, genome.build = genome.build)

# Save binned SNV local mutation rate as bed file
if (!is.null(local_mutrate[[1]])) {
  
  write.table(local_mutrate[[1]], file = local.mutrate.snv.file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
}

# Save binned indel local mutation rate as bed file
if (!is.null(local_mutrate[[2]])) {
  
  write.table(local_mutrate[[2]], file = local.mutrate.indel.file, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  
}

rm(local_mutrate)
gc()

}

## Step 3 ##
nucleotide.stabs.freq.file = paste(output.dir, "nucleotide_stabs_freq.RDS", sep = "")
nucleotide.selected.file = paste(output.dir, "nucleotide_selected.RDS", sep = "")
indels.polyAT.file = paste(output.dir, "indel_polyAT.MAF", sep = "")

## Step 3I ##
if (3.1 %in% run.to) {
  
  print("Select nucleotide context for SNV and filter polyAT from indel")
    
nucleotide_selection = nucleotide.selection(sampled.sites.snv.file = sampled.sites.snv.file, indel.mutations.file = indel.mutations.int, cutoff = cutoff.nucleotide, cores = cores, genome.build = genome.build)

# # Check for need to change threshold
# print(head(nucleotide_selection[[1]][order(nucleotide_selection[[1]]$f, decreasing=TRUE), ]))

# Save lasso selection frequency table as RDS file
if (!is.null(nucleotide_selection[[1]])) {
  
  saveRDS(nucleotide_selection[1:2], file = nucleotide.stabs.freq.file)
  
}

# Save selected nucleotide contexts as RDS file
if (!is.null(nucleotide_selection[[3]]) & length(nucleotide_selection[[3]] != 0)) {
  
  saveRDS(nucleotide_selection[[3]], file = nucleotide.selected.file)
  
}

# Save filtered indel mutations as MAF file
if (!is.null(nucleotide_selection[[4]])) {
  
  if (nrow(nucleotide_selection[[4]]) != 0) {
  
  write.table(nucleotide_selection[[4]], file = indels.polyAT.file, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
  # indel.mutations = indels.polyAT.file
  
  }
  
}

rm(nucleotide_selection)
gc()

}

## Step 3II ##
if (3.2 %in% run.to) {
  
  print("Reselect nucleotide context for SNV using altered threshold")
  
if (!is.null(cutoff.nucleotide.new)) {
  
nucleotide_selection_adjust = nucleotide.selection.adjust(nucleotide.stabs.file = nucleotide.stabs.freq.file, new.cutoff = cutoff.nucleotide.new, top.nucleotide = top.nucleotide)

if (length(nucleotide_selection_adjust) != 0) {
  
saveRDS(nucleotide_selection_adjust, file = nucleotide.selected.file)
  
} else {
  
  file.remove(nucleotide.selected.file)
  
}

rm(nucleotide_selection_adjust)
gc()

}

}

## Step 4 ##
features.stabs.snv.file = paste(output.dir, "features_stabs_snv.RDS", sep = "")
continuous.features.selected.snv.url.file = paste(output.dir, "continuous_features_selected_snv_url.txt", sep = "")
discrete.features.selected.snv.url.file = paste(output.dir, "discrete_features_selected_snv_url.txt", sep = "")
features.stabs.indel.file = paste(output.dir, "features_stabs_indel.RDS", sep = "")
continuous.features.selected.indel.url.file = paste(output.dir, "continuous_features_selected_indel_url.txt", sep = "")
discrete.features.selected.indel.url.file = paste(output.dir, "discrete_features_selected_indel_url.txt", sep = "")
features.sds.file = paste(output.dir, "features_sds.RDS", sep = "")

if (4.1 %in% run.to) {
  
  print("Select epigenetic features for SNV/indel")
  
    ## Step 4I ##
epigenetic_selection = epigenetic.selection(sampled.sites.snv.file = sampled.sites.snv.file, sampled.sites.indel.file = sampled.sites.indel.file, genomic.features.snv = genomic.features.snv, genomic.features.indel = genomic.features.indel, genomic.features = genomic.features, genomic.features.fixed.snv = genomic.features.fixed.snv, genomic.features.fixed.indel = genomic.features.fixed.indel, genomic.features.fixed = genomic.features.fixed, cores = cores, cutoff = cutoff.features, feature.dir = features.dir, genome.build = genome.build)

# Check for need to change threshold for SNVs
# if (!is.null(epigenetic_selection[[1]])) {
  
# print(head(epigenetic_selection[[1]][order(epigenetic_selection[[1]]$f, decreasing = TRUE), ]))

# }

# Save lasso selection frequency table as RDS file 
if (!is.null(epigenetic_selection[[1]])) {
  
  saveRDS(epigenetic_selection[1:2], file = features.stabs.snv.file)
  
}

# Save URLs of selected SNV continuous features as text file
if (!is.null(epigenetic_selection[[3]])) {
  
  if (nrow(epigenetic_selection[[3]]) != 0) {
  
  write.table(epigenetic_selection[[3]], file = continuous.features.selected.snv.url.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  }
  
}

# Save URLs of selected SNV discrete features as text file
if (!is.null(epigenetic_selection[[4]])) {
  
  if (nrow(epigenetic_selection[[4]]) != 0) {
  
  write.table(epigenetic_selection[[4]], file = discrete.features.selected.snv.url.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  }

}

# Check for need to change threshold for indels
# if (!is.null(epigenetic_selection[[5]])) {
#   
# print(head(epigenetic_selection[[5]][order(epigenetic_selection[[4]]$f, decreasing = TRUE), ]))
# 
# }

# Save lasso selection frequency table as RDS file
if (!is.null(epigenetic_selection[[5]])) {
  
  saveRDS(epigenetic_selection[5:6], file = features.stabs.indel.file)
  
}

# Save URLs of selected indel continuous features as text file
if (!is.null(epigenetic_selection[[7]])) {
  
  if (nrow(epigenetic_selection[[7]]) != 0) {
  
  write.table(epigenetic_selection[[7]], file = continuous.features.selected.indel.url.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  }
  
}

# Save URLs of selected indel discrete features as text file
if (!is.null(epigenetic_selection[[8]])) {
  
  if (nrow(epigenetic_selection[[8]]) != 0) {
  
  write.table(epigenetic_selection[[8]], file = discrete.features.selected.indel.url.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  }
  
}

# Save SDs of features as RDS file
# if (!is.null(epigenetic_selection[[9]] | !is.null(epigenetic_selection[[10]]))) {
#   
  saveRDS(epigenetic_selection[9:10], file = features.sds.file)
  
# }

  rm(epigenetic_selection)
  gc()
  
}

if (!file.exists(features.stabs.snv.file)) {
  
  features.stabs.snv.file = NULL
  
}
if (!file.exists(features.stabs.indel.file)) {
  
  features.stabs.indel.file = NULL
  
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

## Step 4II ##
if (4.2 %in% run.to) {
  
  print("Reselect epigenetic features for SNV/indel using altered threshold")
  
  continuous.features.selected.snv.url.file = paste(output.dir, "continuous_features_selected_snv_url.txt", sep = "")
  discrete.features.selected.snv.url.file = paste(output.dir, "discrete_features_selected_snv_url.txt", sep = "")
  continuous.features.selected.indel.url.file = paste(output.dir, "continuous_features_selected_indel_url.txt", sep = "")
  discrete.features.selected.indel.url.file = paste(output.dir, "discrete_features_selected_indel_url.txt", sep = "")

if (!is.null(cutoff.features.new.snv) | !is.null(cutoff.features.new.indel)) {

  epigenetic_selection_adjust = epigenetic.selection.adjust(feature.stabs.snv.file = features.stabs.snv.file , feature.stabs.indel.file = features.stabs.indel.file, continuous.features.selected.snv.url.file = continuous.features.selected.snv.url.file, discrete.features.selected.snv.url.file = discrete.features.selected.snv.url.file, continuous.features.selected.indel.url.file = continuous.features.selected.indel.url.file, discrete.features.selected.indel.url.file = discrete.features.selected.indel.url.file, new.cutoff.snv = cutoff.features.new.snv, new.cutoff.indel = cutoff.features.new.indel, top.features = top.features, features.sds = features.sds.file, genomic.features.snv = genomic.features.snv, genomic.features.indel = genomic.features.indel, genomic.features = genomic.features, feature.dir = features.dir)
  
  # Save URLs of selected SNV continuous features as text file
  if (!is.null(epigenetic_selection_adjust[[1]])) {
    
    if (nrow(epigenetic_selection_adjust[[1]] != 0)) {
      
    write.table(epigenetic_selection_adjust[[1]], file = continuous.features.selected.snv.url.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    } else {
      
      file.remove(continuous.features.selected.snv.url.file)
      
    }
    
  }
  
  # Save URLs of selected SNV discrete features as text file
  if (!is.null(epigenetic_selection_adjust[[2]])) {
    
    if (nrow(epigenetic_selection_adjust[[2]] != 0)) {
    
    write.table(epigenetic_selection_adjust[[2]], file = discrete.features.selected.snv.url.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    } else {
      
      file.remove(discrete.features.selected.snv.url.file)
      
    }
    
  }
  
  # Save URLs of selected indel continuous features as text file
  if (!is.null(epigenetic_selection_adjust[[3]])) {
    
    if (nrow(epigenetic_selection_adjust[[3]] != 0)) {
      
    write.table(epigenetic_selection_adjust[[3]], file = continuous.features.selected.indel.url.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    } else {
      
      file.remove(continuous.features.selected.indel.url.file)
      
    }
    
  }
  
  # Save URLs of selected indel discrete features as text file
  if (!is.null(epigenetic_selection_adjust[[4]])) {
    
    if (nrow(epigenetic_selection_adjust[[4]] != 0)) {
      
    write.table(epigenetic_selection_adjust[[4]], file = discrete.features.selected.indel.url.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    } else {
      
      file.remove(discrete.features.selected.indel.url.file)
      
    }
    
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
  
  rm(epigenetic_selection_adjust)
  gc()
  
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
if (!is.null(continuous.features.selected.snv.url.file)) {
  
if (!file.exists(continuous.features.selected.snv.url.file)) {
  
  continuous.features.selected.snv.url.file = NULL
  
}
  
}
if (!is.null(discrete.features.selected.snv.url.file)) {
  
if (!file.exists(discrete.features.selected.snv.url.file)) {
  
  discrete.features.selected.snv.url.file = NULL
  
}
  
}
if (!is.null(continuous.features.selected.indel.url.file)) {
  
if (!file.exists(continuous.features.selected.indel.url.file)) {
  
  continuous.features.selected.indel.url.file = NULL
  
}
  
}
if (!is.null(discrete.features.selected.indel.url.file)) {
  
if (!file.exists(discrete.features.selected.indel.url.file)) {
  
  discrete.features.selected.indel.url.file = NULL
  
}
  
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
    
    # Define all sites in whole genome
    all.sites = readRDS(all.sites.file)
    all.sites = all.sites[as.character(GenomeInfoDb::seqnames(all.sites)) %in% seqnames]
    all.sites.masked = subtract.regions.from.roi(all.sites, mask.regions, cores = cores)
    # sum(as.numeric(GenomicRanges::width(all.sites.masked)))
    
    # If specified region, redefine all sites to be in specified region
    if (!is.null(region.of.interest)) {
      
      print("specified region")
      
      all.sites = bed.to.granges(region.of.interest)
      all.sites.masked = subtract.regions.from.roi(all.sites, mask.regions, cores = cores)
      all.sites.masked = all.sites.masked[as.character(GenomeInfoDb::seqnames(all.sites.masked)) %in% seqnames]
      
    }
    
    # Read feature file paths
    if (!is.null(continuous.features.selected.snv.url.file)) {
      
      selected.continuous.urls <- read.delim(continuous.features.selected.snv.url.file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
      continuous.selected.features = parallel::mclapply(selected.continuous.urls[ ,2], function(f) {
        
        print(f)
        df = bed.to.granges(as.character(f))
        
      }, mc.cores = cores)
      
      names(continuous.selected.features) = as.character(selected.continuous.urls[ ,1])
      
    } else {
      
      continuous.selected.features = NULL
      
    }
    
    if (!is.null(discrete.features.selected.snv.url.file)) {
      
      selected.discrete.urls <- read.delim(discrete.features.selected.snv.url.file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
      discrete.selected.features = parallel::mclapply(selected.discrete.urls[ ,2], function(f) {
        
        print(f)
        df = bed.to.granges(as.character(f))
        
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
      
      sample.specific.features = read.delim(sample.specific.features.url.file, stringsAsFactors = FALSE)
      rownames(sample.specific.features) = as.character(sample.specific.features$SampleID)
      sample.specific.features = sample.specific.features[which(sample.specific.features$SampleID %in% names(ind.mut.count)), ]
      sample.specific.features$ind.mut.count = ind.mut.count[rownames(sample.specific.features)]
      sample.specific.features = sample.specific.features[ ,-which(colnames(sample.specific.features) == "SampleID")]
      
      sample.specific.features2 = parallel::mclapply(1:ncol(sample.specific.features), FUN = function(x) {
        
        print(colnames(sample.specific.features)[x])
        if(class(sample.specific.features[ ,x]) == "character") {
          
          t = factor(sample.specific.features[ ,x])
          t = model.matrix( ~ t)[ ,-1]
          if (class(t) == "matrix") {
            
          colnames(t) = substr(colnames(t), 2, nchar(colnames(t)))
          colnames(t) = paste(colnames(sample.specific.features)[x], colnames(t), sep = "")
          rownames(t) = rownames(sample.specific.features)
          
          } else {
            
            t = as.data.frame(t)
            colnames(t) = paste(colnames(sample.specific.features)[x], levels(factor(sample.specific.features[,x]))[2], sep = "")
            rownames(t) = rownames(sample.specific.features)
            
          }
          
        } else {
          
          t = as.data.frame(sample.specific.features[ ,x])
          colnames(t) = colnames(sample.specific.features)[x]
          rownames(t) = rownames(sample.specific.features)
          
        }
        return(t)
        
      },mc.cores=cores)
      
      sample.specific.features = do.call(cbind, sample.specific.features2)
      
    } else {
      
      sample.specific.features = as.data.frame(ind.mut.count)
      colnames(sample.specific.features) = "ind.mut.count"
      
    }
    
    # Remove larger objects before tabulating
    sort(sapply(ls(), function(x) { object.size(get(x)) / 10 ^ 6 }))
    rm(list = c("maf.mutations", "maf.ind", "mask.regions", "all.sites", "maf.mutations2"))
    gc(reset = T)
    
    # Tabulate covariates for mutations
    GenomeInfoDb::seqlevels(mut.masked) = as.character(unique(GenomeInfoDb::seqnames(mut.masked)))
    mut.chr = GenomicRanges::split(mut.masked, GenomeInfoDb::seqnames(mut.masked))
    if (genome.build == "Ch37") {
      
    chrs <- names(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[1:23]
    
    } else if (genome.build == "Ch38") {

      chrs <- names(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[1:23]
      
    }
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

        # Extract all nucleotide contexts' positions in specific chromosome
        for (i in 1:nrow(nucleotide.selected)) {
          
          print(paste(chr.interest, nucleotide.selected[i, "sequence"], nucleotide.selected[i, "type"], sep = ":"))
          if (nucleotide.selected$type[i] %in% c("oneMer", "threeMer", "fiveMer")) {
            
            if (genome.build == "Ch37") {
              
            precompute.motif.pos[[paste(nucleotide.selected[i, "type"], nucleotide.selected[i, "sequence"], sep = "")]] <- unique(c(BSgenome::start(Biostrings::matchPattern(nucleotide.selected[i, "sequence"], BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[chr.interest]])), BSgenome::start(Biostrings::matchPattern(as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(nucleotide.selected[i, "sequence"]))), BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[chr.interest]]))))
            
            } else if (genome.build == "Ch38") {
              
              precompute.motif.pos[[paste(nucleotide.selected[i, "type"], nucleotide.selected[i, "sequence"], sep = "")]] <- unique(c(BSgenome::start(Biostrings::matchPattern(nucleotide.selected[i, "sequence"], BSgenome.Hsapiens.UCSC.hg38::Hsapiens[[chr.interest]])), BSgenome::start(Biostrings::matchPattern(as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(nucleotide.selected[i, "sequence"]))), BSgenome.Hsapiens.UCSC.hg38::Hsapiens[[chr.interest]]))))
              
            }
            
          } else if (nucleotide.selected$type[i] == "threeRight") {
            
            if (genome.build == "Ch37") {
              
            a <- unique(BSgenome::start(Biostrings::matchPattern(nucleotide.selected[i, "sequence"], BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[chr.interest]])))
            b <- unique(BSgenome::start(Biostrings::matchPattern(as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(nucleotide.selected[i, "sequence"]))), BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[chr.interest]]))) + 1
            
            } else if (genome.build == "Ch38") {
              
              a <- unique(BSgenome::start(Biostrings::matchPattern(nucleotide.selected[i, "sequence"], BSgenome.Hsapiens.UCSC.hg38::Hsapiens[[chr.interest]])))
              b <- unique(BSgenome::start(Biostrings::matchPattern(as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(nucleotide.selected[i, "sequence"]))), BSgenome.Hsapiens.UCSC.hg38::Hsapiens[[chr.interest]]))) + 1
              
            }
            precompute.motif.pos[[paste(nucleotide.selected[i, "type"], nucleotide.selected[i, "sequence"], sep = "")]] = unique(c(a, b))
            
          } else if (nucleotide.selected$type[i] == "threeLeft") {
            
            if (genome.build == "Ch37") {
              
            a <- unique(BSgenome::start(Biostrings::matchPattern(nucleotide.selected[i, "sequence"], BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[chr.interest]]))) + 1
            b <- unique(BSgenome::start(Biostrings::matchPattern(as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(nucleotide.selected[i, "sequence"]))), BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[chr.interest]])))
            
            } else if (genome.build == "Ch38") {
              
              a <- unique(BSgenome::start(Biostrings::matchPattern(nucleotide.selected[i, "sequence"], BSgenome.Hsapiens.UCSC.hg38::Hsapiens[[chr.interest]]))) + 1
              b <- unique(BSgenome::start(Biostrings::matchPattern(as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(nucleotide.selected[i, "sequence"]))), BSgenome.Hsapiens.UCSC.hg38::Hsapiens[[chr.interest]])))
              
            }
            precompute.motif.pos[[paste(nucleotide.selected[i, "type"], nucleotide.selected[i, "sequence"], sep = "")]] = unique(c(a, b))
            
          }  else if (nucleotide.selected$type[i] == "fiveRight") {
            
            if (genome.build == "Ch37") {
              
            a <- unique(BSgenome::start(Biostrings::matchPattern(nucleotide.selected[i, "sequence"], BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[chr.interest]])))
            b <- unique(BSgenome::start(Biostrings::matchPattern(as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(nucleotide.selected[i, "sequence"]))), BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[chr.interest]]))) + 2
            
            } else if (genome.build == "Ch38") {
              
              a <- unique(BSgenome::start(Biostrings::matchPattern(nucleotide.selected[i, "sequence"], BSgenome.Hsapiens.UCSC.hg38::Hsapiens[[chr.interest]])))
              b <- unique(BSgenome::start(Biostrings::matchPattern(as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(nucleotide.selected[i, "sequence"]))), BSgenome.Hsapiens.UCSC.hg38::Hsapiens[[chr.interest]]))) + 2
              
            }
            precompute.motif.pos[[paste(nucleotide.selected[i, "type"], nucleotide.selected[i, "sequence"], sep = "")]] = unique(c(a, b))
            
          } else if (nucleotide.selected$type[i] == "fiveLeft") {
            
            if (genome.build == "Ch37") {
              
            a <- unique(BSgenome::start(Biostrings::matchPattern(nucleotide.selected[i, "sequence"], BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[chr.interest]]))) + 2
            b <- unique(BSgenome::start(Biostrings::matchPattern(as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(nucleotide.selected[i, "sequence"]))), BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[chr.interest]])))
            
            } else if (genome.build == "Ch38") {
              
              a <- unique(BSgenome::start(Biostrings::matchPattern(nucleotide.selected[i, "sequence"], BSgenome.Hsapiens.UCSC.hg38::Hsapiens[[chr.interest]]))) + 2
              b <- unique(BSgenome::start(Biostrings::matchPattern(as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(nucleotide.selected[i, "sequence"]))), BSgenome.Hsapiens.UCSC.hg38::Hsapiens[[chr.interest]])))
              
            }
            precompute.motif.pos[[paste(nucleotide.selected[i, "type"], nucleotide.selected[i, "sequence"], sep = "")]] = unique(c(a, b))
            
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
      
      rm(list=c("mut.freq","genome.freq.aggregated","precompute.motif.pos"))
      gc()
      
    }
    
  }
  
  ## Step 5b ##
  if ((!is.null(continuous.features.selected.indel.url.file) | !is.null(discrete.features.selected.indel.url.file) | !is.null(sample.indel.features)) & !is.null(sampled.sites.indel.file)) {
    
    sample.specific.features.url.file = sample.indel.features
    indel.mutations.file = indel.mutations.int
    indel.mutations.file2 = indel.mutations
    
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
    
    # Define all sites in whole genome
    all.sites = readRDS(all.sites.file)
    all.sites = all.sites[as.character(GenomeInfoDb::seqnames(all.sites)) %in% seqnames]
    all.sites.masked = subtract.regions.from.roi(all.sites, mask.regions, cores = cores)
    # sum(as.numeric(GenomicRanges::width(all.sites.masked)))
    
    # If specified region, redefine all sites to be in specified region
    if (!is.null(region.of.interest)) {
      
      print("specified region")
      
      all.sites = bed.to.granges(region.of.interest)
      all.sites.masked = subtract.regions.from.roi(all.sites, mask.regions, cores = cores)
      all.sites.masked = all.sites.masked[as.character(GenomeInfoDb::seqnames(all.sites.masked)) %in% seqnames]
      
    }
    
    # Extract all polyA, poly C, poly G and poly T positions in a specific chromosome
    if (genome.build == "Ch37") {
      
    chrs <- names(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[1:23]
    
    } else if (genome.build == "Ch38") {
      
      chrs <- names(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)[1:23]
      
    }
    
    # Read feature file paths
    if (!is.null(continuous.features.selected.indel.url.file)) {
      
      selected.continuous.urls <- read.delim(continuous.features.selected.indel.url.file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
      continuous.selected.features = parallel::mclapply(selected.continuous.urls[ ,2], function(f) {
        
        print(f)
        df = bed.to.granges(as.character(f)) 
        
        }, mc.cores = cores)
      names(continuous.selected.features) = as.character(selected.continuous.urls[ ,1])
      
    } else {
      
      continuous.selected.features = NULL
      
    }
    
    if (!is.null(discrete.features.selected.indel.url.file)) {
      
      selected.discrete.urls <- read.delim(discrete.features.selected.indel.url.file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
      discrete.selected.features = parallel::mclapply(selected.discrete.urls[ ,2], function(f) {
        
        print(f)
        df = bed.to.granges(as.character(f))
        
        }, mc.cores = cores)
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
      
      sample.specific.features = read.delim(sample.specific.features.url.file, stringsAsFactors = FALSE)
      rownames(sample.specific.features) = as.character(sample.specific.features$SampleID)
      sample.specific.features = sample.specific.features[which(sample.specific.features$SampleID %in% names(ind.mut.count)), ]
      sample.specific.features$ind.mut.count = ind.mut.count[rownames(sample.specific.features)]
      sample.specific.features = sample.specific.features[ ,-which(colnames(sample.specific.features) == "SampleID")]
      
      sample.specific.features2 = parallel::mclapply(1:ncol(sample.specific.features), FUN = function(x) {
        
        print(colnames(sample.specific.features)[x])
        if(class(sample.specific.features[ ,x]) == "character") {
          
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
      
    } else {
      
      sample.specific.features = as.data.frame(ind.mut.count)
      colnames(sample.specific.features) = "ind.mut.count"
      
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
      
      if (genome.build == "Ch37") {
        
      polyA <- lapply(chr.interest, function(x) BSgenome::start(Biostrings::matchPattern("AAAAA",
                                                                                         BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[x]])))
      
      } else if (genome.build == "Ch38") {
        
        polyA <- lapply(chr.interest, function(x) BSgenome::start(Biostrings::matchPattern("AAAAA",
                                                                                           BSgenome.Hsapiens.UCSC.hg38::Hsapiens[[x]])))
        
      }
      names(polyA) <- chr.interest
      polyAs = parallel::mclapply(chr.interest, FUN = function(x) {
        
        print(paste(x, "polyA", sep = ":"))
        gr = GenomicRanges::GRanges(x, IRanges::IRanges(polyA[[x]], polyA[[x]] + 5 - 1))
        gr
        
      }, mc.cores = cores)
      polyAs = BiocGenerics::unlist(GenomicRanges::GRangesList(polyAs))
      
      if (genome.build == "Ch37") {
        
      polyC <- lapply(chr.interest, function(x) BSgenome::start(Biostrings::matchPattern("CCCCC",
                                                                                         BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[x]])))
      
      } else if (genome.build == "Ch38") {
        
        polyC <- lapply(chr.interest, function(x) BSgenome::start(Biostrings::matchPattern("CCCCC",
                                                                                           BSgenome.Hsapiens.UCSC.hg38::Hsapiens[[x]])))
        
      }
      names(polyC) <- chr.interest
      polyCs = parallel::mclapply(chr.interest, FUN = function(x) {
        
        print(paste(x, "polyC", sep = ":"))
        gr = GenomicRanges::GRanges(x, IRanges::IRanges(polyC[[x]], polyC[[x]] + 5 - 1))
        gr
        
      }, mc.cores = cores)
      polyCs = BiocGenerics::unlist(GenomicRanges::GRangesList(polyCs))
      
      if (genome.build == "Ch37") {
        
      polyG <- lapply(chr.interest, function(x) BSgenome::start(Biostrings::matchPattern("GGGGG",
                                                                                         BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[x]])))
      
      } else if (genome.build == "Ch38") {
        
        polyG <- lapply(chr.interest, function(x) BSgenome::start(Biostrings::matchPattern("GGGGG",
                                                                                           BSgenome.Hsapiens.UCSC.hg38::Hsapiens[[x]])))
        
      }
      names(polyG) <- chr.interest
      polyGs = parallel::mclapply(chr.interest, FUN = function(x) {
        
        print(paste(x, "polyG", sep = ":"))
        gr = GenomicRanges::GRanges(x, IRanges::IRanges(polyG[[x]], polyG[[x]] + 5 - 1))
        gr
        
      }, mc.cores = cores)
      polyGs = BiocGenerics::unlist(GenomicRanges::GRangesList(polyGs))
      
      if (genome.build == "Ch37") {
        
      polyT <- lapply(chr.interest, function(x) BSgenome::start(Biostrings::matchPattern("TTTTT",
                                                                                         BSgenome.Hsapiens.UCSC.hg19::Hsapiens[[x]])))
      
      } else if (genome.build == "Ch38") {
        
        polyT <- lapply(chr.interest, function(x) BSgenome::start(Biostrings::matchPattern("TTTTT",
                                                                                           BSgenome.Hsapiens.UCSC.hg38::Hsapiens[[x]])))
        
      }
      names(polyT) <- chr.interest
      polyTs = parallel::mclapply(chr.interest, FUN = function(x) {
        
        print(paste(x, "polyT", sep = ":"))
        gr = GenomicRanges::GRanges(x, IRanges::IRanges(polyT[[x]], polyT[[x]] + 5 - 1))
        gr
        
      }, mc.cores = cores)
      polyTs = BiocGenerics::unlist(GenomicRanges::GRangesList(polyTs))
      
      polyAT = c(polyAs, polyTs)
      polyCG = c(polyCs, polyGs)
      
      if (chr.interest %in% names(mut.chr)) {
        
        mut.freq <- mutCovariate.indel.freq.table.muts(continuous.features = continuous.selected.features, discrete.features = discrete.selected.features, sample.specific.features = sample.specific.features, polyAT = polyAT, polyCG = polyCG, sites = mut.chr[[chr.interest]])
        
      } else {
        
        mut.freq = NULL
        
      }
      
      if (chr.interest %in% names(mut.indel.chr)) {
        
        indel.freq <- mutCovariate.indel.freq.table.muts(continuous.features = continuous.selected.features, discrete.features = discrete.selected.features, sample.specific.features = sample.specific.features, polyAT = polyAT, polyCG = polyCG, sites = mut.indel.chr[[chr.interest]])

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
        
        genome.freq <- parallel::mclapply(chunks, function(x) mutCovariate.indel.freq.table.genome(continuous.features = continuous.selected.features, discrete.features = discrete.selected.features, polyAT = polyAT, polyCG = polyCG, sites = all.sites.masked2[x]), mc.cores = cores, mc.preschedule = FALSE, mc.silent = FALSE)
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
      
      rm(list = c("mut.freq", "indel.freq", "genome.freq.aggregated", "polyAs", "polyTs", "polyGs", "polyCs", "polyAT", "polyCG"))
      gc()
      
    }
    
  }
  
}
  
if (5.2 %in% run.to) {
  
  print("Compile feature matrix for snv for all chromosomes")
  
  if ((!is.null(continuous.features.selected.snv.url.file) | !is.null(discrete.features.selected.snv.url.file) | !is.null(sample.snv.features) | !is.null(nucleotide.selected.file)) & !is.null(sampled.sites.snv.file)) {
  
mutCovariate_snv_compile = mutCovariate.snv.compile(mask.regions.file = mask.regions.file, all.sites.file = all.sites.file, snv.mutations.file = snv.mutations.int, sample.specific.features.url.file = sample.snv.features, region.of.interest = region.of.interest, cores = cores, snv.mutations.file2 = snv.mutations, chrom.dir = output.dir, genome.build = genome.build)

saveRDS(mutCovariate_snv_compile, file = mutCovariate.snv.output.file)

if (!debug) {
  
delete.files = Sys.glob(paste(output.dir, "mutCovariate_chr*.RDS", sep = ""))
for (i in delete.files) {
  
  unlink(i)
  
}

}

rm(mutCovariate_snv_compile)
gc()

  }
  
}

if (5.3 %in% run.to) {
  
  print("Convert feature matrix for snv to sparse matrix")
  
  if (file.exists(mutCovariate.snv.output.file)) {
    
mutCovariate_snv_sparse = mutCovariate.snv.sparse(compiled = mutCovariate.snv.output.file)

saveRDS(mutCovariate_snv_sparse[[1]], file = mutCovariate.snv.output.p1)
saveRDS(mutCovariate_snv_sparse[[2]], file = mutCovariate.snv.output.p2)
  
if (!debug) {
  
unlink(mutCovariate.snv.output.file)

}

rm(mutCovariate_snv_sparse)
gc()

  }
  
}

if (5.4 %in% run.to) {
  
  print("Compile feature matrix for indel for all chromosomes")
  
  if ((!is.null(continuous.features.selected.indel.url.file) | !is.null(discrete.features.selected.indel.url.file) | !is.null(sample.indel.features)) & !is.null(sampled.sites.indel.file)) {
    
mutCovariate_indel_compile = mutCovariate.indel.compile(mask.regions.file = mask.regions.file, all.sites.file = all.sites.file, indel.mutations.file = indel.mutations.int, sample.specific.features.url.file = sample.indel.features, region.of.interest = region.of.interest, cores = cores, indel.mutations.file2 = indel.mutations, chrom.dir = output.dir, genome.build = genome.build)

saveRDS(mutCovariate_indel_compile, file = mutCovariate.indel.output.file)

if (!debug) {
  
delete.files = Sys.glob(paste(output.dir, "mutCovariate_indel_chr*.RDS", sep = ""))
for (i in delete.files) {
  
  unlink(i)
  
}

}

rm(mutCovariate_indel_compile)
gc()

  }
  
}

if (5.5 %in% run.to) {
  
  print("Convert feature matrix for indel to sparse matrix")
  
  if (file.exists(mutCovariate.indel.output.file)) {
  
mutCovariate_indel_sparse = mutCovariate.indel.sparse(compiled = mutCovariate.indel.output.file)

saveRDS(mutCovariate_indel_sparse[[1]], file = mutCovariate.indel.output.p1)
saveRDS(mutCovariate_indel_sparse[[2]], file = mutCovariate.indel.output.p2)
  
if (!debug) {
  
unlink(mutCovariate.indel.output.file)

}

rm(mutCovariate_indel_sparse)
gc()

  }
  
}

## Step 6 ##
LRmodel.snv.file = paste(output.dir, "snv-LRmodel", sep = "")
LRmodel.indel.file = paste(output.dir, "indel-LRmodel", sep = "")

if (6 %in% run.to) {
  
  print("Fit prediction model for SNV/indel")
  
  ## Step 6a ##
  if (file.exists(mutCovariate.snv.output.p1)) {
    
  LRmodel = mutLRFit.snv(mutCovariate.table.snv.file = mutCovariate.snv.output.p1, mutCovariate.count.snv.file = mutCovariate.snv.output.p2, continuous.features.selected.snv.url.file = continuous.features.selected.snv.url.file, discrete.features.selected.snv.url.file = discrete.features.selected.snv.url.file, nucleotide.selected.file = nucleotide.selected.file, sample.specific.features.url.file = sample.snv.features, fit.sparse = fit.sparse, drop = drop, output.dir = output.dir)

  if (class(LRmodel)[1] != "list") {
    
save(LRmodel, file = LRmodel.snv.file)
    
  } else {
    
    if (!"unchanged" %in% LRmodel[[2]]) {
      
      saveRDS(LRmodel[[2]], file = nucleotide.selected.file)
      
    }
    if (!"unchanged" %in% LRmodel[[3]]) {
      
      if (!is.null(LRmodel[[3]])) {
        
      write.table(LRmodel[[3]], file = continuous.features.selected.snv.url.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        
      } else {
        
        z <- data.frame(V1 = character(), V2 = character(), stringsAsFactors = FALSE)
        write.table(z, file = continuous.features.selected.snv.url.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        
      }
      
    }
    if (!"unchanged" %in% LRmodel[[4]]) {
      
      if (!is.null(LRmodel[[4]])) {
        
      write.table(LRmodel[[4]], file = discrete.features.selected.snv.url.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        
      } else {
        
        z <- data.frame(V1 = character(), V2 = character(), stringsAsFactors = FALSE)
        write.table(z, file = discrete.features.selected.snv.url.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        
      }
      
    }
    if (!"unchanged" %in% LRmodel[[5]]) {
      
      if (!is.null(LRmodel[[5]])) {
        
      write.table(LRmodel[[5]], file = sample.snv.features, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
        
      } else {
        
        z <- data.frame(V1 = character(), V2 = character(), stringsAsFactors = FALSE)
        write.table(z, file = sample.snv.features, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        
      }
      
    }
    
    LRmodel = LRmodel[[1]]
    save(LRmodel, file = LRmodel.snv.file)
    
  }

  }
  
## Step 6b ##
  if (file.exists(mutCovariate.indel.output.p1)) {
    
LRmodel = mutLRFit.indel(mutCovariate.table.indel.file = mutCovariate.indel.output.p1, mutCovariate.count.indel.file = mutCovariate.indel.output.p2, continuous.features.selected.indel.url.file = continuous.features.selected.indel.url.file, discrete.features.selected.indel.url.file = discrete.features.selected.indel.url.file, sample.specific.features.url.file = sample.indel.features, fit.sparse = fit.sparse, drop = drop, output.dir = output.dir)

if (class(LRmodel)[1] != "list") {
  
  save(LRmodel, file = LRmodel.indel.file)
  
} else {
  
  if (!"unchanged" %in% LRmodel[[2]]) {
    
    if (!is.null(LRmodel[[2]])) {
      
    write.table(LRmodel[[2]], file = continuous.features.selected.indel.url.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      
    } else {
      
      z <- data.frame(V1 = character(), V2 = character(), stringsAsFactors = FALSE)
      write.table(z, file = continuous.features.selected.indel.url.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      
    }
    
  }
  if (!"unchanged" %in% LRmodel[[3]]) {
    
    if (!is.null(LRmodel[[3]])) {
      
    write.table(LRmodel[[3]], file = discrete.features.selected.indel.url.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      
    } else {
      
      z <- data.frame(V1 = character(), V2 = character(), stringsAsFactors = FALSE)
      write.table(z, file = discrete.features.selected.indel.url.file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      
    }
    
  }
  if (!"unchanged" %in% LRmodel[[4]]) {
    
    if (!is.null(LRmodel[[4]])) {
      
    write.table(LRmodel[[4]], file = sample.indel.features, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
      
    } else {
      
      z <- data.frame(V1 = character(), V2 = character(), stringsAsFactors = FALSE)
      write.table(z, file = sample.indel.features, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
      
    }
  }
  
  LRmodel = LRmodel[[1]]
  save(LRmodel, file = LRmodel.indel.file)
  
}

  }
  
  rm(LRmodel)
  gc()
  
}
  
## Step 7 ##
snv.hotspots = paste(output.dir, "snv_hotspots.tsv", sep = "")
snv.hotspots.merged = paste(output.dir,"snv_hotspots_merged.tsv",sep="")
indel.hotspots = paste(output.dir, "indel_hotspots.tsv", sep = "")
indel.hotspots.merged=paste(output.dir,"indel_hotspots_merged.tsv",sep="")
ann.snv.hotspots = paste(output.dir, "snv_hotspots_annotated.tsv", sep = "")
ann.indel.hotspots = paste(output.dir, "indel_hotspots_annotated.tsv", sep = "")

if (7 %in% run.to) {
  
  print("Predict hotspot mutations for snv/indel")
  
  ## Step 7a ##
  if (file.exists(LRmodel.snv.file)) {
    
  results.snv = mutPredict.snv(mask.regions.file = mask.regions.file, nucleotide.selected.file = nucleotide.selected.file, continuous.features.selected.snv.url.file = continuous.features.selected.snv.url.file, discrete.features.selected.snv.url.file = discrete.features.selected.snv.url.file, sample.specific.features.url.file = sample.snv.features,
                         snv.mutations.file = snv.mutations.int, snv.mutations.file2 = snv.mutations, collapse.regions = collapse.regions, region.of.interest = region.of.interest, cores = cores, snv.model.file = LRmodel.snv.file, min.count = min.count.snv, hotspot.size = hotspot.size, genome.size = genome.size, hotspots = hotspots, merge.hotspots = merge.hotspots, output.dir = output.dir,
                         fdr.cutoff = fdr.cutoff, color.line = color.line, color.dots = color.dots, color.muts = color.muts, top.no = top.no,
                         promoter.file = promoter.file, 
                         utr3.file = utr3.file, utr5.file = utr5.file,
                         other.annotations = other.annotations, debug = debug, genome.build = genome.build)

  if(!is.null(results.snv[[1]])) {
    
write.table(results.snv[[1]], file = snv.hotspots, col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")

  }
  
  if (!is.null(results.snv[[2]])) {
    
    write.table(results.snv[[2]], file = snv.hotspots.merged, col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
    
  }
  
  if (!is.null(results.snv[[3]])) {
    
    write.table(results.snv[[3]], file = ann.snv.hotspots, col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
    
  }
  
  }
  
## Step 7b ##
  if (file.exists(LRmodel.indel.file)) {
    
results.indel = mutPredict.indel(mask.regions.file = mask.regions.file, continuous.features.selected.indel.url.file = continuous.features.selected.indel.url.file, discrete.features.selected.indel.url.file = discrete.features.selected.indel.url.file, sample.specific.features.url.file = sample.indel.features,
                             indel.mutations.file = indel.mutations.int, indel.mutations.file2 = indel.mutations, indel.model.file = LRmodel.indel.file, collapse.regions = collapse.regions, region.of.interest = region.of.interest, cores = cores, min.count = min.count.indel, hotspot.size = hotspot.size, genome.size = genome.size, hotspots = hotspots, merge.hotspots = merge.hotspots, output.dir = output.dir,
                             fdr.cutoff = fdr.cutoff, color.line = color.line, color.dots = color.dots, color.muts = color.muts, top.no = top.no,
                             promoter.file = promoter.file, 
                             utr3.file = utr3.file, utr5.file = utr5.file,
                             other.annotations = other.annotations, debug = debug, genome.build = genome.build)

if (!is.null(results.indel[[1]])) {
  
write.table(results.indel[[1]], file = indel.hotspots, col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
  
}

if (!is.null(results.indel[[2]])) {
  
  write.table(results.indel[[2]], file = indel.hotspots.merged, col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
  
}

if (!is.null(results.indel[[3]])) {
  
  write.table(results.indel[[3]], file = ann.indel.hotspots, col.names = TRUE, row.names = TRUE, quote = FALSE, sep = "\t")
  
}

  }
  
  gc()
  
}


## Optional: dilution analysis
if (dilution.analysis) {
  
dilution.test(snv.mutations.file = snv.mutations, mask.regions.file = mask.regions.file, all.sites.file = all.sites.file,
              cores = cores, cutoff.nucleotide = cutoff.nucleotide, cutoff.nucleotide.new = cutoff.nucleotide.new, cutoff.features = cutoff.features,
              cutoff.features.new.snv = cutoff.features.new.snv, genomic.features.snv = genomic.features.snv, genomic.features.indel = genomic.features.indel,
              genomic.features = genomic.features, genomic.features.fixed.snv = genomic.features.fixed.snv,
              genomic.features.fixed.indel = genomic.features.fixed.indel, genomic.features.fixed = genomic.features.fixed,
              feature.dir = features.dir, mutCovariate.table.snv.file = mutCovariate.snv.output.p1,
              mutCovariate.count.snv.file = mutCovariate.snv.output.p2, 
              continuous.features.selected.snv.url.file = continuous.features.selected.snv.url.file, 
              discrete.features.selected.snv.url.file = discrete.features.selected.snv.url.file, 
              nucleotide.selected.file = nucleotide.selected.file, sample.specific.features.url.file = sample.snv.features,
              output.dir = output.dir, genome.build = genome.build)

}

}
