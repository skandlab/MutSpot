#' Runs all or selected steps in MutSpot analysis.
#'
#' @param run.to Numeric vector defining which steps to run, default = 1, 2, 3.1, 3.2, 4.1, 4.2, 5.1, 5.2, 5.3, 5.4, 5.5, 6, 7, 8, 9.1, 9.2, 9.3.
#' @param chromosomes Character vector defining which chromosomes to compute feature matrix on, default = chr1-chrX.
#' @param snv.mutations SNV mutations MAF file.
#' @param indel.mutations Indel mutations MAF file.
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
#' @return Corresponding output from each step in MutSpot analysis.
#' @export

MutSpot = function(run.to = c(1:2, 3.1, 3.2, 4.1, 4.2, 5.1, 5.2, 5.3, 5.4, 5.5, 6:8, 9.1, 9.2, 9.3), chromosomes = c(1:22,"X"), snv.mutations, indel.mutations, mask.regions.file = system.file("extdata", "mask_regions.RDS", package = "MutSpot"), all.sites.file = system.file("extdata", "all_sites.RDS", package = "MutSpot"), region.of.interest = NULL, ratio = 1, sample = T, cores = 1, cutoff.nucleotide = 0.90, cutoff.nucleotide.new = NULL, genomic.features.snv = NULL, genomic.features.indel = NULL, genomic.features = NULL, genomic.features.fixed.snv = NULL, genomic.features.fixed.indel = NULL, genomic.features.fixed = NULL, sample.snv.features = NULL, sample.indel.features = NULL, cutoff.features = 0.75, cutoff.features.new.snv = NULL, cutoff.features.new.indel = NULL, fit.sparse = FALSE, drop = FALSE, min.count.snv = 2, min.count.indel = 2, genome.size = 2533374732, hotspots = TRUE,
                  promoter.file = system.file("extdata", "Ensembl75.promoters.coding.bed", package = "MutSpot"), utr3.file = system.file("extdata", "Ensembl75.3UTR.coding.bed", package = "MutSpot"), utr5.file = system.file("extdata", "Ensembl75.5UTR.coding.bed", package = "MutSpot"), other.annotations = NULL, fdr.cutoff = 0.1, color.line = "red", color.dots = "maroon1", merge.hotspots = TRUE, color.muts = "orange", z.value = FALSE, top.no = 3) {
  
## check format of output directory ##
if (substr(output.dir, nchar(output.dir), nchar(output.dir)) != "/") {
  
  output.dir = paste(output.dir, "/", sep = "")
  
}  

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
local.mutrate.snv.file = paste(output.dir, "localmutrate_snv.bed", sep = "")
local.mutrate.indel.file = paste(output.dir, "localmutrate_indel.bed", sep = "")

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
    
    for (i in paste("chr", chromosomes, sep="")) {
      
mutCovariate.snv(mask.regions.file = mask.regions.file, all.sites.file = all.sites.file, nucleotide.selected.file = nucleotide.selected.file, continuous.features.selected.snv.url.file = continuous.features.selected.snv.url.file, 
                 discrete.features.selected.snv.url.file = discrete.features.selected.snv.url.file, sample.specific.features.url.file = sample.snv.features, snv.mutations.file = snv.mutations.int, region.of.interest = region.of.interest, snv.mutations.file2 = snv.mutations, cores = cores, chr.interest = i)

    }
    
  }
  
  ## Step 5b ##
  if ((!is.null(continuous.features.selected.indel.url.file) | !is.null(discrete.features.selected.indel.url.file) | !is.null(sample.indel.features)) & !is.null(sampled.sites.indel.file)) {
    
    for (i in paste("chr", chromosomes, sep="")) {
      
      mutCovariate.indel(mask.regions.file = mask.regions.file, all.sites.file = all.sites.file, continuous.features.selected.indel.url.file = continuous.features.selected.indel.url.file, 
                         discrete.features.selected.indel.url.file = discrete.features.selected.indel.url.file, sample.specific.features.url.file = sample.indel.features, indel.mutations.file = indel.mutations.int, region.of.interest = region.of.interest, indel.mutations.file2 = indel.mutations, cores = cores, chr.interest = i)
      
    }
    
  }
  
}
  
if (5.2 %in% run.to) {
  
  print("Compile feature matrix for snv for all chromosomes")
  
  if ((!is.null(continuous.features.selected.snv.url.file) | !is.null(discrete.features.selected.snv.url.file) | !is.null(sample.snv.features)) & !is.null(sampled.sites.snv.file)) {
  
mutCovariate_snv_compile = mutCovariate.snv.compile(mask.regions.file = mask.regions.file, all.sites.file = all.sites.file, snv.mutations.file = snv.mutations.int, sample.specific.features.url.file = sample.snv.features, region.of.interest = region.of.interest, cores = cores, snv.mutations.file2 = snv.mutations)

saveRDS(mutCovariate_snv_compile, file = mutCovariate.snv.output.file)

  }
  
}

if (5.3 %in% run.to) {
  
  print("Convert feature matrix for snv to sparse matrix")
  
  if (file.exists(mutCovariate.snv.output.file)) {
    
mutCovariate_snv_sparse = mutCovariate.snv.sparse(compiled = mutCovariate.snv.output.file)

saveRDS(mutCovariate_snv_sparse[[1]], file = mutCovariate.snv.output.p1)
saveRDS(mutCovariate_snv_sparse[[2]], file = mutCovariate.snv.output.p2)
  
  }
  
}

if (5.4 %in% run.to) {
  
  print("Compile feature matrix for indel for all chromosomes")
  
  if ((!is.null(continuous.features.selected.indel.url.file) | !is.null(discrete.features.selected.indel.url.file) | !is.null(sample.indel.features)) & !is.null(sampled.sites.indel.file)) {
    
mutCovariate_indel_compile = mutCovariate.indel.compile(mask.regions.file = mask.regions.file, all.sites.file = all.sites.file, indel.mutations.file = indel.mutations.int, sample.specific.features.url.file = sample.indel.features, region.of.interest = region.of.interest, cores = cores, indel.mutations.file2 = indel.mutations)

saveRDS(mutCovariate_indel_compile, file = mutCovariate.indel.output.file)

}

}

if (5.5 %in% run.to) {
  
  print("Convert feature matrix for indel to sparse matrix")
  
  if (file.exists(mutCovariate.indel.output.file)) {
  
mutCovariate_indel_sparse = mutCovariate.indel.sparse(compiled = mutCovariate.indel.output.file)

saveRDS(mutCovariate_indel_sparse[[1]], file = mutCovariate.indel.output.p1)
saveRDS(mutCovariate_indel_sparse[[2]], file = mutCovariate.indel.output.p2)
  
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
    
    plot_feature_importance(model = LRmodel.snv.file, mutCovariate.table.file = mutCovariate.snv.output.p1, mutCovariate.count.file = mutCovariate.snv.output.p2, continuous.features.selected.url.file = continuous.features.selected.snv.url.file, z.value = z.value, mutation.type = "SNV")
    
  }
  
  ## Step 9b ##
  if (file.exists(LRmodel.indel.file)) {
    
    plot_feature_importance(model = LRmodel.indel.file, mutCovariate.table.file = mutCovariate.indel.output.p1, mutCovariate.count.file = mutCovariate.indel.output.p2, continuous.features.selected.url.file = continuous.features.selected.indel.url.file, z.value = z.value, mutation.type = "indel")
    
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
    
    plot_top_hits(hotspots.file = snv.hotspots.merged, fdr.cutoff = fdr.cutoff, color.muts = color.muts, mutations.file  = snv.mutations, mutation.type = "SNV", top.no = top.no)
    
  } else if (file.exists(snv.hotspots)) {
    
    plot_top_hits(hotspots.file = snv.hotspots, fdr.cutoff = fdr.cutoff, color.muts = color.muts, mutations.file  = snv.mutations, mutation.type = "SNV", top.no = top.no)
    
    
  }
  
  ## Step 9b
  if (file.exists(indel.hotspots.merged)) {
    
    plot_top_hits(hotspots.file = indel.hotspots.merged, fdr.cutoff = fdr.cutoff, color.muts = color.muts, mutations.file = indel.mutations, mutation.type = "indel", top.no = top.no)
    
  } else if (file.exists(indel.hotspots)) {
    
    plot_top_hits(hotspots.file = indel.hotspots, fdr.cutoff = fdr.cutoff, color.muts = color.muts, mutations.file = indel.mutations, mutation.type = "indel", top.no = top.no)
    
  }
  
}

}
