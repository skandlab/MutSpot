#' Mutation hotspot recurrence prediction for SNV.
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
#' @param merge.hotspots To plot overlapping hotspots as 1 hotspot or individual hotspots, default = TRUE.
#' @param output.dir Save plot in given output directory.
#' @param fdr.cutoff FDR cutoff, default = 0.05.
#' @param color.line Color given FDR cutoff, default = red.
#' @param color.dots Color hotspots that passed given FDR cutoff, default = maroon1.
#' @param color.muts Color points, default = orange.
#' @param top.no Number of top hotspots to plot, default = 3.
#' @param promoter.file Promoter regions bed file, depends on genome build, default file = Ensembl75.promoters.coding.bed, Ch37.
#' @param utr3.file 3'UTR regions bed file, depends on genome build, default file = Ensembl75.3UTR.coding.bed, Ch37.
#' @param utr5.file 5'UTR regions bed file, depends on genome build, default file = Ensembl75.5UTR.coding.bed, Ch37.
#' @param other.annotations Text file containing URLs of additional regions to be annotated, default = NULL.
#' @param debug To delete temporary files or not, default = FALSE.
#' @param genome.build Reference genome build, default = Ch37.
#' @return Dataframe containing predicted hotspots significance with hotspots information for SNV/merged/annotated and hotspots manhattan and top hits figures.
#' @export

mutPredict.snv = function(mask.regions.file = system.file("extdata", "mask_regions.RDS", package = "MutSpot"), nucleotide.selected.file, continuous.features.selected.snv.url.file, discrete.features.selected.snv.url.file, sample.specific.features.url.file = NULL, snv.mutations.file, snv.mutations.file2, collapse.regions = FALSE, region.of.interest, cores = 1, snv.model.file, min.count = 2, hotspot.size = 21, genome.size = 2533374732, hotspots = TRUE, merge.hotspots = TRUE, output.dir,
                          fdr.cutoff = 0.05, color.line = "red", color.dots = "maroon1", color.muts = "orange", top.no = 3,
                          promoter.file = system.file("extdata", "Ensembl75.promoters.coding.bed", package = "MutSpot"),
                          utr3.file = system.file("extdata", "Ensembl75.3UTR.coding.bed", package = "MutSpot"), 
                          utr5.file = system.file("extdata", "Ensembl75.5UTR.coding.bed", package = "MutSpot"), 
                          other.annotations = NULL, debug = FALSE, genome.build = "Ch37") {
  

  mutPredict.snv.prepare(mask.regions.file = mask.regions.file, nucleotide.selected.file = nucleotide.selected.file, 
                         continuous.features.selected.snv.url.file = continuous.features.selected.snv.url.file, 
                         discrete.features.selected.snv.url.file = discrete.features.selected.snv.url.file, 
                         sample.specific.features.url.file = sample.specific.features.url.file, 
                         snv.mutations.file = snv.mutations.file, snv.mutations.file2 = snv.mutations.file2, collapse.regions = collapse.regions,
                         region.of.interest = region.of.interest, 
                         cores = cores, snv.model.file = snv.model.file, min.count = min.count, 
                         hotspot.size = hotspot.size, genome.size = genome.size, hotspots = hotspots, output.dir = output.dir, genome.build = genome.build)

  return(mutPredict.snv.run.lr(output.dir = output.dir, merge.hotspots = merge.hotspots, snv.mutations.file = snv.mutations.file,
                               fdr.cutoff = fdr.cutoff, color.line = color.line, color.dots = color.dots,
                               color.muts = color.muts, top.no = top.no,
                               promoter.file = promoter.file,
                               utr3.file = utr3.file, utr5.file = utr5.file, 
                               other.annotations = other.annotations, debug = debug, cores = cores, genome.build = genome.build))
  
}

