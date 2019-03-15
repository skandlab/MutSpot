#' Hotspots annotation.
#'
#' @param hotspots.file Hotspots to be annotated tsv file.
#' @param promoter.file Promoter regions bed file, default file = Ensembl75.promoters.coding.bed.
#' @param utr3.file 3'UTR regions bed file, default file = Ensembl75.3UTR.coding.bed.
#' @param utr5.file 5'UTR regions bed file, default file = Ensembl75.5UTR.coding.bed.
#' @param other.annotations Text file containing URLs of additional regions to be annotated, default = NULL.
#' @return Annotated hotspots.
#' @export

mutAnnotate = function(hotspots.file, promoter.file = system.file("extdata", "Ensembl75.promoters.coding.bed", package = "MutSpot"),
                       utr3.file = system.file("extdata", "Ensembl75.3UTR.coding.bed", package = "MutSpot"), 
                       utr5.file = system.file("extdata", "Ensembl75.5UTR.coding.bed", package = "MutSpot"), 
                       other.annotations = NULL) {
  
  hotspots.ann = read.delim(hotspots.file, header = TRUE, stringsAsFactors = FALSE)
  hotspots.ann$mut_region = rownames(hotspots.ann)
  
  mut.regions <- with(hotspots.ann, GenomicRanges::GRanges(chrom, IRanges::IRanges(start, end), pval = pval, length = length, p.bg = p.bg, k = k, fdr = fdr))
  names(mut.regions) = rownames(hotspots.ann)
  
  # roi.cds <- bed.to.grangeslist(cds.file)
  # roi.transcripts = range(roi.cds)
  # roi.transcripts = BiocGenerics::unlist(roi.transcripts)
  # roi.transcripts = roi.transcripts + 10000
  
  roi.prom <- bed.to.grangeslist(promoter.file)
  roi.3utr <- bed.to.grangeslist(utr3.file)
  roi.5utr <- bed.to.grangeslist(utr5.file)
  # roi.intron <- bed.to.grangeslist(intron.file)

  # Gene hits
  # gene.hits = IRanges::findOverlaps(mut.regions, roi.transcripts)
  # if (length(gene.hits) != 0) {
  #   
  # gene.hits.m = S4Vectors::as.matrix(gene.hits)
  # gene.names = data.frame(gene.hits.m, mut_region = names(mut.regions[gene.hits.m[ ,"queryHits"]]), transcript_id = names(roi.transcripts[gene.hits.m[ ,"subjectHits"]]))
  # gene.names2 = aggregate(transcript_id ~ mut_region, FUN=paste, collapse = "/", data = gene.names)
  # 
  # } else { 
  #   
  #   gene.names2 = data.frame(mut_region = hotspots.ann$mut_region, transcript_id = NA)
  #   
  # }
  # gene.names2$mut_region = as.character(gene.names2$mut_region)
  # colnames(gene.names2)[2] = "gene"
  
  # Promoter hits
  prom.hits = IRanges::findOverlaps(mut.regions, roi.prom)
  if (length(prom.hits) != 0) {
    
  prom.hits.m = S4Vectors::as.matrix(prom.hits)
  prom.names = data.frame(prom.hits.m, mut_region = names(mut.regions[prom.hits.m[ ,"queryHits"]]), transcript_id = names(roi.prom[prom.hits.m[ ,"subjectHits"]]))
  prom.names2 = aggregate(transcript_id ~ mut_region, FUN = paste, collapse = ";", data = prom.names)
  
  } else {
    
    prom.names2 = data.frame(mut_region = hotspots.ann$mut_region, transcript_id = NA)
    
  }
  prom.names2$mut_region = as.character(prom.names2$mut_region)
  colnames(prom.names2)[2] = "promoter"
  
  # 3'UTR hits
  utr3.hits = IRanges::findOverlaps(mut.regions, roi.3utr)
  if (length(utr3.hits) != 0) {
  
  utr3.hits.m = S4Vectors::as.matrix(utr3.hits)
  utr3.names = data.frame(utr3.hits.m, mut_region = names(mut.regions[utr3.hits.m[ ,"queryHits"]]), transcript_id = names(roi.3utr[utr3.hits.m[ ,"subjectHits"]]))
  utr3.names2 = aggregate(transcript_id ~ mut_region, FUN = paste, collapse = ";", data = utr3.names)
  
  } else {
    
    utr3.names2 = data.frame(mut_region = hotspots.ann$mut_region, transcript_id = NA)
    
  }
  utr3.names2$mut_region = as.character(utr3.names2$mut_region)
  colnames(utr3.names2)[2] = "3UTR"
  
  # 5'UTR hits
  utr5.hits = IRanges::findOverlaps(mut.regions, roi.5utr)
  if (length(utr5.hits) != 0) {
  
  utr5.hits.m = S4Vectors::as.matrix(utr5.hits)
  utr5.names = data.frame(utr5.hits.m, mut_region = names(mut.regions[utr5.hits.m[ ,"queryHits"]]), transcript_id = names(roi.5utr[utr5.hits.m[ ,"subjectHits"]]))
  utr5.names2 = aggregate(transcript_id ~ mut_region, FUN = paste, collapse = ";", data = utr5.names)
  
  } else {
    
    utr5.names2 = data.frame(mut_region = hotspots.ann$mut_region, transcript_id = NA)
    
  }
  utr5.names2$mut_region = as.character(utr5.names2$mut_region)
  colnames(utr5.names2)[2] = "5UTR"
  
  # Intron hits  
  # intron.hits = IRanges::findOverlaps(mut.regions, roi.intron)
  # if (length(intron.hits) != 0) {
  # 
  # intron.hits.m = S4Vectors::as.matrix(intron.hits)
  # intron.names = data.frame(intron.hits.m, mut_region = names(mut.regions[intron.hits.m[ ,"queryHits"]]), transcript_id = names(roi.intron[intron.hits.m[ ,"subjectHits"]]))
  # intron.names2 = aggregate(transcript_id ~ mut_region, FUN = paste, collapse = ";", data = intron.names)
  # 
  # } else {
  #   
  #   intron.names2 = data.frame(mut_regions = hotspots.ann$mut_region, transcript_id = NA)
  #   
  # } 
  # intron.names2$mut_region = as.character(intron.names2$mut_region)
  # colnames(intron.names2)[2] = "intron"
  
  # names.list = list(hotspots.ann, gene.names2, prom.names2, utr3.names2, utr5.names2, intron.names2)
  names.list = list(hotspots.ann, prom.names2, utr3.names2, utr5.names2)
  df.out = merged.data.frame = Reduce(function(...) merge(..., by = "mut_region", all.x = T), names.list)

  if (!is.null(other.annotations)) {
    
    other.ann = read.delim(other.annotations, stringsAsFactors = FALSE, header = FALSE)
    
    for (i in 1:nrow(other.ann)) {
    
      print(other.ann[i,1])
      
      roi.other <- bed.to.grangeslist(other.ann[i, 2])
      other.hits = IRanges::findOverlaps(mut.regions, roi.other)
      if (length(other.hits) != 0) {
        
        other.hits.m = S4Vectors::as.matrix(other.hits)
        other.names = data.frame(other.hits.m, mut_region = names(mut.regions[other.hits.m[ ,"queryHits"]]), transcript_id = names(roi.other[other.hits.m[ ,"subjectHits"]]))
        other.names2 = aggregate(transcript_id ~ mut_region, FUN = paste, collapse = ";", data = other.names)
        
      } else {
        
        other.names2 = data.frame(mut_regions = hotspots.ann$mut_region, transcript_id = NA)
        
      } 
      other.names2$mut_region = as.character(other.names2$mut_region)
      colnames(other.names2)[2] = other.ann[i, 1]
        
      names.list = list(df.out, other.names2)
      df.out = merged.data.frame = Reduce(function(...) merge(..., by = "mut_region", all.x = T), names.list)
      
    }
    
  } 
  
  df.out = df.out[order(df.out$pval, decreasing = FALSE), ]
  rownames(df.out) = df.out$mut_region
  df.out = df.out[ ,-1]
  
  return(df.out)
  
}
