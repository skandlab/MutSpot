#' Plot top hotspots.
#'
#' @param hotspots.file Hotspots generated.
#' @param fdr.cutoff FDR cutoff, default = 0.05.
#' @param color.muts Color points, default = orange.
#' @param mutations.file Mutations found in region of interest MAF file.
#' @param mutation.type SNV or indel mutations.
#' @param top.no Number of top hotspots to plot, default = 3.
#' @param output.dir Save plots in given output directory.
#' @return Top hotspots figures.
#' @export

plot_top_hits = function(hotspots.file, fdr.cutoff = 0.05, color.muts ="orange", mutations.file, mutation.type, top.no = 3, output.dir) {
  
  hotspots.plot = hotspots.file
  hotspots.plot$region = rownames(hotspots.plot)
  
    x = hotspots.plot
    top.hits = x[which(x$fdr <= fdr.cutoff), ]
    
  if (nrow(top.hits) >=1 ) {
    
  top.hits.gr = with(top.hits, (GenomicRanges::GRanges(chrom, IRanges::IRanges(start, end))))
  mutations = maf.to.granges(mutations.file)
  
  if (top.no > nrow(top.hits)) {
    
    top.no = nrow(top.hits)
    
  }
  
  # Lollipop plot for each top hit
  for (i in 1:top.no) {
    
    print(paste("Plotting Hit", i, sep = " "))
    
    data = data.frame(position = seq(top.hits$start[i], top.hits$end[i], by = 1))
    
    ovl = IRanges::findOverlaps(top.hits.gr[i], mutations)
    mut.hits = mutations[S4Vectors::subjectHits(ovl)]
    mut.hits = GenomicRanges::as.data.frame(mut.hits)
    mut.hits = aggregate(sid ~ start, mut.hits, FUN = function(k) length(unique(k)))
    colnames(mut.hits) = c("position", "mut.count")
    data = merge(data, mut.hits, by = "position", all = TRUE)
      
    pdf(paste(output.dir, mutation.type, "_hotspot_", i, ".pdf", sep = ""))
    suppressWarnings(print(ggplot2::ggplot(data, ggplot2::aes(x = position, y = mut.count)) +
      ggplot2::geom_segment(ggplot2::aes(x = position, xend = position, y = 0, yend = mut.count), color = "grey") +
      ggplot2::geom_point(color = color.muts, size = 4) +
      ggplot2::theme_light() +
      ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(), panel.grid.minor.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(), panel.grid.minor.y = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black")) +
      ggplot2::xlab("Mutation position") +
      ggplot2::ylab("Number of mutated samples")+
      ggplot2::ggtitle(paste(top.hits$chrom[i], ":",top.hits$start[i], "-", top.hits$end[i], sep = "")) +
      ggplot2::scale_x_continuous(breaks = seq(min(data$position), max(data$position), by = 1), limits = c(data$position[1], data$position[nrow(data)])) +
        ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, max(data[!is.na(data$mut.count),"mut.count"]) + 0.1))+
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))))
    dev.off()
    
  }
  
  } else {
  
    print("No significant hotspots")
    
}

}
