#' Plot manhattan.
#'
#' @param hotspots.file Hotspots generated.
#' @param fdr.cutoff FDR cutoff, default = 0.05.
#' @param color.line Color given FDR cutoff, default = red.
#' @param color.dots Color hotspots that passed given FDR cutoff, default = maroon1.
#' @return Manhattan figure.
#' @export

plot_manhattan = function(hotspots.file, fdr.cutoff = 0.05, color.line = "red", color.dots = "maroon1") {
  
  hotspots.plot = hotspots.file
  hotspots.plot$region = rownames(hotspots.plot)
  
    x = hotspots.plot
    x$chrom = as.character(x$chrom)
    x$CHR = ifelse(x$chrom == "chrX", "23", substr(x$chrom, 4, nchar(x$chrom)))
    x$CHR = as.numeric(x$CHR)
    x$bp = ceiling((x$start + x$end) / 2)
    x$transcript_id = x$region
    
  
  highlight = x[which(x$fdr <= fdr.cutoff), "transcript_id"]
  if (length(highlight) == 0) {
    
    highlight = NULL
    fdr.cutoff = NULL
    
  } else {
    
    fdr.cutoff = mean(c(max(x[which(x$fdr <= fdr.cutoff), "pval"]), min(x[which(x$fdr > fdr.cutoff), "pval"])))
    
  }
  
  # Modified manhattan function
  chr = "CHR"
  bp = "bp"
  p = "pval"
  snp = "transcript_id"
  fdr = "fdr"
  col = c("gray9", "gray49")
  chrlabs = c(1:22, "X")
  logp = TRUE
  annotatePval = NULL
  annotateTop = TRUE
  div = 1
  
  CHR = BP = P = FDR = index = NULL
   
  if (!(chr %in% names(x))) stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(x))) stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(x))) stop(paste("Column", p, "not found!"))
  if (!(fdr %in% names(x))) stop(paste("Column", fdr, "not found!"))
  if (!is.numeric(x[[chr]])) stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers"))
  if (!is.numeric(x[[bp]])) stop(paste(bp, "column should be numeric."))
  if (!is.numeric(x[[p]])) stop(paste(p, "column should be numeric."))
  if (!is.numeric(x[[fdr]])) stop(paste(fdr, "column should be numeric."))
  
  d = data.frame(transcript_id = x[["transcript_id"]], CHR = x[[chr]], BP = x[[bp]], P = x[[p]], FDR = x[[fdr]])
  
  if (!is.null(x[[snp]])) d = transform(d, SNP = x[[snp]])
  
  d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P) & is.numeric(FDR)))
  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    
    d$logp <- -log10(d$P)
    
  } else {
    
    d$logp <- d$P
  }
  
  d$pos = NA
  
  d$index = NA
  ind = 0
  for (i in unique(d$CHR)) {
    
    ind = ind + 1
    d[d$CHR == i, ]$index = ind
    
  }
  
  nchr = length(unique(d$CHR))
  if (nchr == 1) {
    
    d$pos = d$BP
    ticks = floor(length(d$pos)) / 2 + 1
    xlabel = paste('Chromosome', unique(d$CHR), 'position')
    labs = ticks
    
  } else { 
    
    lastbase = 0
    ticks = NULL
    for (i in unique(d$index)) {
      
      if (i == 1) {
        
        d[d$index == i, ]$pos = d[d$index == i, ]$BP
        
      } else {
        
        lastbase = lastbase + tail(subset(d, index == i - 1)$BP, 1)
        d[d$index == i, ]$pos = d[d$index == i, ]$BP + lastbase
        
      }
      ticks = c(ticks, (min(d[d$index == i, ]$pos) + max(d[d$index == i, ]$pos)) / 2 + 1)
      
    }
    xlabel = 'Chromosome'
    labs <- unique(d$CHR)
  }
  
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * (-0.03))
  
  plot(runif(10), runif(10), 
       xlim = c(xmin, xmax), ylim = c(0, ceiling(max(d$logp))), 
       axes = FALSE, # Don't plot the axis 
       type = "n",  # hide the points
       ylab = expression(-log[10](italic(p))), xlab = xlabel)
  
  axis(2, seq(0, ceiling(max(d$logp)), div))
  
  if (!is.null(chrlabs)) {
    
    if (is.character(chrlabs)) {
      
      if (length(chrlabs) == length(labs)) {
        
        labs <- chrlabs
        
      } else {
        
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
        
      }
      
    } else {
      
      warning ("If you're trying to specify chromosome labels, chrlabs must be a character vector")
      
    }
    
  }
  
  if (nchr == 1) {
    
    axis(1, at=ticks, labels = labs, las = 2, cex.axis = 0.7)
    
  } else {
    
    axis(1, at = ticks, labels = labs, las = 2, cex.axis = 0.7)
    
  }
  
  col = rep(col, max(d$CHR))
  
  if (nchr == 1) {
    
    with(d, points(pos, logp, pch = 16, col = col[1]))
    
  } else {
    
    icol = 1
    for (i in unique(d$index)) {
      
      with(d[d$index == unique(d$index)[i], ], points(pos, logp, col = col[icol],  pch = 16))
      icol = icol + 1
      
    }
    
  }
  
  par(xpd = FALSE)
  
  # Highlight snps from a character vector
  if (!is.null(highlight)) {
    
    if (any(!(highlight %in% d$transcript_id))) warning("You're trying to highlight SNPs that don't exist in your results.")
    d.highlight = d[which(d$transcript_id %in% highlight), ]
    with(d.highlight, points(pos, logp, col = color.dots, pch = 16))
    
  }
  
  if (!is.null(fdr.cutoff)) {
    
  abline(h = (-log10(fdr.cutoff)), col = color.line)
    
  }
  
}
