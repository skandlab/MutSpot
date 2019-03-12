#' 100 iterations of LASSO selection for nucleotide contexts.
#' 
#' @param x_data Covariates matrix.
#' @param y_data Response matrix.
#' @param threshold Frequency cutoff/threshold to determine stable nucleotide contexts, ranges from 0.5 to 1, default = 0.90.
#' @param cores Number of cores, default = 1.
#' @return A list containing the frequency of nucleotide contexts selected through LASSO and the contexts that passed the threshold.
#' @export

stability.sel.nucleotide <- function(x_data, y_data, threshold = 0.90, cores = 1) {
  
  seed.no1 = seq(5001, 10000, 50)
  seed.no2 = seq(10001, 15000, 50)
  B = 100
  fraction = 0.5
  verbose = TRUE
  
  if (threshold > 1 | threshold < 0.5) {
    
    stop("threshold has to be in (0.5, 1)")
  
  }
  
  n = nrow(x_data)
  p = ncol(x_data)
  sel.n = floor(fraction * n)
  col.nam = colnames(x_data)
  
  sel.mat <- matrix(FALSE, nrow = B, ncol = p)
  
  oneSample <- function(bi, s1, s2) {
    
    if (verbose & (bi %% ceiling(B / 100)) == 0) {
      
      print(paste("Progress : ", bi, "/", B, sep = ""))
      
    }
    
    set.seed(s1)
    sel <- sample(1:n, sel.n, replace = FALSE)
    x.sel <- x_data[sel, ]
    y.sel <- y_data[sel, ]
    sel.model <- lasso.fit.nucleotide(x.sel, y.sel, s2)
    sel.feat.1se <- sel.model[[1]] 
    sel.coef.1se <- sel.model[[2]]
    
    out.1se <- logical(p)
    out.1se[sel.feat.1se] <- TRUE
    
    return(list(sel.coef.1se = sel.coef.1se, out.1se = out.1se))
    
  }
  
  sel.mat <- matrix(unlist(parallel::mclapply(1:B, FUN = function(k) oneSample(k, seed.no1[k], seed.no2[k]), mc.cores = cores)), nrow = B, byrow = TRUE)
  coef.1se <- sel.mat[ ,1:p]
  freq.1se <- sel.mat[ ,(p + 1):(2 * p)]
  
  freq.1se <- colMeans(freq.1se)
  names(freq.1se) <- col.nam
  sel.current.1se <- which(freq.1se >= threshold)
  
  if (length(sel.current.1se) == 0) {
    
    sel.current.1se <- NULL
    
  }
  colnames(coef.1se) <- col.nam
  
  out <- list(selected.1se = sel.current.1se, coef.1se = coef.1se, threshold = threshold, freq.1se = freq.1se, method = "stability")
  
  return(out)
  
}
