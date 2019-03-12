#' 100 iterations of LASSO selection for epigenetic features.
#' 
#' @param x_data Covariates matrix.
#' @param y_data Response matrix.
#' @param threshold Frequency cutoff/threshold to determine stable nucleotide contexts, ranges from 0.5 to 1, default = 0.75.
#' @param cores Number of cores, default = 1.
#' @return A list containing the frequency of epigenetic features selected through LASSO and the continuous and discrete features that passed the threshold.
#' @export

stability.sel.epigenetic <- function(x_data, y_data, threshold = 0.75, cores = 1) {
  
  fraction = 0.50
  B = 100
  verbose = TRUE
  seed.no1 = seq(5001, 10000, 50)
  seed.no2 = seq(10001, 15000, 50)
  seed.no3 = seq(1, 5000, 50)

  p = ncol(x_data) 
  col.nam = colnames(x_data) 
  
  mut = which(y_data == "1")
  nonmut = which(y_data == "0")
  n1 = length(mut)
  sel.n1 = floor(fraction * n1)
  n2 = length(nonmut)
  sel.n2 = floor(fraction * n2)
  
  oneSample <- function(bi, s1, s2, s3) { 
    
    if (verbose & (bi %% ceiling(B / 100)) == 0) { 
      
      print(paste("Progress : ", bi, "/", B, sep = "")) 
      
    } 
    
    set.seed(s1)
    sel1 <- sample(mut, sel.n1, replace = FALSE) 
    set.seed(s2)
    sel2 <- sample(nonmut, sel.n2, replace=FALSE)
    x.sel <- x_data[c(sel1, sel2), ] 
    y.sel <- y_data[c(sel1, sel2)] 
    
    sel.model <- lasso.fit.epigenetic(x.sel, y.sel, s3) 
    sel.feat.1se <- sel.model[[1]] 
    sel.coef.1se <- sel.model[[2]] 
    out.1se <- logical(p) 
    out.1se[sel.feat.1se] <- TRUE 
    
    return(list(sel.coef.1se = sel.coef.1se, out.1se = out.1se)) 
    
  } 
  
  sel.mat <- matrix(unlist(parallel::mclapply(1:B, FUN = function(k) oneSample(k, seed.no1[k], seed.no2[k], seed.no3[k]), mc.cores = cores)), nrow = B, byrow = TRUE) 
  coef.1se <- sel.mat[ ,1:p] 
  freq.1se <- sel.mat[ ,(p + 1):(2 * p)] 
  freq.1se <- colMeans(freq.1se) 
  names(freq.1se) <- col.nam 
  sel.current.1se <- which(freq.1se >= threshold) 
  if (length(sel.current.1se) == 0) {
    
    sel.current.1se <- NULL 
    
  }
  colnames(coef.1se) <- col.nam 
  
  out <- list(selected.1se = sel.current.1se, coef.1se = coef.1se, threshold = threshold, freq.1se = freq.1se)
  
  return(out) 
  
} 
