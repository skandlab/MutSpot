#' 1 iteration of LASSO selection on a subset of data.
#' 
#' @param x_data Covariates matrix.
#' @param y_data Response matrix.
#' @param s Seed for random rumber.
#' @return Nucleotide contexts selected for 1 iteration of LASSO.
#' @export

lasso.fit.nucleotide <- function(x_data, y_data, s) {
  
  set.seed(s)
  cv.fit = glmnet::cv.glmnet(x_data, as.matrix(y_data), alpha = 1, type.measure = "deviance", nfolds = 10, standardize = FALSE, family = "binomial")
  cv.sel.feat <- as.data.frame(as.matrix(coef(cv.fit, s = "lambda.1se")))
  d.1se <- cv.sel.feat[-1, 1]
  cv.sel.feat = cv.sel.feat[-1, ]
  m.1se <- which(cv.sel.feat != 0)
  
  results = list(a = m.1se, b = d.1se)
  # str(results)
  
  return(results)
  
}
