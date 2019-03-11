#' Prepare covariate sparse matrix and response matrix.
#' 
#' @param compiled Aggregated table RDS file, default = mutCovariate-compile-part1.RDS.
#' @return A list containing covariates sparse matrix and response matrix.
#' @export

mutCovariate.indel.sparse = function(compiled = "mutCovariate-indel-compile-part1.RDS") {
  
aggregated.table = readRDS(compiled)

aggregated.table = as.data.frame(aggregated.table)
counts = aggregated.table[ ,c("mut.count", "tot.count", "indel.length", "nonmut.count")]
counts = counts[ ,-which(colnames(counts) %in% c("tot.count", "indel.length"))]
features = colnames(aggregated.table)
features = features[!features %in% c("mut.count", "tot.count", "indel.length", "nonmut.count")]
aggregated.table = aggregated.table[ ,features]

aggregated.table = apply(aggregated.table, 2, FUN = function(x) Matrix::Matrix(x, sparse = TRUE))
aggregated.table = do.call(cbind, aggregated.table)
colnames(aggregated.table) = features

return(list(aggregated.table, counts))

}
