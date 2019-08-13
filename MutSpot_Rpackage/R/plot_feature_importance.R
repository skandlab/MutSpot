#' Plot feature importance.
#' 
#' @param LRmodel Prediction model.
#' @param mutCovariate.table.file RDS file, SNV covariates sparse matrix.
#' @param mutation.type SNV or indel mutation.
#' @param output.dir Save plot in given output directory.
#' @param feature.names Feature names.
#' @return Variable importance plot.
#' @export

plot_feature_importance = function(LRmodel, mutCovariate.table.file, mutation.type, output.dir, feature.names) {

  if (!"glmnet" %in% class(LRmodel)) {
    
  t = summary(LRmodel)
  df = coef(t)
  df = as.data.frame(df)
  df$feat = rownames(df)
  df = df[-1, ]
  colnames(df) = c("estimate", "std.error", "z.value", "p.value", "feat")

  # Clean up variable name
  df$feat = gsub("`", "", gsub("`1", "", df$feat))
  df$feat = ifelse(!df$feat %in% feature.names, substr(df$feat, 1, nchar(df$feat) - 1), df$feat)
  
  df = df[order(df$z.value, decreasing = TRUE), ]
  df$feat = factor(df$feat, levels = df$feat)
  colnames(df)[3] = "imp"
  
  } else {
    
    mutfreq.aggregated = readRDS(file = mutCovariate.table.file)
    mutfreq.aggregated = as.data.frame(as.matrix(mutfreq.aggregated))
    
    betas = LRmodel$beta
    
    sds <- apply(as.matrix(mutfreq.aggregated), 2, sd)
    cs <- as.matrix(betas)
    df <- cs[ ,1] * sds
    
    df = data.frame(df)
    df$feat = rownames(df)
    colnames(df)[1] = "imp"
    
    # Clean up variable name
    df$feat = gsub("`", "", gsub("`1", "", df$feat))
    df$feat = ifelse(!df$feat %in% feature.names, substr(df$feat, 1, nchar(df$feat) - 1), df$feat)
    
    df = df[order(df$imp, decreasing = TRUE), ]
    df$feat = factor(df$feat, levels = df$feat)

  }
  
  pdf(paste(output.dir, mutation.type, "_features.pdf", sep = ""))
  print(ggplot2::ggplot(df, ggplot2::aes(x = feat, y = imp)) + ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          axis.line = ggplot2::element_line(colour = "black"))+
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
      ggplot2::xlab("Features") + 
      ggplot2::ylab("Variable importance"))
  dev.off()
  
}
