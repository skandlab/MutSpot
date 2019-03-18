#' Plot feature importance.
#' 
#' @param model Prediction model.
#' @param mutCovariate.table.file RDS file, SNV covariates sparse matrix.
#' @param mutCovariate.count.file RDS file, SNV response matrix.
#' @param continuous.features.selected.url.file Text file containing URLs of SNV continuous features selected for model.
#' @param z.value To use z-value for plot or coefficients, default = FALSE.
#' @param mutation.type SNV or indel mutation.
#' @param output.dir Save plot in given output directory.
#' @return Variable importance plot.
#' @export

plot_feature_importance = function(model, mutCovariate.table.file, mutCovariate.count.file, continuous.features.selected.url.file, z.value = FALSE, mutation.type, output.dir) {

    load(model)
  if (!"glmnet" %in% class(LRmodel)) {
    
  t = summary(LRmodel)
  df = coef(t)
  df = as.data.frame(df)
  df$feat = rownames(df)
  df = df[-1, ]
  colnames(df) = c("estimate", "std.error", "z.value", "p.value", "feat")

  
  df = df[order(df$z.value, decreasing = TRUE), ]
  df$feat = factor(df$feat, levels = df$feat)
  colnames(df)[3] = "imp"
  
  } else {
    
    mutfreq.aggregated = readRDS(file = mutCovariate.table.file)
    mutfreq.aggregated = as.data.frame(as.matrix(mutfreq.aggregated))
    
    if (z.value) {
      
      mutfreq.aggregated2 = readRDS(file = mutCovariate.count.file)
      mutfreq.aggregated = cbind(mutfreq.aggregated2, mutfreq.aggregated)
      
      if (!is.null(continuous.features.selected.url.file)) {
        
        selected.continuous.urls <- read.delim(continuous.features.selected.url.file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
        
      }
      
      for(i in colnames(mutfreq.aggregated)[which(!colnames(mutfreq.aggregated) %in% c("mut.count", "nonmut.count", "sid", selected.continuous.urls[ ,1]))]) {
        
        mutfreq.aggregated[ ,i] = as.character(mutfreq.aggregated[ ,i])
        
      }
      
      LRmodel2 <- stats::glm(cbind(mut.count,nonmut.count) ~ ., family = binomial(logit), data = mutfreq.aggregated, x = F, y = F)
      
      t = summary(LRmodel2)
      df = coef(t)
      df = as.data.frame(df)
      df$feat = rownames(df)
      df = df[-1, ]
      colnames(df) = c("estimate", "std.error", "z.value", "p.value", "feat")
      
      
      df = df[order(df$z.value, decreasing = TRUE), ]
      df$feat = factor(df$feat, levels = df$feat)
      colnames(df)[3] = "imp"
      
    } else {
      
    betas = LRmodel$beta
    
    sds <- apply(as.matrix(mutfreq.aggregated), 2, sd)
    cs <- as.matrix(betas)
    df <- cs[ ,1] * sds
    
    df = data.frame(df)
    df$feat = rownames(df)
    colnames(df)[1] = "imp"
    df = df[order(df$imp, decreasing = TRUE), ]
    df$feat = factor(df$feat, levels = df$feat)

    }
    
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
