###
### select features most correlated with mutation rate using the site and regional models
### July 13 2016
### qsub -cwd -pe OpenMP 16 -l mem_free=100G,h_rt=48:00:00 -q medium.q -b y R CMD BATCH --no-save --no-restore lasso_feature_selection.R lasso_feature_selection.intron.Rout
###

#library(fitdistrplus)
library(stats)
library(glmnet)
library(parallel)
library(MASS)

cores= 8

lasso.firstq <- function (x_data, y_data, q, off_data=NULL, family) {
  if (family == "binomial") {
    fit <- glmnet(x_data, y_data, dfmax = q,alpha=1,standardize=TRUE,family="binomial")
  } else if (family == "poisson"){
    fit <- glmnet(x_data, y_data, dfmax = q,alpha=1,offset=off_data,standardize=TRUE,family="poisson")
  }
  #Type "nonzero" returns a list of the indices of the nonzero coefficients for each value of s
  m <- predict(fit, type = "nonzero")
  delta <- q - unlist(lapply(m, length))
  delta[delta < 0] <- Inf
  take <- which.min(delta)
  m=m[[take]]
  # Type "coefficients" computes the coefficients at the requested values for s.
  d <- predict(fit, type = "coefficients")
  d <- d[,take][-1]
  return(list(m=m,d=d))
}

stability.firstq <- function (x_data, y_data, off_data=NULL, family, EV = 1, threshold = 0.75, B = 100, fraction = 0.5, ncores = 1, verbose =TRUE) {
  if (threshold > 1 | threshold < 0.5)
    stop("threshold has to be in (0.5, 1)")
  if (family != "binomial" && family != "poisson")
    stop ("family has to be either binomial or poisson")
  n <- nrow(x_data)
  p <- ncol(x_data)
  col.nam <- colnames(x_data)
  q <- ceiling(sqrt(EV * p * (2 * threshold - 1)))
  sel.mat <- matrix(FALSE, nrow = B, ncol = p)
  sel.n <- floor(fraction * n)
  
  oneSample <- function(bi) {
    if (verbose & (bi %% ceiling(B/20)) == 0) {
      print(paste("Progress : ",bi,"/",B,sep=""))
    }
    
    sel <- sample(1:n, sel.n, replace = FALSE)
    x.sel <- x_data[sel, ]
    y.sel <- x_data[sel]
    if (! is.null(off_data)) {
      off.sel <- off_data[sel]
    } else {
      off.sel <- NULL
    }
    #v2
    sel.model <- lasso.firstq(x_data=x.sel, y_data=y.sel, q=q, off_data=off.sel, family=family)
    sel.feat <- sel.model[[1]]
    sel.coef <- sel.model[[2]]
    #out <- logical(ncol(x))
    out <- logical(p)
    out[sel.feat] <- TRUE
    return(list(sel.coef=sel.coef,out=out))
  }
  
  sel.mat <- matrix(unlist(mclapply(1:B, oneSample, mc.cores = ncores)), nrow = B, byrow = TRUE)
  coef <- sel.mat[,1:p]
  freq <- sel.mat[,(p+1):(2*p)]
  
  freq <- colMeans(freq)
  names(freq) <- col.nam
  out <- list()
  sel.current <- which(freq >= threshold)
  names(sel.current) <- col.nam[sel.current]
  if (length(sel.current) == 0)
    sel.current <- NULL
  colnames(coef) <- col.nam
  #   neg.sign <- apply(coef,2,function(x) sum(sign(x)==-1))
  #   names(neg.sign) <- col.nam
  #   pos.sign <- apply(coef,2,function(x) sum(sign(x)==1))
  #   names(pos.sign) <- col.nam
  
  out <- sel.current
  out <- list(selected = sel.current, coef = coef, EV = EV, threshold = threshold, freq = freq, q = q, method = "stability", call = match.call())
  return(out)
}


## stability selection from full model
lasso.fit <- function(x_data, y_data, off_data=NULL, family){
  if (family == "binomial") {
    cv.fit=cv.glmnet(x_data,y_data,alpha=1,type.measure="deviance",nfolds=10,standardize=TRUE,family="binomial")
    fit.min = glmnet(x_data, y_data, standardize = TRUE, alpha = 1, family = "binomial", lambda = cv.fit$lambda.min)
    fit.1se = glmnet(x_data, y_data, standardize = TRUE, alpha = 1, family = "binomial", lambda = cv.fit$lambda.1se)
  } else if (family == "poisson"){
    cv.fit=cv.glmnet(x_data,y_data,offset=off_data,alpha=1,type.measure="deviance",nfolds=10,standardize=TRUE,family="poisson")
    fit.min = glmnet(x_data, y_data, offset = off_data, standardize = TRUE, alpha = 1, family = "poisson", lambda = cv.fit$lambda.min)
    fit.1se = glmnet(x_data, y_data, offset = off_data, standardize = TRUE, alpha = 1, family = "poisson", lambda = cv.fit$lambda.1se)
  }
  m.min <- predict(fit.min, type = "nonzero")
  d.min <- predict(fit.min, type = "coefficients")[-1]
  m.1se <- predict(fit.1se, type = "nonzero")
  d.1se <- predict(fit.1se, type = "coefficients")[-1]
  results=list(a = m.min, b = d.min, c = m.1se, d = d.1se)
  #str(results)
  return(results)
}

stability.sel <- function(x_data, y_data, off_data=NULL, family, B, fraction, threshold, verbose = TRUE, EV=1, ncores = 1){
  #B = 100, fraction = 0.5, threshold=0.75
  if (family != "binomial" && family != "poisson")
    stop ("family has to be either binomial or possion")
  
  n = nrow(x_data)
  p = ncol(x_data)
  sel.n = floor(fraction*n)
  col.nam = colnames(x_data)
  
  oneSample <- function(bi) {
    if (verbose & (bi %% ceiling(B/100)) == 0) {
      print(paste("Progress : ",bi,"/",B,sep=""))
    }
    
    sel <- sample(1:n, sel.n, replace = FALSE)
    x.sel <- x_data[sel,]
    y.sel <- y_data[sel]
    off.se <- c()
    if (! is.null(off_data)) {
      off.sel <- off_data[sel]
    } else {
      off.sel <- NULL
    }    
    sel.model <- lasso.fit(x_data=x.sel, y_data=y.sel, off_data=off.sel, family=family)
    #sel.model is a dataframe with four columns
    sel.feat.min <- sel.model[[1]]
    sel.coef.min <- sel.model[[2]]
    sel.feat.1se <- sel.model[[3]]
    sel.coef.1se <- sel.model[[4]]
    sel.feat.min = as.vector(sel.feat.min[,1])
    sel.feat.1se = as.vector(sel.feat.1se[,1])
    
    out.min <- logical(p)
    out.min[sel.feat.min] <- TRUE
    out.1se <- logical(p)
    out.1se[sel.feat.1se] <- TRUE
    return(list(sel.coef.min=sel.coef.min,out.min=out.min,sel.coef.1se=sel.coef.1se,out.1se=out.1se))
  }
  
  sel.mat <- matrix(unlist(mclapply(1:B, oneSample, mc.cores = ncores)), nrow = B, byrow = TRUE)
  #print (dim(sel.mat))
  coef.min <- sel.mat[,1:p]
  freq.min <- sel.mat[,(p+1):(2*p)]
  coef.1se <- sel.mat[,(2*p+1):(3*p)]
  freq.1se <- sel.mat[,(3*p+1):(4*p)]
  
  freq.min <- colMeans(freq.min)
  names(freq.min) <- col.nam
  out <- list()
  sel.current.min <- which(freq.min >= threshold)
  names(sel.current.min) <- col.nam[sel.current.min]
  if (length(sel.current.min) == 0)
    sel.current.min <- NULL
  colnames(coef.min) <- col.nam
  
  freq.1se <- colMeans(freq.1se)
  names(freq.1se) <- col.nam
  sel.current.1se <- which(freq.1se >= threshold)
  names(sel.current.1se) <- col.nam[sel.current.1se]
  if (length(sel.current.1se)==0)
    sel.current.1se <- NULL
  colnames(coef.1se) <- col.nam
  
  out <- list(selected.min = sel.current.min, coef.min = coef.min, selected.1se = sel.current.1se, coef.1se = coef.1se, EV = EV, threshold = threshold, freq.min = freq.min, freq.1se = freq.1se, method = "stability", call = match.call())
  return(out)
}

###
### regional model
###
# feat.roi=readRDS('/mnt/projects/guoy1/wgs/Gastric_Cancer/data/features/feat.matrix.intron_RF.rds')
# #remove segway features
# feat.roi=feat.roi[,-c(21,22,31:36, 91:96, 98:122)]
# stabs.regional<-stability.sel(as.matrix(feat.roi[,4:ncol(feat.roi)]),feat.roi$mut.count, off_data=log(feat.roi$roi.width),family="poisson", EV=1,threshold=0.75,B=100,fraction=0.5,ncores=cores,verbose=TRUE)
# saveRDS(stabs.regional, file='/mnt/projects/guoy1/wgs/Gastric_Cancer/data/features/lasso_stability_region_intron_RF_noSegwayDeconvBroad.rds')
#roi.1se=data.frame(coef=colMeans(stabs.regional$coef.1se[,stabs.regional$selected.1se]),freq=stabs.regional$freq.1se[stabs.regional$selected.1se])
#roi.1se[order(roi.1se$freq, decreasing=T),]
#roi.min=data.frame(coef=colMeans(stabs.regional$coef.min[,stabs.regional$selected.min]),freq=stabs.regional$freq.min[stabs.regional$selected.min])
#roi.min[order(roi.min$freq, decreasing=T),]


# feat.site=readRDS('/mnt/projects/guoy1/wgs/Gastric_Cancer/data/features/feat.matrix.all_nonMSI_sites.rds')
# #remove decovoluted expression and histone features, repeated Broad histone marks and segway features
# feat.site=feat.site[,-c(19,20,29:34,89:94, 96:120)]
feat.site=readRDS('/mnt/projects/guoy1/wgs/Gastric_Cancer/data/features/feat.matrix.all_nonMSI_sites_prefiltered.rds')
feat.site=feat.site[,-c(10:15,70:75, 77:101)]
stabs.site<-stability.sel(as.matrix(feat.site[,-1]),feat.site[,1], family="binomial", EV=1,threshold=0.75,B=100,fraction=0.5,ncores=cores,verbose=TRUE)
saveRDS(stabs.site, file='/mnt/projects/guoy1/wgs/Gastric_Cancer/data/features/lasso_stability_site_all_nonMSI_noSegwayDeconvBroad_prefiltered.rds')
#site.1se=data.frame(coef=colMeans(stabs.site$coef.1se[,stabs.site$selected.1se]),freq=stabs.site$freq.1se[stabs.site$selected.1se])
#site.1se[order(site.1se$freq, decreasing=T),]
#site.min=data.frame(coef=colMeans(stabs.site$coef.min[,stabs.site$selected.min]),freq=stabs.site$freq.min[stabs.site$selected.min])
#site.min[order(site.min$freq, decreasing=T),]