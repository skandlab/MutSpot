###
### run logistic regression from covariate table;
### Sept 9 2016
### qsub -cwd -pe OpenMP 1 -l mem_free=200G,h_rt=48:00:00 -q medium.q -b y R CMD BATCH mutCovariate_LRFit_hotspot.R
###

setwd("/mnt/projects/guoy1/wgs/Gastric_Cancer/")
source('~/Gastric_Cancer/scripts/mutrec-lib.R')
library(speedglm)
#library(biglm)

## mutation count of each individuals
maf <- maf.to.granges('/mnt/projects/guoy1/wgs/Gastric_Cancer/data/gastric_RF_nonMSI_prefiltered.MAF')
maf.ind=split(maf, maf$sid)
ind.mut.count=sapply(maf.ind, length)
nind=length(ind.mut.count)

###
### Logistic regression 
###
# convert data frame into the input format for logistic regression
mutfreq.aggregated=readRDS(file='/mnt/projects/guoy1/wgs/Gastric_Cancer/data/covariate_genome_freq_table_nonMSI_prefiltered.rds')
#encode sid as mutation count
# this is necessary as some individuals have 0 mututation at CTCF sites ==> fitted probabilities numerically 0 or 1
mutfreq.aggregated$sites.sid=ind.mut.count[mutfreq.aggregated$sites.sid]
mutfreq.aggregated=mutfreq.aggregated[ ,!(names(mutfreq.aggregated) %in% c("tot.count"))]

# runs speeedglm 
# run logistic regression 
#LRmodel<-speedglm(cbind(mut.count,nonmut.count) ~ ., family = binomial(logit), data=mutfreq.aggregated, x=F, y=F,sparse=F)
#LRmodel<-glm(cbind(mut.count,nonmut.count) ~ ., family = binomial(logit), data=mutfreq.aggregated, x=F, y=F)
#save(LRmodel, file="/mnt/projects/guoy1/wgs/Gastric_Cancer/data/LRmodel_gastric_genome_nonMSI_prefiltered")

## remove features that are not predicitive in the SNV hotspot model ("H3K27me3_E094","H2BK15ac_meta","H3K23ac_meta")
mutfreq.aggregated2=mutfreq.aggregated[,!names(mutfreq.aggregated)%in% c("H3K27me3_E094","H2BK15ac_meta","H3K23ac_meta")]
LRmodel<-glm(cbind(mut.count,nonmut.count) ~ ., family = binomial(logit), data=mutfreq.aggregated2, x=F, y=F)
save(LRmodel, file="/mnt/projects/guoy1/wgs/Gastric_Cancer/data/LRmodel_gastric_genome_nonMSI_prefiltered_trunc")
