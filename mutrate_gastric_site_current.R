###
### compute site specific epigenetic features of 10 bp windows around each site in regulatory region
### include all mutated sites and sample non-mutated sites at regular intervals
### July 13 2016
### qsub -cwd -pe OpenMP 16 -l mem_free=100G,h_rt=48:00:00 -q medium.q -b y R CMD BATCH --no-save --no-restore mutrate_gastric_site_current.R 
### 

setwd("/mnt/projects/guoy1/wgs/Gastric_Cancer/")
source('~/Gastric_Cancer/scripts/mutrec-lib.R')
library(BSgenome.Hsapiens.UCSC.hg19)
library(poibin)
cores=16

chrOrder<-c(paste("chr",1:22,sep=""),"chrX")
seqi = seqinfo(Hsapiens)[seqnames(Hsapiens)[1:23]]


## Read in mutated and non-mutated promoter sites
sites=readRDS('/mnt/projects/guoy1/wgs/Gastric_Cancer/data/features/all_nonMSI_sites_prefiltered_downsampled.rds')
seqlevels(sites)=as.character(unique(seqnames(sites)))
# prom.sites=readRDS('/mnt/projects/guoy1/wgs/Gastric_Cancer/data/features/prom_sites_RF_downsampled.rds')
# prom.sites$type ="Prom"
# intron.sites=readRDS('/mnt/projects/guoy1/wgs/Gastric_Cancer/data/features/intron_sites_RF_downsampled.rds')
# intron.sites$type ="Intron"
# enh.sites=readRDS('/mnt/projects/guoy1/wgs/Gastric_Cancer/data/features/enh_sites_RF_downsampled.rds')
# enh.sites$type ="Enh"
# utr3.sites=readRDS('/mnt/projects/guoy1/wgs/Gastric_Cancer/data/features/UTR3_sites_RF_downsampled.rds')
# utr3.sites$type ="UTR3"
# utr5.sites=readRDS('/mnt/projects/guoy1/wgs/Gastric_Cancer/data/features/UTR5_sites_RF_downsampled.rds')
# utr5.sites$type ="UTR5"
# sites=c(prom.sites,intron.sites,enh.sites,utr3.sites,utr5.sites)
# seqlevels(sites)=as.character(unique(seqnames(sites)))
# prom<-sites$type=="Prom" 
# intron<-sites$type=="Intron" 
# enh<-sites$type=="Enh" 
# utr3<-sites$type=="UTR3" 
# utr5<-sites$type=="UTR5" 

###
### Get Tan Expression data
###
exp.summarize.bins <- function(exp,rois) {
  locus.m=do.call(rbind, strsplit(as.character(exp$locus),":|-"))
  locus=data.frame('chr'=locus.m[,1],'start'=as.numeric(locus.m[,2]), 'end'=as.numeric(locus.m[,3])) 
  fpkm.gr <- with(locus, GRanges(chr, IRanges(start, end), mcols=exp[,-(1:2)]))
  names(mcols(fpkm.gr))= names(exp[,-(1:2)])
  fpkm.gr=fpkm.gr[seqnames(fpkm.gr) %in% seqnames(seqi)]
  seqlevels(fpkm.gr)=as.character(unique(seqnames(fpkm.gr)))
  hits=findOverlaps(rois, fpkm.gr)
  hits.m=as.matrix(hits)
  rois.fpkm=data.frame(mut_region=hits.m[,"queryHits"], mcols(fpkm.gr[hits.m[,"subjectHits"]]))
  # find regions with no overlaps with expression data
  no.hit=which(!(c(1:length(rois))%in%hits.m[,1]))
  no.hit.m= data.frame(mut_region=no.hit, matrix(0, nrow=length(no.hit), ncol=ncol(exp)-2))
  colnames(no.hit.m)=colnames(rois.fpkm)  
  rois.fpkm=rbind(rois.fpkm, no.hit.m)
  # for regions that overlap more than 1 transcripts, take the max value for each feature
  rois.fpkm.aggregated= aggregate(.~mut_region, FUN=max, data=rois.fpkm)
  # order regions by row names
  rois.fpkm.aggregated=rois.fpkm.aggregated[order(rois.fpkm.aggregated$mut_region),]
  rois.fpkm.aggregated[,-1]  
}

# decov.expr.df <- read.table('/mnt/projects/skanderupamj/wgs/tan_gastric/bw_50bp/nnls/expr_dcv_may25.tsv', header=T)
# decov.expr.df$gene=rownames(decov.expr.df)
# decov.expr.df=decov.expr.df[,c("gene","stroma.avg","cancer.avg","tumor.med","normal.med")]
# decov.expr.df$stroma.avg[decov.expr.df$stroma.avg<0]<-0
# decov.expr.df$cancer.avg[decov.expr.df$cancer.avg<0]<-0
# cellLine.expr.df= read.table('/mnt/projects/guoy1/wgs/Gastric_Cancer/data/Tan_cellLine_expression_normalized.txt', header=T, sep="\t")
# expression.df=merge(cellLine.expr.df, decov.expr.df, by="gene", all=F)
# expression.features = exp.summarize.bins(expression.df, sites)

###
### Fetch local mutation rate
###
local_mutrate=bigwig.summarize.bins("/mnt/projects/guoy1/wgs/Gastric_Cancer/data/mutrate.local.100kb.RF.nonMSI.bigWig", sites, 'mean')


###
### Fetch Replication data
###
mean_rep_time=bigwig.summarize.bins("/mnt/projects/guoy1/wgs/Gastric_Cancer/Mixed_model/RepTimeWig/wgEncodeUwRepliSeqWaveSignalMean.bigWig", sites, 'mean')


###
### Fetch ensembl reg. build data
###
ensembl.urls <- read.table('/mnt/projects/guoy1/wgs/Gastric_Cancer/data/ensembl_all_urls_local_dnase.txt')
ensembl.features = mclapply(ensembl.urls[,2],function(f) {bigwig.summarize.bins(f,sites)}, mc.cores=cores, mc.preschedule = FALSE, mc.silent = FALSE)
ensembl.features = do.call(cbind,ensembl.features)
colnames(ensembl.features) = as.character(ensembl.urls[,1])
# clean up cache in /tmp
# unlink('/tmp/udcCache/http/ngs.sanger.ac.uk/production/ensembl/regulation/hg19',recursive=TRUE)


###
### Fetch Roadmap Epigenome data
###
roadmap.urls <- read.table('/mnt/projects/guoy1/wgs/Gastric_Cancer/data/roadmap_all_urls.txt', header=F, sep="\t")
roadmap.gastric.urls=roadmap.urls[roadmap.urls[,2]%in% c("E094", "E092", "E110", "E111"),] #27 features
roadmap.gastric.features= mclapply(roadmap.gastric.urls[,3],function(f) {bigwig.summarize.bins(f,sites, 'mean')}, mc.cores=4, mc.preschedule = FALSE, mc.silent = FALSE)
roadmap.gastric.features= do.call(cbind,roadmap.gastric.features)
colnames(roadmap.gastric.features) <- paste(roadmap.gastric.urls[,1], roadmap.gastric.urls[,2], sep="_")

roadmap.meta.urls <- read.table('/mnt/projects/guoy1/wgs/Gastric_Cancer/data/roadmap_meta_urls.txt', header=F, sep="\t")
gastric.marks=c("H3K27ac", "H3K36me3", "H3K4me3","DNase","H3K27me3","H3K9me3","H3K4me1","H3K9ac")
roadmap.meta.urls=roadmap.meta.urls[!roadmap.meta.urls[,1]%in% gastric.marks,]
roadmap.meta.features=  mclapply(roadmap.meta.urls[,3],function(f) {bigwig.summarize.bins(f,sites, 'mean')}, mc.cores=2, mc.preschedule = FALSE, mc.silent = FALSE)
roadmap.meta.features= do.call(cbind,roadmap.meta.features)
colnames(roadmap.meta.features) <- paste(roadmap.meta.urls[,1], roadmap.meta.urls[,2], sep="_")
#sort( sapply(ls(),function(x){object.size(get(x))}))

###
### Tan Epigenetic data updated
###
tan.epi.urls <- read.table('/mnt/projects/guoy1/wgs/Gastric_Cancer/data/Tan_Epi_urls_update3.txt', header=F, sep="\t")
tan.epi.features = mclapply(tan.epi.urls[,3],function(f) {bigwig.summarize.bins(f,sites)}, mc.cores=cores, mc.preschedule = FALSE, mc.silent = FALSE)
tan.epi.features= do.call(cbind,tan.epi.features)
colnames(tan.epi.features) <- paste(tan.epi.urls[,1], tan.epi.urls[,2], sep="_")
tan.epi.features[tan.epi.features<0]=0

###
### Fetch Broad histone data
###
histone.features = c()
histone.urls.lines <- readLines('./data/broad_histone_urls.txt')
histone.urls.names = unlist(lapply(histone.urls.lines,function(x) {strsplit(x,"\t")[[1]][1]}))
histone.urls = lapply(histone.urls.lines,function(x) {strsplit(x,"\t")[[1]][-1]})
names(histone.urls) = histone.urls.names
histone.urls=histone.urls[!names(histone.urls)%in% gastric.marks]
# summarize over features and multiple cell lines
for (i in 1:length(histone.urls)) {
  feature.name <- names(histone.urls)[i]
  print(paste0('>>',feature.name))
  #run in parallel
  scores = mclapply(histone.urls[[i]],function(f) {bigwig.summarize.bins(f, sites)}, mc.cores=cores, mc.preschedule = FALSE, mc.silent = FALSE)
  scores = do.call(cbind,scores)
  # average/median over all files / cell types
  mean.scores = apply(scores,1,mean)
  histone.features = cbind(histone.features,mean.scores)
}
# names
colnames(histone.features) <- paste(names(histone.urls), 'Broad', sep="_")


# make feature matrix for epigenomics data
# feat.matrix=data.frame("mut"=sites$mut, "local.mutrate"=local_mutrate, "mean.rep.time"=mean_rep_time, expression.features, tan.epi.features, roadmap.gastric.features, roadmap.meta.features, histone.features, ensembl.features)
# No gene expression feature for hotspot
feat.matrix=data.frame("mut"=sites$mut, "local.mutrate"=local_mutrate, "mean.rep.time"=mean_rep_time, tan.epi.features, roadmap.gastric.features, roadmap.meta.features, histone.features, ensembl.features)
saveRDS(feat.matrix, file='/mnt/projects/guoy1/wgs/Gastric_Cancer/data/features/feat.matrix.all_nonMSI_sites_prefiltered.rds')
# write Rdata objects for each roi
# feat.matrix.prom=feat.matrix[prom,]
# feat.matrix.intron=feat.matrix[intron,]
# feat.matrix.enh=feat.matrix[enh,]
# feat.matrix.5utr=feat.matrix[utr5,]
# feat.matrix.3utr=feat.matrix[utr3,]
# saveRDS(feat.matrix.prom, file='/mnt/projects/guoy1/wgs/Gastric_Cancer/data/features/feat.matrix.prom_site_RF.rds')
# saveRDS(feat.matrix.intron, file='/mnt/projects/guoy1/wgs/Gastric_Cancer/data/features/feat.matrix.intron_site_RF.rds')
# saveRDS(feat.matrix.enh, file='/mnt/projects/guoy1/wgs/Gastric_Cancer/data/features/feat.matrix.enh_site_RF.rds')
# saveRDS(feat.matrix.5utr, file='/mnt/projects/guoy1/wgs/Gastric_Cancer/data/features/feat.matrix.5utr_site_RF.rds')
# saveRDS(feat.matrix.3utr, file='/mnt/projects/guoy1/wgs/Gastric_Cancer/data/features/feat.matrix.3utr_site_RF.rds')
