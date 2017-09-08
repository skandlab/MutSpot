###
### Make mutation counts and covariate table to input into logistic regression
### Sept 13 2016
### qsub -cwd -pe OpenMP 16 -l mem_free=480G,h_rt=48:00:00 -q medium.q -b y R CMD BATCH mutrec_logistic_hotspot.R
### 

setwd("/mnt/projects/guoy1/wgs/Gastric_Cancer/")
source('~/Gastric_Cancer/scripts/mutrec-lib.R')
library(stringr)
library(speedglm)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings) #for reverseComplement function

cores = 16
chrOrder<-c(paste("chr",1:22,sep=""),"chrX")
seqi = seqinfo(Hsapiens)[seqnames(Hsapiens)[1:23]]
seqnames=seqnames(seqinfo(Hsapiens))[1:23]

####
#### functions
####
# store signals overlapping regions of interest as RLELists
subset.bigwig <- function(bigWigUrl, bins) {
  print(paste0('>>',bigWigUrl))
  ss.cov = import(as.character(bigWigUrl),as='RleList')
  ss.cov = ss.cov[intersect(names(ss.cov),seqlevels(bins))]  
  stopifnot(is(bins, "GRanges"))
  stopifnot(is(ss.cov, "RleList"))
  stopifnot(identical(sort(seqlevels(bins)), sort(names(ss.cov))))
  bins.per.chrom <- split(ranges(bins), seqnames(bins))
  means.list <- lapply(seqlevels(bins),
                       function(seqname) {
                         # subset RLE list with Granges, set the rest of the values to zero
                         rle.subset=Rle(0, length(ss.cov[[seqname]]))
                         rle.subset[bins.per.chrom[[seqname]]] <- ss.cov[[seqname]][bins.per.chrom[[seqname]]]
                         rle.subset
                       })
  # set compress=F for sum(elementLengths(x)) >2^32 to prevent integer overflow
  means.list<-RleList(means.list, compress=F)
  names(means.list)<-seqlevels(bins)
  means.list
}

#get motif and nucleotide context matrix
find.motif=function(seq, motifs){
  
  #get 1mer, 3mer and 5mer
  seq$oneMer=substr(seq$dna,26,26)
  seq$threeMer=substr(seq$dna,25,27)
  seq$fiveMer=substr(seq$dna,24,28) 
  
  #get rev comp for 1mer, 3mer and 5mer
  seq$onerc<-as.character(reverseComplement(DNAStringSet(seq$oneMer))) #rev comp for onemer
  seq$threerc<-as.character(reverseComplement(DNAStringSet(seq$threeMer))) #rev comp for threemer
  seq$fiverc<-as.character(reverseComplement(DNAStringSet(seq$fiveMer))) #rev comp for fivemer
  
  #get pairs for 1mer, 3mer and 5mer
  seq$one.pair=ifelse(seq$oneMer %in% c("A","G"),seq$oneMer,seq$onerc) #oneMer and onerc pair
  seq$three.pair=ifelse(substr(seq$threeMer,2,2) %in% c("A","G"),seq$threeMer,seq$threerc) #threeMer and threerc pair
  seq$five.pair=ifelse(substr(seq$fiveMer,3,3) %in% c("A","G"),seq$fiveMer,seq$fiverc) #fiveMer and fiverc pair
  
  #get right and left flanks for 3mer and 5mer
  seq$three.right=substr(seq$three.pair,2,3) #right/3' flank of three.pair
  seq$three.left=substr(seq$three.pair,1,2) #left/5' flank of three.pair
  seq$five.right=substr(seq$five.pair,3,5) #right/3' flank of five.pair
  seq$five.left=substr(seq$five.pair,1,3) #left/5' flank of five.pair
  
  #get nucleotide context features
  seq$tM=ifelse(seq$three.pair %in% motifs$threeMer,seq$three.pair,"0")
  seq$tM=factor(seq$tM, levels=c("0", motifs$threeMer))
  #n+1 levels in tM - n three.pair selected by lambda.1se
  seq$oM=ifelse(seq$one.pair %in% motifs$oneMer,seq$one.pair,"0")
  seq$oM=factor(seq$oM, levels=c("0", motifs$oneMer))
  #n+! levels in oM - n one.pair selected by lambda.1se
  seq$fM=ifelse(seq$five.pair %in% motifs$fiveMer,seq$five.pair,"0")
  seq$fM=factor(seq$fM, levels=c("0", motifs$fiveMer))
  #n levels in fM - n five.pair selected by lambda.1se
  seq$tR=ifelse(seq$three.right %in% motifs$threeRight,seq$three.right,"0")
  seq$tR=factor(seq$tR, levels=c("0", motifs$threeRight))
  #n+1 levels in tR - n three.right selected by lambda.1se
  seq$tL=ifelse(seq$three.left %in% motifs$threeLeft,seq$three.left,"0")
  seq$tL=factor(seq$tL, levels=c("0", motifs$threeLeft))
  #n+1 levels in tL - n three.left selected by lambda.1se
  seq$fR=ifelse(seq$five.right %in% motifs$fiveRight,seq$five.right,"0")
  seq$fR=factor(seq$fR, levels=c("0", motifs$fiveRight))
  #n+1 levels in fR - n five.right selected by lambda.1se
  seq$fL=ifelse(seq$five.left %in% motifs$fiveLeft,seq$five.left,"0")
  seq$fL=factor(seq$fL, levels=c("0", motifs$fiveLeft))
  #n+1 levels in fL - n five.left selected by lambda.1se
  
  #get motif features
  if (!is.null(motifs$motif$logo)){
    for (i in 1:(length(motifs$motif$logo))){
      a=motifs$motif$a[i]
      motif.logo=motifs$motif$logo[[i]]
      motif.name=motifs$motif$name[i]
      seq[,motif.name]=sapply(seq$dna,function(x){
        if(sum(str_detect(substr(x,(26-a),(26+a)),motif.logo))>=1)
        {1} else {0}
      })
    }
  }
  
  t=NULL
  seq=seq[,20:ncol(seq)]
  colnames(seq)[1:7]=c("threeMer","oneMer","fiveMer","threeRight","threeLeft","fiveRight","fiveLeft")
  # remove nucleotide context categories that are not selected  
  seq=seq[,names(seq) %in% c(names(motifs),motifs$motif$name)]
  return(seq)
}


## get feature matrix given roi 
get.features <- function(roi, features, motifs){
  ## Get a list of sites for each region
  sites=lapply(roi, function(f) {GRanges(seqnames(f), IRanges(seq(start(f),end(f),1), seq(start(f),end(f),1)))})
  sites=do.call(c, sites)
  
  # remove sites with N in 5mers
  seq = GRanges(seqnames(sites), IRanges(ranges(sites)@start-25, ranges(sites)@start+25))
  seq=Views(Hsapiens,seq)
  seq=as.data.frame(seq)   
  if (sum(grepl("N",seq$dna))!=0){
    sites.del=which(grepl("N",seq$dna))
    seq=seq[-sites.del,] #remove those with N in 5mers
    sites=sites[-sites.del]
  }
  
  ## For each roi, make a covariate table. nrows= roi length, ncols=number of covariates
  roi.features1 <- lapply(features,function(f) {as.numeric(f[[as.character(unique(seqnames(sites)))]][ranges(sites)])})
  
  ## discretize all features
  roi.features=sapply(1:length(roi.features1),function(i) {
    if (class(discrete[[i]])=="numeric") {
      out=ifelse(roi.features1[[i]]>discrete[i],1,0)
    }else{
      feat.bins=.bincode(roi.features1[[i]], discrete[[i]][[1]])
      out=sapply(feat.bins, function(x) discrete[[i]][[2]][x])
    }})
  #print (roi.features)
  roi.features=data.frame(roi.features)
  colnames(roi.features)=names(roi.features1)
  
  roi.motif<-find.motif(seq, motifs)
  out=cbind(roi.features,roi.motif)
  return (out)
}

## calculate regional p-values for regions of interest, based on background mutation rates predicted by logistic regression
run.mutrec.lr <- function(roi, maf, model, features, motifs, discrete, min.count=min.count, genome.size=3*10**9, cores=1, only.samples=NULL,debug=FALSE) {
  
  # first check that roi names are unique
  if (any(duplicated(names(roi)))) {
    print("Error: regions with duplcated names")
    return()
  }
  
  print(paste("Using ",length(roi)," ROIs from input",sep=""))
  
  if(debug) {
    print(roi)
  }
  
  # filter maf if sample list specified
  if(!is.null(only.samples)) {
    print("Filtering MAF based on supplied sample list")
    maf = maf[which(maf$sid %in% only.samples)]
  }
  
  # we will overlap with maf many times, so we create a GIntervalTree
  print(">> Intersecting ROIs and MAF ...")
  #maf.gtree = GIntervalTree(maf)
  maf.ovl <- findOverlaps(roi,maf,ignore.strand=TRUE)
  maf.ovl.m = as.matrix(maf.ovl)
  
  if (debug) {
    print(maf.ovl.m)
  }
  
  ## 1) select rois mutated in min.count.sample
  nsamples = length(unique(maf$sid))
  #roi.count.mutated = tapply(elementMetadata(maf[subjectHits(maf.ovl),'sid'])[,1],names(roi)[queryHits(maf.ovl)],function(muts) length(unique(muts)))
  # encode sample names as integers to increase the efficiency of the code
  roi.count.mutated = tapply(as.numeric(maf$sid[maf.ovl.m[,2]]),names(roi)[maf.ovl.m[,1]],function(s) length(unique(s)))
  roi.mutated = names(roi.count.mutated)[roi.count.mutated>=min.count]
  print(paste("Using ",length(roi.mutated)," ROIs mutated in >=",min.count," samples",sep=""))
  
  # run analysis for all mutated roi
  rois.to.run = roi.mutated #[1:24]
  
  if(debug) {
    print(rois.to.run)
  }
  
  results=c()
  
  if (length(rois.to.run) == 0) {
    return(results)
  }
  
  process.single.roi <- function(x) {
    # print progress for every 1% of data
    roi.progress <- which(rois.to.run == x)
    if ((roi.progress %% ceiling(length(rois.to.run)/100)) == 0) {
      print(paste("Progress : ",roi.progress,"/",length(rois.to.run),sep=""))
    }
    # make sure we are only using first hit
    # should not have multiple hits in new version, where we check for duplicate roi names
    x.idx = which(names(roi) == x)[1] # in roi
    
    # extract features for all sites in roi[x], note roi here is a GRange object, not GRangeList
    roi.feat=get.features(unlist(reduce(roi[x.idx])), features, motifs)
    x.len = nrow(roi.feat)
    
    # vector of patient IDs
    sid=unique(maf$sid)
    
    # compute probability of mutation in region for all samples; vector of 212 probabilities
    p.bg=sapply (sid, function(s) {
      sid=rep(s, x.len)      
      # encode sid as the mutation count of each individual    
      sid=rep(ind.mut.count[[s]], x.len)
      df=cbind(roi.feat, sites.sid=sid)
      # compute background mutation rate foreach site in each individual
      p=predict(object=model,newdata=df, type="response")  
      #print(p[1:10])
      # compute background mutation rate of the region in each individual
      p.roi=1-prod(1-p)      
      return(p.roi)})
    
    # k, number of samples where x is mutated
    if (debug) {
      print(x)
      print(x.idx)
      print(maf.ovl.m[,1]==x.idx)
    }
    
    k = length(unique(maf$sid[maf.ovl.m[maf.ovl.m[,1]==x.idx,2]]))
    # Poisson binomial
    # pval = 1-ppoibin(k-1,p.dist,method='RNA') # fast
    # pval = 1-ppoibin(k-1,p_dist,method='DFT-CF') # exact
    # trick to get other tail of distr.
    # only RF provides estimates for low p-values
    pval = ppoibin(length(p.bg)-k,1-p.bg,method='RF') # fast
    # pval<-tryCatch(
    #   ppoibin(length(p.bg)-(k-1),1-p.bg,method='RF'),
    #   warning=function(x) {
    #     message(x)
    #     message(p.bg)
    #     return(NA)
    #   },
    #   error=function(x) {
    #     message(x)
    #     message(p.bg)
    #     return(NA)
    #   })    
    c(x,pval,x.len,mean(p.bg),k)
  }
  
  results = mclapply(rois.to.run,function(xr) {process.single.roi(xr)}, mc.cores=cores, mc.preschedule = FALSE, mc.silent = FALSE)
  
  #error checking and remove entries with error
  results.error=sapply(results, function(x) {length(x)==1})
  print(rois.to.run[which(results.error)])
  print(results[results.error])
  results=results[!results.error]
  
  results = do.call(rbind,results)
  rownames(results) = results[,1]
  # drop=FALSE, force matrix if single row / vector
  results = results[,-1,drop=FALSE]
  class(results) <- "numeric"
  # include ignored regions
  fdr = p.adjust(results[,1],method='BH')*genome.size/nrow(results)
  fdr=sapply(fdr, function(x) min(x,1))  
  results = cbind(results, fdr)
  colnames(results) = c('pval','length','p.bg','k','fdr')
  results = results[order(results[,'pval']),,drop=FALSE]
}

###
###   Main
###
min.count=3

mappability=import("/mnt/projects/guoy1/wgs/Gastric_Cancer/data/wgEncodeCrgMapabilityAlign75mer.bigWig")
## define reads that can map to more than 1 genomic location as non-mappable
nonmappable=mappability[mappability$score<1,]
# convert zero-based coordinates to one-based coordinates
nonmappable=shift(nonmappable,1)
nonmappable= reduce(nonmappable)
nonmappable=nonmappable[seqnames(nonmappable) %in% seqnames(seqi)]
seqlevels(nonmappable)=as.character(unique(seqnames(nonmappable)))

## mask CDS and ig loci
roi.cds <- bed.to.granges('/mnt/projects/guoy1/wgs/Gastric_Cancer/roi_Ensembl75/Ensembl75.CDS.bed')
roi.cds.ext <- reduce(roi.cds + 5) # extend each region with +/- 5 bases and get all non-overlapping regions
#immunoglobulin loci
ig.loci <- bed.to.granges('/mnt/projects/guoy1/wgs/Gastric_Cancer/data/ig_loci.bed')
ig.loci <- reduce(ig.loci + 10**5) # extend each region with 100kb and get all non-overlapping regions
# combine mask regions
# one of the ig.loci located on sequence chr14 is out-of-bound, trim at this stage as there is no seqinfo associated with ig.loci and roi.cds.ext
mask.regions=reduce(trim(c(ig.loci,roi.cds.ext,nonmappable))) 

maf <- maf.to.granges('/mnt/projects/guoy1/wgs/Gastric_Cancer/data/gastric_RF_nonMSI_prefiltered.MAF')
maf=maf[seqnames(maf) %in% seqnames(seqi)]
seqlevels(maf)=as.character(unique(seqnames(maf)))
# remove 5 samples with oxidative damage
maf=maf[!maf$sid%in%c("tan2001206", "tan20021007", "tan980319", "tan2000986", "tan980436")]
maf.ind=split(maf, maf$sid)
ind.mut.count=sapply(maf.ind, length)

# mask mutations
maf.masked <-maf[-subjectHits(findOverlaps(mask.regions,maf))]
dupl = duplicated(maf.masked)
maf.uniq= maf.masked[!dupl,]

# extend each mutation with +/- 10 bases
mut.regions=maf.uniq+10
names(mut.regions)=paste("Mutation", c(1:length(mut.regions)), sep="")

rm(list=c("maf.uniq", "mask.regions", "nonmappable", "roi.cds.ext", "roi.cds", "ig.loci"))
gc()

load(file='/mnt/projects/guoy1/wgs/Gastric_Cancer/data/LRmodel_gastric_genome_nonMSI_prefiltered_trunc')

# sel.motif=list()
# ## selected motifs
# sel.motif$threeMer=c("AAG","AGG","AGT","AAA","TGG","CGT")
# sel.motif$threeLeft=c("CA","GA","CG")
# sel.motif$threeRight=c("AC","AA","GC")
# sel.motif$fiveLeft=c("AGA","TCA","TGG")
# sel.motif$fiveRight=c("GCA","AAG","AGT")
# motifs.count=5

sel.motif=list()
## selected motifs
sel.motif$oneMer=c("A")
sel.motif$threeMer=c("AAG","AAC")
sel.motif$threeLeft=c("CA","AA","CG")
sel.motif$threeRight=c("GA")
sel.motif$fiveLeft=c("TAA","AAG","TTG")
sel.motif$fiveRight=c("AAG","AGT","AGA")
motifs.count=6

## Read in covariate data
selected.urls <- read.table('/mnt/projects/guoy1/wgs/Gastric_Cancer/data/hotspot_selected_urls_nonMSI_prefiltered_trunc.txt')
selected.features = mclapply(selected.urls[,2],function(f) {subset.bigwig(f, mut.regions)}, mc.cores=4, mc.preschedule = FALSE, mc.silent = FALSE)
names(selected.features) = as.character(selected.urls[,1])

## discretization
# read in replication timing bins
repBins=readRDS(file='/mnt/projects/guoy1/wgs/Gastric_Cancer/data/repTime_bins.rds')
repBinMean=readRDS(file='/mnt/projects/guoy1/wgs/Gastric_Cancer/data/repTime_bin_mean.rds')
# read in local mutation rate bins
localMutRateBins=readRDS(file='/mnt/projects/guoy1/wgs/Gastric_Cancer/data/localMutRate_nonMSI_bins.rds')
localMutRateBinMean=readRDS(file='/mnt/projects/guoy1/wgs/Gastric_Cancer/data/localMutRate_nonMSI_bin_mean.rds')

# list of discretization cutoffs
discrete=list(
  localMutRate=list(localMutRateBins, unlist(localMutRateBinMean)),
  repTime=list(repBins, unlist(repBinMean)),
  0.25,4,4,4,4,4,4,4,0,0
)

# Run hotspot recurrence analysis
genome.size= 2533374732
mut.rec <- run.mutrec.lr(mut.regions, maf, model=LRmodel, features=selected.features, motifs=sel.motif, discrete=discrete, min.count=min.count,genome.size=genome.size, cores=cores)
mut.regions2=mut.regions[names(mut.regions) %in% rownames(mut.rec)]
mut.rec.hotspot= data.frame(chrom=as.character(seqnames(mut.regions2[rownames(mut.rec)])), start=start(mut.regions2[rownames(mut.rec)]), end=end(mut.regions2[rownames(mut.rec)]), mut.rec)

write.table(mut.rec.hotspot,'/mnt/projects/guoy1/wgs/Gastric_Cancer/mutrec_lr_results/LRmodel_hotspot_nonMSI_prefiltered-5.tsv',sep="\t",col.names=TRUE,quote=FALSE,row.names=TRUE)
gc()

# sig.hotspots=c("Mutation7658", "Mutation7659", "Mutation30017", "Mutation30018", "Mutation4946", "Mutation27136","Mutation30275", "Mutation30276", "Mutation47189",
#                "Mutation47190", "Mutation10820", "Mutation6783", "Mutation7841")
# muts=lapply(sig.hotspots, function(f) {
#   ovl=findOverlaps(mut.regions[f], maf)
#   maf[subjectHits(ovl)]})
# names(muts)=sig.hotspots
# muts
