###
### Make genome-wide mutation counts and covariate table for hotspot analysis
### Sept 8 2016
### qsub -cwd -pe OpenMP 16 -l mem_free=480G,h_rt=48:00:00 -q medium.q -b y R CMD BATCH --no-save --no-restore mutCovariate_hotspot.R mutCovariate_hotspot_nonMSI_prefiltered.Rout
###

### 2 bugs fixed Dec 6 2016
## bug1: sid= rep(unique(mut.masked$sid),nrow(genome.freq.aggregated)) => sid= rep(unique(mut.masked$sid),each=nrow(genome.freq.aggregated))
# replicate each sid n times instead of replicating the entire sid vector n times
## bug2: in find.motif function, use seq$dna instead of seq[6], as the dna sequence not always in the 6th column
## optimised the code by using tile() function to get a list of individual sites instead of lapply over all regions
## optimised the find.motifs() function by vectorizing operations and using the reverseComplement() function

setwd("/mnt/projects/guoy1/wgs/Gastric_Cancer/")
library(data.table)
source('~/Gastric_Cancer/scripts/mutrec-lib.R')
library("lsr")
library(BSgenome.Hsapiens.UCSC.hg19)
library("plyr") #count function
library(stringr)
library(Biostrings) #for reverseComplement function
cores=16

chrOrder<-c(paste("chr",1:22,sep=""),"chrX")
seqi = seqinfo(Hsapiens)[seqnames(Hsapiens)[1:23]]
seqnames=seqnames(seqinfo(Hsapiens))[1:23]

## output file name
outfile = '/mnt/projects/guoy1/wgs/Gastric_Cancer/data/covariate_genome_freq_table_nonMSI_prefiltered.rds'

###
### import nonmappable regions
###
mappability=import("/mnt/projects/guoy1/wgs/Gastric_Cancer/data/wgEncodeCrgMapabilityAlign75mer.bigWig")
## define reads that can map to more than 1 genomic location as non-mappable
nonmappable=mappability[mappability$score<1,]
# convert zero-based coordinates to one-based coordinates
nonmappable=shift(nonmappable,1)
nonmappable= reduce(nonmappable)
nonmappable=nonmappable[seqnames(nonmappable) %in% seqnames(seqi)]
seqlevels(nonmappable)=as.character(unique(seqnames(nonmappable)))

# trim ends of each chromosome that are N's
mappability.grl=split(mappability, seqnames(mappability))
ranges=lapply(mappability.grl, range)
ranges=GRangesList(ranges)
ranges=unlist(ranges)
extrSeq=Views(Hsapiens,ranges) # check that all N's at chromosome ends are trimmed
# trim the ends
start(ranges)=pmax(26, start(ranges))

## mask CDS and ig loci
roi.cds <- bed.to.granges('/mnt/projects/guoy1/wgs/Gastric_Cancer/roi_Ensembl75/Ensembl75.CDS.bed')
roi.cds.ext <- reduce(roi.cds + 5) # extend each region with +/- 5 bases and get all non-overlapping regions
#immunoglobulin loci
ig.loci <- bed.to.granges('/mnt/projects/guoy1/wgs/Gastric_Cancer/data/ig_loci.bed')
ig.loci <- reduce(ig.loci + 10**5) # extend each region with 100kb and get all non-overlapping regions
# combine mask regions
# one of the ig.loci located on sequence chr14 is out-of-bound, trim at this stage as there is no seqinfo associated with ig.loci and roi.cds.ext
mask.regions=reduce(trim(c(ig.loci,roi.cds.ext,nonmappable))) 

# Create Genomic Ranges for the chromosomes 1:22 and chrX
all.sites=ranges[seqnames(ranges)%in% seqnames(seqi)]
#all.sites=ranges[seqnames(ranges)%in% c("chr21","chr22")]
seqlevels(all.sites)=as.character(unique(seqnames(all.sites)))
sum(as.numeric(width(all.sites))) #2951331845
all.sites.masked=subtract.regions.from.roi(all.sites, mask.regions, cores=cores) #3961533
sum(as.numeric(width(all.sites.masked))) #2648698330

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
  0.25,4,4,4,4,4,4,4,4,4,0,4,0
)

# gastric
# sel.motif=list()
# ## selected motifs
# sel.motif$threeMer=c("AAG","AGG","AGT","AAA","TGG","CGT")
# sel.motif$threeLeft=c("CA","GA","CG")
# sel.motif$threeRight=c("AC","AA","GC")
# sel.motif$fiveLeft=c("AGA","TCA","TGG")
# sel.motif$fiveRight=c("GCA","AAG","AGT")
# motifs.count=5

# gastric non-MSI
sel.motif=list()
## selected motifs
sel.motif$oneMer=c("A")
sel.motif$threeMer=c("AAG","AAC")
sel.motif$threeLeft=c("CA","AA","CG")
sel.motif$threeRight=c("GA")
sel.motif$fiveLeft=c("TAA","AAG","TTG")
sel.motif$fiveRight=c("AAG","AGT","AGA")
motifs.count=6

# read covariate file paths
selected.urls <- read.table('/mnt/projects/guoy1/wgs/Gastric_Cancer/data/hotspot_selected_urls_nonMSI_prefiltered.txt')
#selected.urls=selected.urls[1:3,]

###
### functions
###
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


###
### function that creates a frequency table of mutated and non-mutated given a region of interest
###
freq.table.genome= function(features, motifs,discrete, roi) {
  print (as.character(seqnames(roi[1])))
  ## For each region, make a covariate table. nrows= region length, ncols=number of covariates
  # faster, all the positions have to come from the same chromosome
  selected.features.site1 <- lapply(features,function(f) {as.numeric(f[[as.character(unique(seqnames(roi)))]][ranges(roi)])})
  
  ## discretize all features
  selected.features.site=sapply(1:length(features),function(i) {
    if (class(discrete[[i]])=="numeric") {
      out=ifelse(selected.features.site1[[i]]>discrete[i],1,0)
    }else{
      feat.bins=.bincode(selected.features.site1[[i]], discrete[[i]][[1]])
      out=sapply(feat.bins, function(x) discrete[[i]][[2]][x])
    }})
  selected.features.site=data.frame(selected.features.site)
  colnames(selected.features.site)=names(selected.features.site1)
  
  ## Get a list of sites for each promoter
  sites=tile(roi, width=1)
  sites=unlist(sites)  
  seq=sites+25
  seq=Views(Hsapiens,seq)
  seq=as.data.frame(seq)
  
  #add motifs to selected.features.site
  feat.motif=find.motif(seq, motifs) #a dataframe with 4 columns/motifs
  
  ##add motifs to selected.features.site
  selected.features.site=cbind(selected.features.site,feat.motif) #a dataframe with 19 columns/features
  
  freq.table=count(selected.features.site)
  
  rm(list=c("selected.features.site", "sites", "seq"))
  gc()
  freq.table
}


freq.table.muts= function(features, motifs,discrete, muts) {
  print (as.character(seqnames(muts[1])))
  ## Get a list of sites
  sites=muts
  seq = muts+25
  # strip meta columns from seq
  seq= GRanges(seqnames(seq), ranges(seq))
  seq = Views(Hsapiens, seq)
  seq = as.data.frame(seq)

  ## For each region, make a covariate table. nrows= region length, ncols=number of covariates
  # faster, all the positions have to come from the same chromosome
  selected.features.site1 <- lapply(features,function(f) {as.numeric(f[[as.character(unique(seqnames(sites)))]][ranges(sites)])})
  ## discretize all features
  selected.features.site=sapply(1:length(features),function(i) {
    if (class(discrete[[i]])=="numeric") {
      out=ifelse(selected.features.site1[[i]]>discrete[i],1,0)
    }else{
      feat.bins=.bincode(selected.features.site1[[i]], discrete[[i]][[1]])
      out=sapply(feat.bins, function(x) discrete[[i]][[2]][x])
    }})
  selected.features.site=data.frame(selected.features.site)
  colnames(selected.features.site)=names(selected.features.site1)
  
  #add motifs to selected.features.site
  feat.motif=find.motif(seq, motifs) #a dataframe with 4 columns/motifs
  
  ##add motifs to selected.features.site
  selected.features.site=cbind(selected.features.site,feat.motif, sites$sid) #a dataframe with a column for each feature
  
  freq.table=count(selected.features.site)
  rm(list=c("selected.features.site", "sites", "seq"))
  gc()
  freq.table
}

###
### Get Mutations
###
maf.gastric <- maf.to.granges('/mnt/projects/guoy1/wgs/Gastric_Cancer/data/gastric_RF_nonMSI_prefiltered.MAF')
# find mutations in regions of interest
mut.masked <-maf.gastric[subjectHits(findOverlaps(all.sites.masked,maf.gastric))]

###
### Read in covariate data
###
selected.features = lapply(selected.urls[,2],function(f) {import(as.character(f),as='RleList')})
names(selected.features) = as.character(selected.urls[,1])

## remove larger objects before tabulating
sort( sapply(ls(),function(x){object.size(get(x))/10^6}))
rm(list=c("mappability", "maf.gastric","mask.regions", "nonmappable", "roi.cds", "roi.cds.ext"))
gc(reset=T)

## tabulate covariates for mutations
# split mut.masked by chromosome number
seqlevels(mut.masked)=as.character(unique(seqnames(mut.masked)))
mut.chr=split(mut.masked, seqnames(mut.masked))
mut.chr=mut.chr[1:23]
mut.freq<-mclapply(mut.chr, function(x) freq.table.muts(features=selected.features,motifs=sel.motif, discrete=discrete, x), mc.cores=cores, mc.preschedule = FALSE, mc.silent = FALSE)
#mut.freq
#saveRDS(mut.freq, file=outfile)
# faster alternative to do.call rbind?
mut.freq1=rbindlist(mut.freq) #a data.table and data.frame
mut.freq1=data.frame(mut.freq1, check.names=F) #convert to a dataframe, check.names=F to ensure the '-' in sample names are not replaced with '.'
# sum the frequencies of mutations across all regions of interest
mutfreq.aggregated=aggregate(mut.freq1$freq, by=mut.freq1[,1:(ncol(mut.freq1)-1)], FUN=sum)
rm(list=c("mut.chr", "mut.freq","mut.freq1"))

## divide data into chuncks of n regions each, with all regions in the same chromosome
len=sapply(split(all.sites.masked, seqnames(all.sites.masked)), length)
len=len[len!=0]
len2=sapply(1:length(len), function(i) {sum(len[1:i])})
len2=c(0,len2)
chunk.size=5000
chunks<-lapply(1:length(len), function(j) {lapply(1:ceiling(len[j]/chunk.size), function(i) ((len2[j]+(i-1)*chunk.size+1):(len2[j]+min((i)*chunk.size, len[j]))))})
chunks=do.call(c, chunks)
genome.freq<-mclapply(chunks, function(x) freq.table.genome(features=selected.features,motifs=sel.motif, discrete=discrete, roi=all.sites.masked[x]), mc.cores=cores, mc.preschedule = FALSE, mc.silent = FALSE)
# faster alternative to do.call rbind?
genome.freq=rbindlist(genome.freq) #a data.table and data.frame
genome.freq=data.frame(genome.freq, check.names=F) #convert to a dataframe, check.names=F to ensure the '-' in sample names are not replaced with '.'
# sum the frequencies of mutations across all regions of interest
genome.freq.aggregated=aggregate(genome.freq$freq, by=genome.freq[,1:(ncol(genome.freq)-1)], FUN=sum)
rm(list=c("genome.freq"))

nind=length(unique(mut.masked$sid))
genome.freq.aggregated2=genome.freq.aggregated[rep(1:nrow(genome.freq.aggregated),nind),]
sid= rep(unique(mut.masked$sid),each=nrow(genome.freq.aggregated))
genome.freq.aggregated2=data.frame(genome.freq.aggregated2[,1:(ncol(genome.freq.aggregated2)-1)], sites.sid=sid, tot.count=genome.freq.aggregated2[,ncol(genome.freq.aggregated2)])
aggregated.table=merge(mutfreq.aggregated,genome.freq.aggregated2, by=colnames(genome.freq.aggregated2)[1:(ncol(genome.freq.aggregated2)-1)], all=T)
names(aggregated.table)[(ncol(aggregated.table)-1):ncol(aggregated.table)]=c("mut.count", "tot.count")
aggregated.table[is.na(aggregated.table)]=0
aggregated.table$nonmut.count=aggregated.table$tot.count-aggregated.table$mut.count
saveRDS(aggregated.table, file=outfile)

debug=list(mutfreq.aggregated, genome.freq.aggregated, genome.freq.aggregated2, aggregated.table)
saveRDS(debug, file="/mnt/projects/guoy1/wgs/Gastric_Cancer/data/covariate_genome_freq_table_nonMSI_prefiltered_debug.rds")
