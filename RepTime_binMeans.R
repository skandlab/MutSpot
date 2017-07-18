###
### Bin replication timing according to their quantiles and calculate mean values for each bin
### July 7 2016
### qsub -cwd -pe OpenMP 4 -l mem_free=16G,h_rt=2:00:00 -q medium.q -b y R CMD BATCH RepTime_binMeans.R
### 

setwd("/mnt/projects/guoy1/wgs/Gastric_Cancer/")
source('~/Gastric_Cancer/scripts/mutrec-lib.R')
library(BSgenome.Hsapiens.UCSC.hg19)

seqi = seqinfo(Hsapiens)[seqnames(Hsapiens)[1:23]]
seqnames=seqnames(seqinfo(Hsapiens))[1:23]

# integer overflow on attempt to concatenate reptime rlelist; as the reptime data is in 1000bp bins, convert rlelist to granges of 1000bp width
rep.time = import("/mnt/projects/guoy1/wgs/Gastric_Cancer/Mixed_model/RepTimeWig/wgEncodeUwRepliSeqWaveSignalMean.bigWig",as='RleList')
binsize = 1*1e3
trim=runLength(rep.time)
tail = sapply (trim, function(y)  tail(y, n=1))
head = sapply (trim, function(y)  y[[1]])

tail=tail[seqnames]
seqlengths(seqi)=seqlengths(seqi)-tail
genome.bins <- tileGenome(seqi, tilewidth=binsize, cut.last.tile.in.chrom=TRUE)
genome.bins.ls=split(genome.bins, seqnames(genome.bins))

head=head[names(genome.bins.ls)]
genome.bins.ls= mapply ( function(x,y) trim(shift(x, y)), genome.bins.ls, head)
genome.bins.ls= GRangesList(genome.bins.ls)
genome.bins2= unlist(genome.bins.ls)
invalid.range=which(width(genome.bins2)==0)
genome.bins2=genome.bins2[-invalid.range]

# number of bins to divide replication timing into
nbins=8
rep.time=bigwig.summarize.bins("/mnt/projects/guoy1/wgs/Gastric_Cancer/Mixed_model/RepTimeWig/wgEncodeUwRepliSeqWaveSignalMean.bigWig",genome.bins2)
RepBins=quantile(rep.time, seq(0,1, 1/nbins))
# take the floor for the minimun value to make sure everything falls under the range when bining
RepBins[1]=floor(RepBins[1])

#RepBins=c(floor(-3.229303), 19.310651, 27.845977, 37.270827, 45.996883, 53.899962, 61.491663, 69.232374, 83.230598)
rep.time.bin=.bincode(rep.time, RepBins)
bin.mean=lapply(1:nbins, function(x) mean(rep.time[rep.time.bin==x]))
saveRDS(RepBins, file='/mnt/projects/guoy1/wgs/Gastric_Cancer/data/repTime_bins.rds')
saveRDS(bin.mean, file='/mnt/projects/guoy1/wgs/Gastric_Cancer/data/repTime_bin_mean.rds')