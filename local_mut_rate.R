###
### Calculate local mutation rate for 100kb nonoverlapping bins across the genome
### mask CDS regions, nonmappable regions and ig loci
### July 12 2016
### 

setwd("/mnt/projects/guoy1/wgs/Gastric_Cancer/")
source('~/Gastric_Cancer/scripts/mutrec-lib.R')
library(BSgenome.Hsapiens.UCSC.hg19)

cores=16
## read in mutation data
maf.gastric <- maf.to.granges('/mnt/projects/guoy1/wgs/Gastric_Cancer/data/gastric_RF.MAF')
# will only use chr 1:22 and X
seqi = seqinfo(Hsapiens)[intersect(seqnames(seqinfo(Hsapiens))[1:23],as.character(seqnames(maf.gastric)))]
nind = length(levels(maf.gastric$sid))

## Tile genome into equally sized bins
binsize = 100*1e3
genome.bins <- tileGenome(seqi, tilewidth=binsize, cut.last.tile.in.chrom=TRUE)
nbin = length(genome.bins)
names(genome.bins)=paste("n", seq(1:nbin), sep="")
genome.bins.grl<- split(genome.bins, names(genome.bins))

### mask CDS regions, nonmappable regions and ig loci
mappability=import("/mnt/projects/guoy1/wgs/Gastric_Cancer/data/wgEncodeCrgMapabilityAlign75mer.bigWig")
nonmappable=mappability[mappability$score<1,]
# convert zero-based coordinates to one-based coordinates
nonmappable=shift(nonmappable,1)
nonmappable= reduce(nonmappable)
nonmappable=nonmappable[seqnames(nonmappable) %in% seqnames(seqi)]
seqlevels(nonmappable)=as.character(unique(seqnames(nonmappable)))
#CDS
roi.cds <- bed.to.granges('/mnt/projects/guoy1/wgs/Gastric_Cancer/roi_Ensembl75/Ensembl75.CDS.bed')
roi.cds.ext <- reduce(roi.cds + 5) # extend each region with +/- 5 bases and get all non-overlapping regions
#immunoglobulin loci
ig.loci <- bed.to.granges('/mnt/projects/guoy1/wgs/Gastric_Cancer/data/ig_loci.bed')
ig.loci <- reduce(ig.loci + 10**5) # extend each region with 100kb and get all non-overlapping regions
# combine mask regions
# one of the ig.loci located on sequence chr14 is out-of-bound, trim at this stage as there is no seqinfo associated with ig.loci and roi.cds.ext
mask.regions=reduce(trim(c(ig.loci,roi.cds.ext,nonmappable))) 

genome.bins.grl.masked= subtract.regions.from.roi(genome.bins.grl, mask.regions, cores=cores)

## calculate mutation rate for each bin
genome.bins.length= sum(width(genome.bins.grl.masked))
genome.bins.length=data.frame(name=names(genome.bins.length),length=genome.bins.length)
maf.ovl <- findOverlaps(genome.bins.grl.masked, maf.gastric)
maf.ovl.m = as.matrix(maf.ovl)
genome.bins.mutcount=tapply(maf.ovl.m[,2], names(genome.bins.grl.masked)[maf.ovl.m[,1]], function(s) length(s))
genome.bins.mutcount=data.frame(name=names(genome.bins.mutcount),mutcount=genome.bins.mutcount)
genome.bins.mutrate=merge(genome.bins.length,genome.bins.mutcount,all=T)
genome.bins.mutrate[is.na(genome.bins.mutrate)]=0
genome.bins.mutrate$mut.rate=genome.bins.mutrate$mutcount/genome.bins.mutrate$length/nind

mean.mutrate=sum(genome.bins.mutrate$mutcount)/sum(as.numeric(genome.bins.mutrate$length))/nind
genome.bins.mutrate2=rep(mean.mutrate, length(genome.bins.grl))
names(genome.bins.mutrate2)=names(genome.bins)
genome.bins.mutrate2[as.character(genome.bins.mutrate$name)]=genome.bins.mutrate$mut.rate
values(genome.bins)=data.frame(mutrate=genome.bins.mutrate2[names(genome.bins)])

# adjust start coordinate to be 0-based
df <- data.frame(seqnames=seqnames(genome.bins), starts=as.integer(start(genome.bins)-1), ends=end(genome.bins),mutrate=genome.bins$mutrate)
write.table(df,file="/mnt/projects/guoy1/wgs/Gastric_Cancer/data/mutrate.local.100kb.RF.bed",quote=F, sep="\t", row.names=F, col.names=F)

# convert bed to bigwig format
#/mnt/software/bin/bedGraphToBigWig /mnt/projects/guoy1/wgs/Gastric_Cancer/data/mutrate.local.100kb.RF.bed /mnt/projects/guoy1/wgs/Gastric_Cancer/Mixed_model/RepTimeWig/hg19.chrom.sizes /mnt/projects/guoy1/wgs/Gastric_Cancer/data/mutrate.local.100kb.RF.bigWig

# divide local mutation rate in to n bins and calculate the mean rate in each bin
nbins=10
mutrateBins=quantile(df$mutrate, seq(0,1, 1/nbins))
mutrateBins[1]=mutrateBins[1]-1e-20
mutrateBins[nbins]=mutrateBins[nbins]*1.1
bin.mean=lapply(1:nbins, function(x) mean(df$mutrate[mutrate.bin==x]))
saveRDS(mutrateBins, file='/mnt/projects/guoy1/wgs/Gastric_Cancer/data/localMutRate_bins.rds')
saveRDS(bin.mean, file='/mnt/projects/guoy1/wgs/Gastric_Cancer/data/localMutRate_bin_mean.rds')
