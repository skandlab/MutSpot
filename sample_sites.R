###
### sample all mutated sites and non-mutated sites at regular intervals
### July 19 2016
### qsub -cwd -pe OpenMP 16 -l mem_free=100G,h_rt=48:00:00 -q medium.q -b y R CMD BATCH sample_sites.R 
###
setwd("/mnt/projects/guoy1/wgs/Gastric_Cancer/")
source('~/Gastric_Cancer/scripts/mutrec-lib.R')
library(BSgenome.Hsapiens.UCSC.hg19)
cores=16

chrOrder<-c(paste("chr",1:22,sep=""),"chrX")
seqi = seqinfo(Hsapiens)[seqnames(Hsapiens)[1:23]]

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

# Create Genomic Ranges for the chromosomes 1:22 and chrX
all.sites=ranges[seqnames(ranges)%in% seqnames(seqi)]
sum(as.numeric(width(all.sites))) #2951331845

## mask CDS and ig loci
roi.cds <- bed.to.granges('/mnt/projects/guoy1/wgs/Gastric_Cancer/roi_Ensembl75/Ensembl75.CDS.bed')
roi.cds.ext <- reduce(roi.cds + 5) # extend each region with +/- 5 bases and get all non-overlapping regions
#immunoglobulin loci
ig.loci <- bed.to.granges('/mnt/projects/guoy1/wgs/Gastric_Cancer/data/ig_loci.bed')
ig.loci <- reduce(ig.loci + 10**5) # extend each region with 100kb and get all non-overlapping regions
# combine mask regions
# one of the ig.loci located on sequence chr14 is out-of-bound, trim at this stage as there is no seqinfo associated with ig.loci and roi.cds.ext
mask.regions=reduce(trim(c(ig.loci,roi.cds.ext,nonmappable))) 

###
### Get Mutations
###
maf.gastric <- maf.to.granges('/mnt/projects/guoy1/wgs/Gastric_Cancer/data/gastric_RF.MAF')
maf.gastric=maf.gastric[seqnames(maf.gastric)%in%seqnames(seqi)]
length(maf.gastric) #4189020
npatients = length(levels(maf.gastric$sid))
# extend each mutations by 5 bp on each side
# reduce merges adjacent mutations, use unique instead
maf.gastric <-unique(maf.gastric)
length(maf.gastric) #4130700
maf.gastric= subtract.regions.from.roi(maf.gastric, mask.regions, cores=cores) #3947629

###
### Get a sample of non-mutated bases
###
# target number of sites to samples, take into account of larger masked regions, and that mutated sites tend not be in masked regions
nsites= length(maf.gastric)*1.12
# number of sites to sample per chromosome
nsites.chrom=round(width(all.sites)/sum(as.numeric(width(all.sites)))*nsites)
# sample all sites 
all.sites.samples=lapply(1:23, function(i) {
     pos=sample(start(all.sites)[i]:end(all.sites)[i],nsites.chrom[i])
    gr=GRanges(names(all.sites)[i], IRanges(pos, pos))})
all.sites.samples=GRangesList(all.sites.samples)
all.sites.samples=unlist(all.sites.samples) #4342392
# mask selected sites that are mutated or in nonmapple regions
mask.regions=reduce(c(maf.gastric, mask.regions))
nonmut.sample=subtract.regions.from.roi(all.sites.samples, mask.regions, cores=cores) #3961533
nonmut.sample$mut=0
maf.gastric$mut=1
sampled.sites=sort(c(nonmut.sample, maf.gastric))
saveRDS(sampled.sites, file='/mnt/projects/guoy1/wgs/Gastric_Cancer/data/features/all_sites_downsampled.rds')

###
### Get sampled sites in all regions of interest
###
sampled.sites=readRDS('/mnt/projects/guoy1/wgs/Gastric_Cancer/data/features/all_sites_downsampled.rds')
prom= bed.to.granges("./roi_Ensembl75/roi.prom.mappable.reduced.bed")
prom.sites <-sampled.sites[subjectHits(findOverlaps(prom,sampled.sites))] #101505
sum(prom.sites$mut) #42068
saveRDS(prom.sites, file='/mnt/projects/guoy1/wgs/Gastric_Cancer/data/features/prom_sites_RF_downsampled.rds')

intron= bed.to.granges("./roi_Ensembl75/roi.intron.mappable.reduced.bed")
intron.sites <-sampled.sites[subjectHits(findOverlaps(intron,sampled.sites))] #544323
sum(intron.sites$mut) #226528
saveRDS(intron.sites, file='/mnt/projects/guoy1/wgs/Gastric_Cancer/data/features/intron_sites_RF_downsampled.rds')

utr5= bed.to.granges("./roi_Ensembl75/roi.5UTR.mappable.reduced.bed")
utr5.sites <-sampled.sites[subjectHits(findOverlaps(utr5,sampled.sites))] #21470
sum(utr5.sites$mut) #9062
saveRDS(utr5.sites, file='/mnt/projects/guoy1/wgs/Gastric_Cancer/data/features/UTR5_sites_RF_downsampled.rds')

utr3= bed.to.granges("./roi_Ensembl75/roi.3UTR.mappable.reduced.bed")
utr3.sites <-sampled.sites[subjectHits(findOverlaps(utr3,sampled.sites))] #76369
sum(utr3.sites$mut) #32072
saveRDS(utr3.sites, file='/mnt/projects/guoy1/wgs/Gastric_Cancer/data/features/UTR3_sites_RF_downsampled.rds')

enh= bed.to.granges("/mnt/projects/guoy1/wgs/Gastric_Cancer/roi/roi.enh.bed")
enh=reduce(enh)
enh.sites <-sampled.sites[subjectHits(findOverlaps(enh,sampled.sites))] #27535
sum(enh.sites$mut) #11344
saveRDS(enh.sites, file='/mnt/projects/guoy1/wgs/Gastric_Cancer/data/features/enh_sites_RF_downsampled.rds')