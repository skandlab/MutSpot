#' Dilution analysis.
#' 
#' @param snv.mutations.file SNV mutations MAF file.
#' @param mask.regions.file Regions to mask in genome, for example, non-mappable regions/immunoglobin loci/CDS regions RDS file, default file = mask_regions.RDS.
#' @param all.sites.file All sites in whole genome RDS file, default file = all_sites.RDS.
#' @param cores Number of cores, default = 1.
#' @param cutoff.nucleotide Frequency cutoff/threshold to determine nucleotide contexts used in prediction model, ranges from 0.5 to 1, default = 0.90.
#' @param cutoff.nucleotide.new Updated frequency cutoff/threshold to determine nucleotide contexts used in prediction model, ranges from 0.5 to 1.
#' @param cutoff.features Frequency cutoff/threshold to determine epigenetic features used in prediction model, ranges from 0.5 to 1, default = 0.75.
#' @param cutoff.features.new.snv Updated frequency cutoff/threshold to determine SNV epigenetic features used in prediction model, ranges from 0.5 to 1.
#' @param genomic.features.snv Text file containing URLs of potential continuous and discrete SNV epigenetic features to select from, default = NULL.
#' @param genomic.features.indel Text file containing URLs of potential continuous and discrete indel epigenetic features to select from, default = NULL.
#' @param genomic.features Text file containing URLs of potential continuous and discrete SNV and indel epigenetic features to select from, default = NULL.
#' @param genomic.features.fixed.snv Text file containing URLs of fixed continuous and discrete SNV epigenetic features, default = NULL.
#' @param genomic.features.fixed.indel Text file containing URLs of fixed continuous and discrete indel epigenetic features, default = NULL.
#' @param genomic.features.fixed Text file containing URLs of fixed continuous and discrete SNV and indel epigenetic features, default = NULL.
#' @param feature.dir Directory containing binned feature bed files.
#' @param mutCovariate.table.snv.file RDS file, SNV covariates sparse matrix.
#' @param mutCovariate.count.snv.file RDS file, SNV response matrix.
#' @param continuous.features.selected.snv.url.file Text file containing URLs of SNV continuous features selected for model.
#' @param discrete.features.selected.snv.url.file Text file containing URLs of SNV discrete features selected for model.
#' @param nucleotide.selected.file Nucleotide context selected for model RDS file.
#' @param sample.specific.features.url.file Text file containing sample specific SNV features, default = NULL.
#' @param output.dir Save plots in given output directory.
#' @return McFadden's R2 for different dilution size.
#' @export

dilution.test=function(snv.mutations.file,mask.regions.file = system.file("extdata", "mask_regions.RDS", package = "MutSpot"),
                       all.sites.file = system.file("extdata", "all_sites.RDS", package = "MutSpot"),
                       cores=1,cutoff.nucleotide=0.9,cutoff.nucleotide.new=NULL,cutoff.features = 0.75,cutoff.features.new.snv=NULL,
                       genomic.features.snv=NULL,genomic.features.indel=NULL,genomic.features=NULL,genomic.features.fixed.snv=NULL,
                       genomic.features.fixed.indel=NULL,genomic.features.fixed=NULL,feature.dir,
                       mutCovariate.table.snv.file,
                       mutCovariate.count.snv.file, 
                       continuous.features.selected.snv.url.file, 
                       discrete.features.selected.snv.url.file, 
                       nucleotide.selected.file, 
                       sample.specific.features.url.file=NULL,output.dir) {
  
## Sample sites (20%, 40%, 60%, 80%)
dilution.size= c(0.2,0.4,0.6,0.8)

# Chr1-ChrX
chrOrder <- c(paste("chr", 1:22, sep=""), "chrX")
seqi = GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[GenomeInfoDb::seqnames(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[1:23]]

# Define masked region i.e. CDS, immunoglobulin loci and nonmappable
mask.regions = readRDS(mask.regions.file)
mask.regions = mask.regions[as.character(GenomeInfoDb::seqnames(mask.regions)) %in% as.character(GenomeInfoDb::seqnames(seqi))]

# Define all sites in whole genome
all.sites = readRDS(all.sites.file)
all.sites = all.sites[as.character(GenomeInfoDb::seqnames(all.sites)) %in% as.character(GenomeInfoDb::seqnames(seqi))]

# Define SNV mutations
maf.snv.mutations <- maf.to.granges(snv.mutations.file)
maf.snv.mutations = maf.snv.mutations[as.character(GenomeInfoDb::seqnames(maf.snv.mutations)) %in% as.character(GenomeInfoDb::seqnames(seqi))]
maf.snv.mutations2 = maf.snv.mutations

for (size in dilution.size) {
  
  print(paste("Sample ",size*100, "% sites",sep=""))
maf.snv.mutations=maf.snv.mutations2
max.sites = length(maf.snv.mutations2)*size*2*1.12
min.sites = 4000 * 1.12

#if (length(maf.snv.mutations) > max.sites / 2) {
  
  downsample.snv = TRUE
  print(paste("Downsample SNV mutations as number of SNV mutations exceeded ", max.sites, sep = ""))
  
if (length(maf.snv.mutations) < min.sites / 2) {
  
  ratio.snv = ceiling((min.sites - length(maf.snv.mutations)) / length(maf.snv.mutations))
  print(paste("Ratio of number of mutated sites to non-mutated sites for SNV is 1:", ratio.snv, sep = ""))
  
} else {
  
  ratio.snv = 1
  
}

nsites.snv = max.sites / 2

t = GenomicRanges::split(maf.snv.mutations, GenomeInfoDb::seqnames(maf.snv.mutations))
nsites.snv.chrom = round(unlist(lapply(t, FUN=function(x) sum(as.numeric(GenomicRanges::width(x))))) / sum(unlist(lapply(t, FUN = function(x) sum(as.numeric(GenomicRanges::width(x)))))) * nsites.snv)
seed.rand.snv = seq(1:length(t)) * 4 

# Downsample sites
downsampled.snv.sites = parallel::mclapply(1:length(t), function(i) {
  
  pop = IRanges::tile(t[[i]], width = 1)
  pop = BiocGenerics::unlist(pop)
  set.seed(seed.rand.snv[i])
  pos = sample(GenomicRanges::start(pop), nsites.snv.chrom[i])
  
  if (length(pos) > 0) {
    
    gr = GenomicRanges::GRanges(unique(as.character(GenomeInfoDb::seqnames(t[[i]]))), IRanges::IRanges(pos, pos))
    return(gr)
    
  } else {
    
    return(NULL)
    
  }
  
}, mc.cores = cores)

downsampled.snv.sites[sapply(downsampled.snv.sites, is.null)] <- NULL
downsampled.snv.sites = suppressWarnings(do.call(c, downsampled.snv.sites))
maf.snv.mutations = downsampled.snv.sites

npatients.snv = length(unique(maf.snv.mutations$sid))
maf.snv.mutations <- unique(maf.snv.mutations)

# Remove SNV mutations in masked region
maf.snv.mutations = subtract.regions.from.roi(maf.snv.mutations, mask.regions, cores=cores)

# Target number of sites to sample, take into account of larger masked regions, and that mutated sites tend not be in masked regions
nsites.snv = length(maf.snv.mutations)*(ratio.snv + ratio.snv * 0.12)

if (nsites.snv < c(10000 * 1.12)) {
  
  nsites.snv = 10000 * 1.12
  
}

all.sites.snv = all.sites

# Number of sites to sample per chromosome
nsites.snv.chrom = round(GenomicRanges::width(all.sites.snv) / sum(as.numeric(GenomicRanges::width(all.sites.snv))) * nsites.snv)
seed.rand.snv = seq(1:length(all.sites.snv)) * 4

# Sample sites
all.sites.snv.samples = parallel::mclapply(1:length(all.sites.snv), function(i) {
  
  set.seed(seed.rand.snv[i])
  pos = sample(GenomicRanges::start(all.sites.snv)[i]:GenomicRanges::end(all.sites.snv)[i], nsites.snv.chrom[i])
  if (length(pos) > 0) {
    
    gr = GenomicRanges::GRanges(as.character(GenomeInfoDb::seqnames(all.sites.snv)[i]), IRanges::IRanges(pos, pos))
    return (gr)
    
  } else { 
    
    return(NULL)
    
  }
  
}, mc.cores = cores)

all.sites.snv.samples[sapply(all.sites.snv.samples, is.null)] <- NULL
if (length(all.sites.snv.samples) == 0) {
  
  filtered.snv.mutations = NULL
  sampled.snv.sites = NULL
  
} else {
  
  # all.sites.snv.samples = suppressWarnings(do.call(getMethod(c, "GenomicRanges"), all.sites.snv.samples))
  all.sites.snv.samples = suppressWarnings(do.call(c, all.sites.snv.samples))
  
  # Mask selected sites that are mutated or in masked region
  GenomicRanges::mcols(maf.snv.mutations2) = NULL
  mask.snv.regions = GenomicRanges::reduce(c(maf.snv.mutations2, mask.regions))
  nonmut.snv.sample = subtract.regions.from.roi(all.sites.snv.samples, mask.snv.regions, cores = cores)
  
  if (length(nonmut.snv.sample) != 0) {
    
    nonmut.snv.sample$mut = 0
    maf.snv.mutations$mut = 1
    sampled.snv.sites = sort(c(nonmut.snv.sample, maf.snv.mutations))
    
  } else {
    
    sampled.snv.sites = NULL
    filtered.snv.mutations = NULL
    
  }
  
  # print(table(sampled.snv.sites$mut))
saveRDS(sampled.snv.sites,file=paste(output.dir,"sampled_",size,".RDS",sep=""))

}
}

## Nucleotide selection
if (!is.null(cutoff.nucleotide.new)) {
  n.cutoff=cutoff.nucleotide.new
} else {
  
  n.cutoff=cutoff.nucleotide
}

for (size in dilution.size) {
  
  print(paste("Nucleotide selection for ", size*100,"% sites",sep=""))
  
  nucleotide.test = nucleotide.selection(sampled.sites.snv.file=paste("sampled_",size,".RDS",sep=""), indel.mutations.file=NULL, cutoff = n.cutoff, cores = cores) #1min
  # print(nucleotide.test[[3]])
  saveRDS(nucleotide.test,file=paste(output.dir,"nucleotide_",size,".RDS",sep=""))
  
}

## Feature selection
if (!is.null(cutoff.features.new.snv)) {
  f.cutoff=cutoff.features.new.snv
} else {
  
  f.cutoff=cutoff.features
}

for (size in dilution.size) {
  
  print(paste("Epigenomic features selection for ", size*100,"% sites",sep=""))
  
  epig.test=epigenetic.selection(sampled.sites.snv.file=paste("sampled_",size,".RDS",sep=""), sampled.sites.indel.file=NULL, genomic.features.snv = genomic.features.snv, genomic.features.indel = NULL, genomic.features = genomic.features, genomic.features.fixed.snv = genomic.features.fixed.snv, genomic.features.fixed.indel = NULL, genomic.features.fixed = genomic.features.fixed, cores = 4, cutoff = 0.75, feature.dir=feature.dir) #4min
#    print(epig[[3]])
#   print(epig[[4]])
  saveRDS(epig.test,file=paste(output.dir,"epig_",size,".RDS",sep=""))
  
}

## Model Fit and McFadden R2
mutfreq.aggregated = readRDS(file = mutCovariate.table.snv.file)
# Define 2-column response matrix
mutfreq.aggregated2 = readRDS(file = mutCovariate.count.snv.file)

mutfreq.aggregated = as.data.frame(as.matrix(mutfreq.aggregated))
mutfreq.aggregated = cbind(mutfreq.aggregated2, mutfreq.aggregated)

if (!is.null(continuous.features.selected.snv.url.file)) {
  
  selected.continuous.urls <- read.delim(continuous.features.selected.snv.url.file, header = FALSE, stringsAsFactors = FALSE)
  selected.continuous.urls = selected.continuous.urls[ ,1]
  
} else {
  
  selected.continuous.urls = NULL
  
}

if (!is.null(sample.specific.features.url.file)) {
  
  sample.specific.urls <- read.delim(sample.specific.features.url.file, stringsAsFactors = FALSE)
  mutfreq.aggregated=mutfreq.aggregated[,-which(colnames(mutfreq.aggregated)%in%colnames(sample.specific.urls))]
  
}  
  
continuous.sample.specific = NULL

for(i in colnames(mutfreq.aggregated)[which(!colnames(mutfreq.aggregated) %in% c("mut.count", "nonmut.count", "ind.mut.count", selected.continuous.urls, continuous.sample.specific))]) {
  
  mutfreq.aggregated[ ,i] = as.character(mutfreq.aggregated[ ,i])
  
}
mutfreqz=mutfreq.aggregated

mcfadden.r2=as.list(numeric(5))
names(mcfadden.r2)=c(paste("subsample",dilution.size,sep=""),"subsample1")

for (size in dilution.size) {
  
  print(paste("Model fit for ",size*100,"% sites",sep=""))
  n.feat=readRDS(paste(output.dir,"nucleotide_",size,".RDS",sep=""))[[3]]
  n.feat=gsub("one","oneMer",n.feat)
  n.feat=gsub("three","threeMer",n.feat)
  n.feat=gsub("five","fiveMer",n.feat)
  n.feat=gsub("threeMer.left","threeLeft",n.feat)
  n.feat=gsub("threeMer.right","threeRight",n.feat)
  n.feat=gsub("fiveMer.right","fiveRight",n.feat)
  n.feat=gsub("fiveMer.left","fiveLeft",n.feat)
  
  e.feat=readRDS(paste(output.dir,"epig_",size,".RDS",sep=""))[3:4]
  e.feat=do.call(rbind,e.feat)
  e.feat=e.feat[,1]
  
  cols.keep=c("ind.mut.count","mut.count","nonmut.count",as.character(n.feat),as.character(e.feat))
  mutfreq.aggregated=mutfreqz[,cols.keep]
  LRmodel <- stats::glm(cbind(mut.count, nonmut.count) ~ ., family = binomial(logit), data = mutfreq.aggregated, x = F, y = F)
  saveRDS(LRmodel,file=paste(output.dir,"model",size,".RDS",sep=""))
  mut.count=mutfreq.aggregated$mut.count
  nonmut.count=mutfreq.aggregated$nonmut.count
  mcfadden.r2[[paste("subsample",size,sep="")]]=pscl::pR2(LRmodel)["McFadden"]
  
}

print(paste("Model fit for 100% sites",sep=""))
mutfreq.aggregated=mutfreqz
LRmodel <- stats::glm(cbind(mut.count, nonmut.count) ~ ., family = binomial(logit), data = mutfreq.aggregated, x = F, y = F)
saveRDS(LRmodel,file=paste(output.dir,"model",size,".RDS",sep=""))
mut.count=mutfreq.aggregated$mut.count
nonmut.count=mutfreq.aggregated$nonmut.count
mcfadden.r2[["subsample1"]]=pscl::pR2(LRmodel)["McFadden"]

## Figure outputs

## McFadden's R2
print("Plot McFadden's R2")
subsample=c(20,40,60,80,100)
rsq=unlist(mcfadden.r2)
ymin=round(min(rsq)-0.001,digits=3)
ymax=round(max(rsq)+0.001,digits=3)
# Plot runtime axis
pdf(paste(output.dir,"dilution_R2.pdf",sep=""))
plot(subsample, rsq, pch = 16, axes = FALSE, ylim = c(ymin,ymax), xlab = "", ylab = "", type = "b", col = "black")
axis(2, at=pretty(range(ymin,ymax), 5), col = "black", col.axis = "black",cex.axis=0.8,las=0)
mtext("McFadden's R2", side = 2, line = 2.5)
box()
par(new = TRUE)
# Plot sample size axis
mtext("Percentage of sampled sites / %", side = 1, col = "black", line = 2.5)  
axis(1, at=subsample, col = "black", col.axis = "black",cex.axis=0.6,las=2)
dev.off()

## No. of features selected
print("Plot number of features selected")
# nucleotide selection
n.feats=parallel::mclapply(dilution.size,FUN=function(size) {
  
  n.feat=readRDS(paste(output.dir,"nucleotide_",size,".RDS",sep=""))[[3]]
  return(length(n.feat))
},mc.cores=cores)
names(n.feats)=names(mcfadden.r2)[1:4]
n.feats=c(n.feats,length(readRDS(nucleotide.selected.file)))
names(n.feats)[5]="subsample1"


e.feats=parallel::mclapply(dilution.size,FUN=function(size) {
  
  e.feat=readRDS(paste(output.dir,"epig_",size,".RDS",sep=""))[3:4]
  e.feat=do.call(rbind,e.feat)
  e.feat=e.feat[,1]
  return(length(e.feat))
},mc.cores=cores)
names(e.feats)=names(n.feats)[1:4]

e.feats=c(e.feats,ncol(mutfreq.aggregated)-3-n.feats[[5]])
names(e.feats)[5]="subsample1"

nucleo=unlist(n.feats)
epig=unlist(e.feats)
ymax=max(c(nucleo,epig))

pdf(paste(output.dir,"number_features_selected_line.pdf",sep=""))
par(mar = c(5, 5, 2, 5))
# Plot runtime axis
plot(subsample, nucleo, pch = 16, axes = FALSE, ylim = c(0, ymax), xlab = "", ylab = "", type = "b", col = "black")
axis(2, pretty(range(seq(0:ymax)), 5), col = "black", las = 1)
mtext("No. of features selected", side = 2, line = 2.5)
box()

par(new = TRUE)

# Plot memory usage axis
plot(subsample, epig, pch = 15, xlab = "", ylab = "", ylim = c(0, 8), axes = FALSE, type = "b", col = "red", lty = 2)

# Plot sample size axis
mtext("Percentage of sampled sites / %", side = 1, col = "black", line = 2.5)  
axis(1, at=subsample, col = "black", col.axis = "black",cex.axis=0.9,las=2)

# Add Legend
legend("bottomright", legend = c("Nucleotide context", "Epigenomic feature"), text.col = c("black", "red"), lty = c(1, 2), pch = c(16, 15), col = c("black", "red"), cex = 0.8)
dev.off()

## Selected frequency heatmap
print("Plot stability frequency of features selected")
feats=colnames(mutfreq.aggregated)
feats=feats[!feats %in% c("mut.count","nonmut.count","ind.mut.count")]
epig=feats[which(!grepl("Mer",feats) & !grepl("Left",feats) & !grepl("Right",feats))]
nucleo=feats[-which(!grepl("Mer",feats) & !grepl("Left",feats) & !grepl("Right",feats))]

df=NULL
for(size in dilution.size) {
  
nucleotide=readRDS(paste(output.dir,"nucleotide_",size,".RDS",sep=""))[[1]]  
nucleotide[,1]=as.character(nucleotide[,1])
nucleotide[,1]=gsub("one","oneMer",nucleotide[,1])
nucleotide[,1]=gsub("three","threeMer",nucleotide[,1])
nucleotide[,1]=gsub("five","fiveMer",nucleotide[,1])
nucleotide[,1]=gsub("threeMer.left","threeLeft",nucleotide[,1])
nucleotide[,1]=gsub("threeMer.right","threeRight",nucleotide[,1])
nucleotide[,1]=gsub("fiveMer.right","fiveRight",nucleotide[,1])
nucleotide[,1]=gsub("fiveMer.left","fiveLeft",nucleotide[,1])
nucleotide=nucleotide[which(nucleotide[,1] %in% nucleo),]
rownames(nucleotide)=NULL
nucleotide$size=as.character(size*100)
colnames(nucleotide)[2]="freq"
 df=rbind(df,nucleotide)
}

nucleotide=readRDS(paste(output.dir,"nucleotide_stabs_freq.RDS",sep=""))[[1]]  
nucleotide[,1]=as.character(nucleotide[,1])
nucleotide[,1]=gsub("one","oneMer",nucleotide[,1])
nucleotide[,1]=gsub("three","threeMer",nucleotide[,1])
nucleotide[,1]=gsub("five","fiveMer",nucleotide[,1])
nucleotide[,1]=gsub("threeMer.left","threeLeft",nucleotide[,1])
nucleotide[,1]=gsub("threeMer.right","threeRight",nucleotide[,1])
nucleotide[,1]=gsub("fiveMer.right","fiveRight",nucleotide[,1])
nucleotide[,1]=gsub("fiveMer.left","fiveLeft",nucleotide[,1])
nucleotide=nucleotide[which(nucleotide[,1] %in% nucleo),]
rownames(nucleotide)=NULL
nucleotide$size="100"
colnames(nucleotide)[2]="freq"
df=rbind(df,nucleotide)

df2=df
df2=df2[order(df2$freq,df2$size),]

df=NULL
for(size in dilution.size) {
  
  epigenomic=readRDS(paste(output.dir,"epig_",size,".RDS",sep=""))[[1]]
  epigenomic=epigenomic[which(epigenomic[,1] %in% epig),]
  rownames(epigenomic)=NULL
  epigenomic$size=as.character(size*100)
  colnames(epigenomic)[2]="freq"
  df=rbind(df,epigenomic)
}

epigenomic=readRDS(paste(output.dir,"features_stabs_snv.RDS",sep=""))[[1]]
epigenomic=epigenomic[which(epigenomic[,1] %in% epig),]
rownames(epigenomic)=NULL
epigenomic$size="100"
colnames(epigenomic)[2]="freq"
df=rbind(df,epigenomic)
df=df[order(df$freq,df$size),]

df=rbind(df,df2)
df$feature=factor(df$feature,levels = unique(df$feature))
df$size=factor(df$size,levels=subsample)

pdf(paste(output.dir,"stability_freq_heatmap.pdf",sep=""))
print(ggplot2::ggplot(df,aes(x=size,y=feature))+geom_tile(aes(fill=freq),colour="grey")+xlab("Percentage of sampled sites / %")+
  ylab(NULL)+
  geom_text(aes(label = freq)) +
  scale_fill_gradient(low = "white",high = "steelblue",limits=c(0,1))+theme_minimal() +theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                             panel.grid.minor = element_blank(),axis.line=element_line()))
dev.off()



 }

