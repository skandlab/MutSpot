#' Mutation hotspot recurrence prediction for SNV.
#' 
#' @param mask.regions.file Regions to mask in genome, for example, non-mappable regions/immunoglobin loci/CDS regions RDS file, default file = mask_regions.RDS.
#' @param nucleotide.selected.file Nucleotide context selected for model RDS file.
#' @param continuous.features.selected.snv.url.file Text file containing URLs of SNV continuous features selected for model.
#' @param discrete.features.selected.snv.url.file Text file containing URLs of SNV discrete features selected for model.
#' @param sample.specific.features.url.file Text file containing URLs of sample specific features, default = NULL.
#' @param snv.mutations.file SNV mutations found in region of interest MAF file.
#' @param snv.mutations.file2 SNV mutations MAF file.
#' @param region.of.interest Region of interest bed file, default = NULL.
#' @param cores Number of cores, default = 1.
#' @param snv.model.file SNV model.
#' @param min.count Minimum number of mutated samples in each hotspot, default = 2.
#' @param genome.size Total number of hotspots to run analysis on, default = 2533374732.
#' @param hotspots To run hotspot analysis or region-based analysis, default = TRUE.
#' @param merge.hotspots To plot overlapping hotspots as 1 hotspot or individual hotspots, default = TRUE.
#' @return Dataframe containing predicted hotspots significance with hotspots information for SNV.
#' @export

mutPredict.snv = function(mask.regions.file = system.file("extdata", "mask_regions.RDS", package = "MutSpot"), nucleotide.selected.file, continuous.features.selected.snv.url.file, discrete.features.selected.snv.url.file, sample.specific.features.url.file = NULL, snv.mutations.file, snv.mutations.file2, region.of.interest, cores = 1, snv.model.file, min.count = 2, genome.size = 2533374732, hotspots = TRUE, merge.hotspots = TRUE) {
  
# Chr1-X
chrOrder <- c(paste("chr", 1:22, sep = ""), "chrX")
seqi = GenomicRanges::seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[GenomicRanges::seqnames(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[1:23]]
seqnames = GenomicRanges::seqnames(GenomicRanges::seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens))[1:23]

# Define masked region i.e. CDS, immunoglobulin loci and nonmappable
mask.regions = readRDS(mask.regions.file)
mask.regions = mask.regions[as.character(GenomicRanges::seqnames(mask.regions)) %in% seqnames]

# Define SNV mutations in region of interest
maf.snv <- maf.to.granges(snv.mutations.file)
maf.snv = maf.snv[as.character(GenomicRanges::seqnames(maf.snv)) %in% seqnames]
GenomeInfoDb::seqlevels(maf.snv) = as.character(unique(GenomicRanges::seqnames(maf.snv)))

# Define SNV sample mutation count based on full SNV mutations file
maf.snv2 <- maf.to.granges(snv.mutations.file2)
maf.snv2 = maf.snv2[as.character(GenomicRanges::seqnames(maf.snv2)) %in% seqnames]
maf.ind.snv = GenomicRanges::split(maf.snv2, maf.snv2$sid)
ind.mut.count.snv = sapply(maf.ind.snv, length)

# Define sample-specific features e.g. CIN index, COSMIC signatures
if (!is.null(sample.specific.features.url.file)) {
  
  sample.specific.features = read.delim(sample.specific.features.url.file,stringsAsFactors = FALSE)
  if (nrow(sample.specific.features)!=0) {
  rownames(sample.specific.features)=as.character(sample.specific.features$SampleID)
  sample.specific.features=sample.specific.features[which(sample.specific.features$SampleID %in% names(ind.mut.count.snv)),]
  sample.specific.features$sample.count=ind.mut.count.snv[rownames(sample.specific.features)]
  sample.specific.features=sample.specific.features[,-which(colnames(sample.specific.features)=="SampleID")]
  
  sample.specific.features2=parallel::mclapply(1:ncol(sample.specific.features), FUN=function(x) {
    
    print(colnames(sample.specific.features)[x])
    if(class(sample.specific.features[,x])=="character") {
      
      t=factor(sample.specific.features[,x])
      t=model.matrix(~t)[,-1]
      colnames(t)=substr(colnames(t),2,nchar(colnames(t)))
      colnames(t)=paste(colnames(sample.specific.features)[x],colnames(t),sep="")
      rownames(t)=rownames(sample.specific.features)
      
    } else {
      
      t=as.data.frame(sample.specific.features[,x])
      colnames(t)=colnames(sample.specific.features)[x]
      rownames(t)=rownames(sample.specific.features)
      
    }
    return(t)
    
  },mc.cores=cores)
  
  sample.specific.features=do.call(cbind,sample.specific.features2)
  
  sample.specific.urls <- read.delim(sample.specific.features.url.file, stringsAsFactors = FALSE)
  continuous.sample.specific=NULL
  for (j in 1:ncol(sample.specific.urls)) {
    
    if (class(sample.specific.urls[,j]) != "character"){
      
      continuous.sample.specific=c(continuous.sample.specific,colnames(sample.specific.urls)[j])
      
    }
    
  }
  
} else {
  sample.specific.features=as.data.frame(ind.mut.count.snv)
  colnames(sample.specific.features)="sample.count"
  continuous.sample.specific = "sample.count"
}
  
} else {
  
  sample.specific.features=as.data.frame(ind.mut.count.snv)
  colnames(sample.specific.features)="sample.count"
  continuous.sample.specific = "sample.count"
  
}

# Remove masked regions from SNV mutations
maf.masked.snv <-maf.snv[-S4Vectors::subjectHits(IRanges::findOverlaps(mask.regions, maf.snv))]
dupl.snv = duplicated(maf.masked.snv)
maf.uniq.snv = maf.masked.snv[!dupl.snv, ]

if (length(maf.uniq.snv) == 0) {
  
  mut.rec.hotspot <- data.frame(chrom=character(),
                   start=integer(),
                   end=integer(),
                   pval=numeric(),
                   length=integer(),
                   p.bg=numeric(),
                   k=integer(),
                   fdr=numeric(),
                   stringsAsFactors=FALSE)
  
  mut.rec.hotspot2=NULL
  
} else {
  
# Extend each mutation with +/- 10 bases to define hotspots
mut.regions = maf.uniq.snv + 10
names(mut.regions) = paste("Mutation", c(1:length(mut.regions)), sep = "")

rm(list = c("maf.uniq.snv", "mask.regions"))
gc()

k=glmnet::glmnet(x=matrix(sample(c(1,0),200,replace=TRUE),ncol=2),y=matrix(sample(c(1,0),100,replace=TRUE),ncol=1))
rm(k)

load(file = snv.model.file)
LRmodel.snv = LRmodel

# Define selected nucleotide contexts
if (!is.null(nucleotide.selected.file)) {
  
nucleotide.selected = readRDS(nucleotide.selected.file)
if (!is.null(nucleotide.selected)){
nucleotide.selected = as.data.frame(nucleotide.selected)
colnames(nucleotide.selected) = "sequence"
nucleotide.selected$type = "type"
nucleotide.selected$type = ifelse(grepl("one", nucleotide.selected$sequence), "oneMer", nucleotide.selected$type)
nucleotide.selected$type = ifelse(grepl("three", nucleotide.selected$sequence), "threeMer", nucleotide.selected$type)
nucleotide.selected$type = ifelse(grepl("three.right", nucleotide.selected$sequence), "threeRight", nucleotide.selected$type)
nucleotide.selected$type = ifelse(grepl("three.left", nucleotide.selected$sequence), "threeLeft", nucleotide.selected$type)
nucleotide.selected$type = ifelse(grepl("five", nucleotide.selected$sequence), "fiveMer", nucleotide.selected$type)
nucleotide.selected$type = ifelse(grepl("five.right", nucleotide.selected$sequence), "fiveRight", nucleotide.selected$type)
nucleotide.selected$type = ifelse(grepl("five.left", nucleotide.selected$sequence), "fiveLeft", nucleotide.selected$type)
nucleotide.selected$sequence = gsub("one", "", nucleotide.selected$sequence)
nucleotide.selected$sequence = gsub("three", "", nucleotide.selected$sequence)
nucleotide.selected$sequence = gsub("five", "", nucleotide.selected$sequence)
nucleotide.selected$sequence = gsub("*.*right", "", nucleotide.selected$sequence)
nucleotide.selected$sequence = gsub("*.*left", "", nucleotide.selected$sequence)

sel.motif = list()

if ("oneMer" %in% nucleotide.selected$type){ sel.motif$oneMer = nucleotide.selected[which(nucleotide.selected$type == "oneMer"), "sequence"] }
if ("threeMer" %in% nucleotide.selected$type){ sel.motif$threeMer = nucleotide.selected[which(nucleotide.selected$type == "threeMer"), "sequence"] }
if ("threeRight" %in% nucleotide.selected$type){ sel.motif$threeRight = nucleotide.selected[which(nucleotide.selected$type == "threeRight"),"sequence"] }
if ("threeLeft" %in% nucleotide.selected$type){ sel.motif$threeLeft = nucleotide.selected[which(nucleotide.selected$type == "threeLeft"), "sequence"] }
if ("fiveMer" %in% nucleotide.selected$type){ sel.motif$fiveMer = nucleotide.selected[which(nucleotide.selected$type == "fiveMer"), "sequence"] }
if ("fiveRight" %in% nucleotide.selected$type){ sel.motif$fiveRight = nucleotide.selected[which(nucleotide.selected$type == "fiveRight"), "sequence"] }
if ("fiveLeft" %in% nucleotide.selected$type){ sel.motif$fiveLeft = nucleotide.selected[which(nucleotide.selected$type == "fiveLeft"), "sequence"] }
} else {
  nucleotide.selected = sel.motif = NULL
  
}
} else {
  
  nucleotide.selected = sel.motif = NULL
  
  }

# Read feature file paths
if (!is.null(continuous.features.selected.snv.url.file)) {
  
selected.continuous.urls.snv <- read.delim(continuous.features.selected.snv.url.file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)

if (nrow(selected.continuous.urls.snv)!=0) {
continuous.selected.features.snv = parallel::mclapply(selected.continuous.urls.snv[ ,2], function(f) {
  
  print(f)
  df = read.delim(as.character(f), stringsAsFactors = FALSE, header = FALSE)
  with(df, GenomicRanges::GRanges(V1, IRanges::IRanges(V2, V3),score = V4))
  
  }, mc.cores = cores)
names(continuous.selected.features.snv) = as.character(selected.continuous.urls.snv[ ,1])
} else {
  continuous.selected.features.snv = NULL
  
}
} else {
  
  continuous.selected.features.snv = NULL

}

if (!is.null(discrete.features.selected.snv.url.file)) {
  
selected.discrete.urls.snv <- read.delim(discrete.features.selected.snv.url.file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
if (nrow(selected.discrete.urls.snv)!=0) {

discrete.selected.features.snv = parallel::mclapply(selected.discrete.urls.snv[ ,2], function(f) {
  
  print(f)
  df = read.delim(as.character(f), stringsAsFactors = FALSE, header = FALSE)
  with(df, GenomicRanges::GRanges(V1, IRanges::IRanges(V2, V3)))
  
  }, mc.cores = cores)
names(discrete.selected.features.snv) = as.character(selected.discrete.urls.snv[ ,1])
} else {
  discrete.selected.features.snv = NULL
  
}
} else {
  
  discrete.selected.features.snv = NULL
  
  }

continuous.sample.specific=c(continuous.sample.specific,names(continuous.selected.features.snv))

# Run hotspot recurrence analysis
if (is.null(region.of.interest)) {
  
  mut.rec <- mutPredict.snv.run.lr(roi = mut.regions, maf.snv = maf.snv, maf.snv2 = maf.snv2, model.snv = LRmodel.snv, continuous.features.snv = continuous.selected.features.snv, discrete.features.snv = discrete.selected.features.snv, motifs = sel.motif, nucleotide.selected = nucleotide.selected, sample.specific.features = sample.specific.features, continuous.sample.specific=continuous.sample.specific, min.count = min.count, genome.size = genome.size, cores = cores)

  mut.regions2 = mut.regions[names(mut.regions) %in% rownames(mut.rec)]
  mut.rec.hotspot = data.frame(chrom = as.character(GenomicRanges::seqnames(mut.regions2[rownames(mut.rec)])), start = GenomicRanges::start(mut.regions2[rownames(mut.rec)]), end = GenomicRanges::end(mut.regions2[rownames(mut.rec)]), mut.rec)
  
} else {
  
  # Redefine hotspots if not whole genome analysis
  regions = read.delim(region.of.interest, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  regions = with(regions, GenomicRanges::GRanges(V1, IRanges::IRanges(V2, V3)))
  regions = regions[as.character(GenomicRanges::seqnames(regions)) %in% as.character(GenomicRanges::seqnames(seqi))]
  names(regions) = paste("Region", c(1:length(regions)), sep = "")
  
  # Define masked region i.e. CDS, immunoglobulin loci and nonmappable
  mask.regions = readRDS(mask.regions.file)
  
  # Remove masked regions from region of interest
  maf.masked.regions <-regions[-S4Vectors::subjectHits(IRanges::findOverlaps(mask.regions, regions))]
  
  if (hotspots) {
    
  mut.rec <- mutPredict.snv.run.lr(roi = mut.regions, maf.snv = maf.snv, maf.snv2 = maf.snv2, model.snv = LRmodel.snv, continuous.features.snv = continuous.selected.features.snv, discrete.features.snv = discrete.selected.features.snv, motifs = sel.motif, nucleotide.selected = nucleotide.selected, sample.specific.features = sample.specific.features, continuous.sample.specific = continuous.sample.specific, min.count = min.count, genome.size = sum(GenomicRanges::width(maf.masked.regions)), cores = cores)
  
  mut.regions2 = mut.regions[names(mut.regions) %in% rownames(mut.rec)]
  mut.rec.hotspot = data.frame(chrom = as.character(GenomicRanges::seqnames(mut.regions2[rownames(mut.rec)])), start = GenomicRanges::start(mut.regions2[rownames(mut.rec)]), end = GenomicRanges::end(mut.regions2[rownames(mut.rec)]), mut.rec)
  
  } else {
    
    mut.rec <- mutPredict.snv.run.lr(roi = maf.masked.regions, maf.snv = maf.snv, maf.snv2 = maf.snv2, model.snv = LRmodel.snv, continuous.features.snv = continuous.selected.features.snv, discrete.features.snv = discrete.selected.features.snv, motifs = sel.motif, nucleotide.selected = nucleotide.selected, sample.specific.features = sample.specific.features, continuous.sample.specific=continuous.sample.specific, min.count = min.count, genome.size = length(maf.masked.regions), cores = cores)
    
    mut.regions2 = maf.masked.regions[names(maf.masked.regions) %in% rownames(mut.rec)]
    mut.rec.hotspot = data.frame(chrom = as.character(GenomicRanges::seqnames(mut.regions2[rownames(mut.rec)])), start = GenomicRanges::start(mut.regions2[rownames(mut.rec)]), end = GenomicRanges::end(mut.regions2[rownames(mut.rec)]), mut.rec)
    
  }
  
}

# merge overlapping hotspots in mut.rec.hotspot, assign smallest pvalue and recalculate k
if (merge.hotspots) {
  
  print("Merge overlapping hotspots in final results...")
  mut.rec.hotspot2=mut.rec.hotspot
  mut.rec.hotspot2$ID=as.character(rownames(mut.rec.hotspot2))
  mut.rec.hotspot2=with(mut.rec.hotspot2,GenomicRanges::GRanges(chrom,IRanges::IRanges(start,end),pval=pval, length=length,p.bg=p.bg,k=k,fdr=fdr,ID=ID))
  hotspots=IRanges::reduce(mut.rec.hotspot2)
  hotspots$hs=paste("hs",1:length(hotspots),sep="")
  ovl.mut=IRanges::findOverlaps(maf.snv,hotspots)
  hotspots2=hotspots[S4Vectors::subjectHits(ovl.mut)]
  hotspots2$sample=maf.snv[S4Vectors::queryHits(ovl.mut)]$sid
  hotspots2=GenomicRanges::as.data.frame(hotspots2)
  hotspots2=aggregate(sample~hs,hotspots2,FUN=function(k) length(unique(k)))
  colnames(hotspots2)[2]="k"
  rownames(hotspots2)=hotspots2$hs
  hotspots$k=0
  for (i in 1:length(hotspots)) {
    # print(i)
    hotspots$k[i]=hotspots2[which(hotspots2$hs==hotspots$hs[i]),"k"]
    
  }

  ovl=IRanges::findOverlaps(mut.rec.hotspot2,hotspots)
  mut.rec.hotspot2=mut.rec.hotspot2[S4Vectors::queryHits(ovl)]
  mut.rec.hotspot2$hs=hotspots[S4Vectors::subjectHits(ovl)]$hs
  mut.rec.hotspot2$region.start=IRanges::start(hotspots[S4Vectors::subjectHits(ovl)])
  mut.rec.hotspot2$region.end=IRanges::end(hotspots[S4Vectors::subjectHits(ovl)])
  mut.rec.hotspot2$new.k=hotspots[S4Vectors::subjectHits(ovl)]$k
  mut.rec.hotspot2=GenomicRanges::as.data.frame(mut.rec.hotspot2)
  mut.rec.hotspot2=mut.rec.hotspot2[order(mut.rec.hotspot2$pval,decreasing=FALSE),]
  mut.rec.hotspot2=mut.rec.hotspot2[!duplicated(mut.rec.hotspot2$hs),]
  mut.rec.hotspot2=mut.rec.hotspot2[,c("seqnames","region.start","region.end","pval","length","p.bg","new.k","fdr","ID")]
  mut.rec.hotspot2$length=mut.rec.hotspot2$region.end-mut.rec.hotspot2$region.start+1
  colnames(mut.rec.hotspot2)=c(colnames(mut.rec.hotspot),"ID")
  rownames(mut.rec.hotspot2)=mut.rec.hotspot2$ID
  mut.rec.hotspot2=mut.rec.hotspot2[,-ncol(mut.rec.hotspot2)]

  
  
  
  
} else {
  mut.rec.hotspot2=NULL
}


}

return(list(mut.rec.hotspot,mut.rec.hotspot2))

}
