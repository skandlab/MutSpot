library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(BSgenome)
library(stringi)
library(seqinr)
library(ggplot2)
library(stringr)
library(Biostrings)

###
# prep fasta files for dreme
# for 50bp
###
# get sequences from Hsapiens
gr=readRDS("enh_sites_RF_downsampled.rds") # for enhancers 27535
# gr=readRDS("UTR3_sites_RF_downsampled.rds") # for 3UTR 76369
# gr=readRDS("UTR5_sites_RF_downsampled.rds") # for 5UTR 21470
# gr=readRDS("prom_sites_RF_downsampled.rds") # for promoters 101505
# gr=readRDS("intron_sites_RF_downsampled.rds") # for introns 544323

table(gr$mut)

unique(width(gr))
range(start(gr)) 

# get 50 bp regions
gr2<-gr+25

# get seq for all sites in gr2
extrSeq=Views(Hsapiens,gr2)

dnaSeq=data.frame(seqnames=as.character(seqnames(gr2)),start=start(gr2),end=end(gr2),width=width(gr2),strand=as.character(strand(gr2)),mut=gr2$mut)
z=as(extrSeq,"DNAStringSet") 
dnaSeq=cbind(dnaSeq,as.character(z))
colnames(dnaSeq)[7]="dna"

# remove N from dnaString
N.string=which(stri_detect_fixed(dnaSeq$dna,"N")) 
table(dnaSeq[N.string,]$mut)
if(!is.null(N.string)){
  N.dnaSeq=dnaSeq[N.string,]
  dnaSeq=dnaSeq[-N.string,]
}

table(dnaSeq$mut)

mut=dnaSeq[which(dnaSeq$mut==1),] 
mut$dna=as.character(mut$dna)
non.mut=dnaSeq[which(dnaSeq$mut==0),]
non.mut$dna=as.character(non.mut$dna)

# write file in fasta format to be input files for dreme
name=c(1:nrow(mut))
name=paste("ID",name,sep="")
d=as.list(mut$dna)
names(d)=name
write.fasta(sequences=d,names=names(d),"mutated_final_enh.fa")
# write.fasta(sequences=d,names=names(d),"mutated_final_3utr.fa")
# write.faste(sequences=d,names=names(d),"mutated_final_5utr.fa")
# write.fasta(sequences=d,names=names(d),"mutated_final_prom.fa")
# write.fasta(sequences=d,names=names(d),"mutated_final_intron.fa")

name=c(1:nrow(non.mut))
name=paste("ID",name,sep="")
d=as.list(non.mut$dna)
names(d)=name
write.fasta(sequences=d,names=names(d),"non_mutated_final_enh.fa")
# write.fasta(sequences=d,names=names(d),"non_mutated_final_3utr.fa")
# write.fasta(sequences=d,names=names(d),"non_mutated_final_5utr.fa")
# write.fasta(sequences=d,names=names(d),"non_mutated_final_prom.fa")
# write.fasta(sequences=d,names=names(d),"non_mutated_final_intron.fa")

save.image(file="dreme_final_enh.RData")
save.image(file="dreme_final_3utr.RData")
save.image(file="dreme_final_5utr.RData")
save.image(file="dreme_final_prom.RData")
save.image(file="dreme_final_intron.RData")

###
# run DREME on 50bp region
###
# qsub -cwd -pe OpenMP 1 -l mem_free=4G,h_rt=300:00:00 /mnt/software/bin/dreme -oc /mnt/projects/changmm/data/dreme/results_enh -p /mnt/projects/changmm/data/dreme/mutated_final_enh.fa -n /mnt/projects/changmm/data/dreme/non_mutated_final_enh.fa -dna 
# qsub -cwd -pe OpenMP 1 -l mem_free=4G,h_rt=300:00:00 /mnt/software/bin/dreme -oc /mnt/projects/changmm/data/dreme/results_enhb -p /mnt/projects/changmm/data/dreme/non_mutated_final_enh.fa -n /mnt/projects/changmm/data/dreme/mutated_final_enh.fa -dna 

# qsub -cwd -pe OpenMP 1 -l mem_free=4G,h_rt=300:00:00 /mnt/software/bin/dreme -oc /mnt/projects/changmm/data/dreme/results_3utr -p /mnt/projects/changmm/data/dreme/mutated_final_3utr.fa -n /mnt/projects/changmm/data/dreme/non_mutated_final_3utr.fa -dna
# qsub -cwd -pe OpenMP 1 -l mem_free=4G,h_rt=300:00:00 /mnt/software/bin/dreme -oc /mnt/projects/changmm/data/dreme/results_3utrb -p /mnt/projects/changmm/data/dreme/non_mutated_final_3utr.fa -n /mnt/projects/changmm/data/dreme/mutated_final_3utr.fa -dna 

# qsub -cwd -pe OpenMP 1 -l mem_free=4G,h_rt=300:00:00 /mnt/software/bin/dreme -oc /mnt/projects/changmm/data/dreme/resuts_5utr -p /mnt/projects/changmm/data/dreme/mutated_final_5utr.fa -n /mnt/projects/changmm/data/dreme/non_mutated_final_5utr.fa -dna
# qsub -cwd -pe OpenMP 1 -l mem_free=4G,h_rt=300:00:00 /mnt/software/bin/dreme -oc /mnt/projects/changmm/data/dreme/results_5utrb -p /mnt/projects/changmm/data/dreme/non_mutated_final_5utr.fa -n /mnt/projects/changmm/data/dreme/mutated_final_5utr.fa -dna

# qsub -cwd -pe OpenMP 1 -l mem_free=4G,h_rt=300:00:00 /mnt/software/bin/dreme -oc /mnt/projects/changmm/data/dreme/results_prom -p /mnt/projects/changmm/data/dreme/mutated_final_prom.fa -n /mnt/projects/changmm/data/dreme/non_mutated_final_prom.fa -dna
# qsub -cwd -pe OpenMP 1 -l mem_free=4G,h_rt=300:00:00 /mnt/software/bin/dreme -oc /mnt/projects/changmm/data/dreme/results_promb -p /mnt/projects/changmm/data/dreme/non_mutated_final_prom.fa -n /mnt/projects/changmm/data/dreme/mutated_final_prom.fa -dna

# qsub -cwd -pe OpenMP 1 -l mem_free=4G,h_rt=300:00:00 /mnt/software/bin/dreme -oc /mnt/projects/changmm/data/dreme/results_intron -p /mnt/projects/changmm/data/dreme/mutated_final_intron.fa -n /mnt/projects/changmm/data/dreme/non_mutated_final_intron.fa -dna
# qsub -cwd -pe OpenMP 1 -l mem_free=4G,h_rt=300:00:00 /mnt/software/bin/dreme -oc /mnt/projects/changmm/data/dreme/results_intronb -p /mnt/projects/changmm/data/dreme/non_mutated_final_intron.fa -n /mnt/projects/changmm/data/dreme/mutated_final_intron.fa -dna

########################################################################################################################

###
# plot freq figures for motifs found in dreme
# get positions of mutated motifs in mutated
###
# enhancers
B=c("C","G","T")
R=c("A","G")
# motif 1: BRCG - B:C,G,T; R:A,G
brcg.motif=NULL
for (i in B){
  for (j in R){
    brcg.motif=c(brcg.motif,paste(i,j,"CG",sep=""))
  }
}
brcg=apply(mut,1,function(x) {
  stri_locate_all_fixed(x[7],brcg.motif)
})
brcg=lapply(brcg,FUN=function(x) do.call(rbind,x))
brcg=do.call(rbind,brcg)
brcg=brcg[!is.na(brcg[,1]),]
brcg=as.data.frame(brcg)
# 4527 2
mean(brcg$start) # 24.33223
mean(brcg$end) # 27.33223

# motif 1.rc: CGR'B' - R':T,C; B':G,C,A
brcg.motif.rc=lapply(brcg.motif,FUN=function(x) {
  t=chartr("ATCG","TAGC",x)
  t=reverse(t)
  t
})
brcg.motif.rc=unlist(brcg.motif.rc)
brcg.rc=apply(mut,1,function(x) {
  stri_locate_all_fixed(x[7],brcg.motif.rc)
})
brcg.rc=lapply(brcg.rc,FUN=function(x) do.call(rbind,x))
brcg.rc=do.call(rbind,brcg.rc)
brcg.rc=brcg.rc[!is.na(brcg.rc[,1]),]
brcg.rc=as.data.frame(brcg.rc)
# 4447 2
mean(brcg.rc$start) # 24.85159
mean(brcg.rc$end) # 27.85159

# motif 2: BCG - B:C,G,T
bcg.motif=NULL
for (i in B){
  bcg.motif=c(bcg.motif,paste(i,"CG",sep=""))
}
bcg=apply(mut,1,function(x) {
  stri_locate_all_fixed(x[7],bcg.motif)
})
bcg=lapply(bcg,FUN=function(x) do.call(rbind,x))
bcg=do.call(rbind,bcg)
bcg=bcg[!is.na(bcg[,1]),]
bcg=as.data.frame(bcg)
# 8389 2
mean(bcg$start) # 24.89343
mean(bcg$end) # 26.89343

# motif 2.rc: CGB' - B':G,C,A
bcg.motif.rc=lapply(bcg.motif,FUN=function(x) {
  t=chartr("ATCG","TAGC",x)
  t=reverse(t)
  t
})
bcg.motif.rc=unlist(bcg.motif.rc)
bcg.rc=apply(mut,1,function(x) {
  stri_locate_all_fixed(x[7],bcg.motif.rc)
})
bcg.rc=lapply(bcg.rc,FUN=function(x) do.call(rbind,x))
bcg.rc=do.call(rbind,bcg.rc)
bcg.rc=bcg.rc[!is.na(bcg.rc[,1]),]
bcg.rc=as.data.frame(bcg.rc)
# 8480 2
mean(bcg.rc$start) # 25.04564
mean(bcg.rc$end) # 27.04564

Freq.a=table(brcg$start)
Freq.a=data.frame(Freq.a,motif="BRCG")
colnames(Freq.a)=c("position","freq","motif")

tab=table(brcg.rc$start)
tab=data.frame(tab,motif="BRCG(RC)")
colnames(tab)=c("position","freq","motif")
Freq.a=rbind(Freq.a,tab)
# 96 3

Freq.b=table(bcg$start)
Freq.b=data.frame(Freq.b,motif="BCG")
colnames(Freq.b)=c("position","freq","motif")

tab=table(bcg.rc$start)
tab=data.frame(tab,motif="BCG(RC)")
colnames(tab)=c("position","freq","motif")
Freq.b=rbind(Freq.b,tab)
# 98 3

# 3UTR
# B=c("C","G","T")
# H=c("A","C","T")
# R=c("A","G")
# K=c("G","T")
# # motif 1: CGBH - B:C,G,T; H:A,C,T
# cgbh.motif=NULL
# for (i in B){
#   for (j in H){
#     cgbh.motif=c(cgbh.motif,paste("CG",i,j,sep=""))
#   }
# }
# cgbh=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],cgbh.motif)
# })
# cgbh=lapply(cgbh,FUN=function(x) do.call(rbind,x))
# cgbh=do.call(rbind,cgbh)
# cgbh=cgbh[!is.na(cgbh[,1]),]
# cgbh=as.data.frame(cgbh)
# # 13795 2
# mean(cgbh$start) # 24.85335
# mean(cgbh$end) # 27.85335
# 
# # motif 1.rc: H'B'CG - H':T,G,A; B':G,C,A
# cgbh.motif.rc=lapply(cgbh.motif,FUN=function(x) {
#   t=chartr("ATCG","TAGC",x)
#   t=reverse(t)
#   t
# })
# cgbh.motif.rc=unlist(cgbh.motif.rc)
# cgbh.rc=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],cgbh.motif.rc)
# })
# cgbh.rc=lapply(cgbh.rc,FUN=function(x) do.call(rbind,x))
# cgbh.rc=do.call(rbind,cgbh.rc)
# cgbh.rc=cgbh.rc[!is.na(cgbh.rc[,1]),]
# cgbh.rc=as.data.frame(cgbh.rc)
# # 13876 2
# mean(cgbh.rc$start) # 24.34289
# mean(cgbh.rc$end) # 27.34289
# 
# # motif 2: ACACACRC - R:A,G
# acacacrc.motif=NULL
# for (i in R){
#   acacacrc.motif=c(acacacrc.motif,paste("ACACAC",i,"C",sep=""))
# }
# acacacrc=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],acacacrc.motif)
# })
# acacacrc=lapply(acacacrc,FUN=function(x) do.call(rbind,x))
# acacacrc=do.call(rbind,acacacrc)
# acacacrc=acacacrc[!is.na(acacacrc[,1]),]
# acacacrc=as.data.frame(acacacrc)
# # 368 2
# mean(acacacrc$start) # 21.11141
# mean(acacacrc$end) # 28.11141
# 
# # motif 2.rc: GR'GTGTGT - R':T,C
# acacacrc.motif.rc=lapply(acacacrc.motif,FUN=function(x) {
#   t=chartr("ATCG","TAGC",x)
#   t=reverse(t)
#   t
# })
# acacacrc.motif.rc=unlist(acacacrc.motif.rc)
# acacacrc.rc=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],acacacrc.motif.rc)
# })
# acacacrc.rc=lapply(acacacrc.rc,FUN=function(x) do.call(rbind,x))
# acacacrc.rc=do.call(rbind,acacacrc.rc)
# acacacrc.rc=acacacrc.rc[!is.na(acacacrc.rc[,1]),]
# acacacrc.rc=as.data.frame(acacacrc.rc)
# # 372 2
# mean(acacacrc.rc$start) # 20.57796
# mean(acacacrc.rc$end) # 27.57796
# 
# # motif 3: CGH - H:A,C,T
# cgh.motif=NULL
# for (i in H){
#   cgh.motif=c(cgh.motif,paste("CG",i,sep=""))
# }
# cgh=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],cgh.motif)
# })
# cgh=lapply(cgh,FUN=function(x) do.call(rbind,x))
# cgh=do.call(rbind,cgh)
# cgh=cgh[!is.na(cgh[,1]),]
# cgh=as.data.frame(cgh)
# # 17495 2
# mean(cgh$start) # 25.04819
# mean(cgh$end) # 27.04819
# 
# # motif 3.rc: H'CG - H':T,G,A
# cgh.motif.rc=lapply(cgh.motif,FUN=function(x) {
#   t=chartr("ATCG","TAGC",x)
#   t=reverse(t)
#   t
# })
# cgh.motif.rc=unlist(cgh.motif.rc)
# cgh.rc=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],cgh.motif.rc)
# })
# cgh.rc=lapply(cgh.rc,FUN=function(x) do.call(rbind,x))
# cgh.rc=do.call(rbind,cgh.rc)
# cgh.rc=cgh.rc[!is.na(cgh.rc[,1]),]
# cgh.rc=as.data.frame(cgh.rc)
# # 17321 2
# mean(cgh.rc$start) # 25.01039
# mean(cgh.rc$end) # 27.01039
# 
# # motif 4: KAAAAAAA - K:G,T
# kaaaaaaa.motif=NULL
# for (i in K){
#   kaaaaaaa.motif=c(kaaaaaaa.motif,paste(i,"AAAAAAA",sep=""))
# }
# kaaaaaaa=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],kaaaaaaa.motif)
# })
# kaaaaaaa=lapply(kaaaaaaa,FUN=function(x) do.call(rbind,x))
# kaaaaaaa=do.call(rbind,kaaaaaaa)
# kaaaaaaa=kaaaaaaa[!is.na(kaaaaaaa[,1]),]
# kaaaaaaa=as.data.frame(kaaaaaaa)
# # 876 2
# mean(kaaaaaaa$start) # 22.62329
# mean(kaaaaaaa$end) # 29.62329
# 
# # motif 4.rc: TTTTTTTK' - K':C,A
# kaaaaaaa.motif.rc=lapply(kaaaaaaa.motif,FUN=function(x) {
#   t=chartr("ATCG","TAGC",x)
#   t=reverse(t)
#   t
# })
# kaaaaaaa.motif.rc=unlist(kaaaaaaa.motif.rc)
# kaaaaaaa.rc=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],kaaaaaaa.motif.rc)
# })
# kaaaaaaa.rc=lapply(kaaaaaaa.rc,FUN=function(x) do.call(rbind,x))
# kaaaaaaa.rc=do.call(rbind,kaaaaaaa.rc)
# kaaaaaaa.rc=kaaaaaaa.rc[!is.na(kaaaaaaa.rc[,1]),]
# kaaaaaaa.rc=as.data.frame(kaaaaaaa.rc)
# # 843 2
# mean(kaaaaaaa.rc$start) # 22.93832
# mean(kaaaaaaa.rc$end) # 29.93832
# 
# # motif 5: ATATATAT
# atatatat.motif=c("ATATATAT")
# atatatat=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],atatatat.motif)
# })
# atatatat=lapply(atatatat,FUN=function(x) do.call(rbind,x))
# atatatat=do.call(rbind,atatatat)
# atatatat=atatatat[!is.na(atatatat[,1]),]
# atatatat=as.data.frame(atatatat)
# # 387 2
# mean(atatatat$start) # 22.00517
# mean(atatatat$end) # 29.00517

# Freq.a=table(cgbh$start)
# Freq.a=data.frame(Freq.a,motif="CGBH")
# colnames(Freq.a)=c("position","freq","motif")
# 
# tab=table(cgbh.rc$start)
# tab=data.frame(tab,motif="CGBH(RC)")
# colnames(tab)=c("position","freq","motif")
# Freq.a=rbind(Freq.a,tab)
# # 96 3
# 
# Freq.b=table(acacacrc$start)
# Freq.b=data.frame(Freq.b,motif="ACACACRC")
# colnames(Freq.b)=c("position","freq","motif")
# 
# tab=table(acacacrc.rc$start)
# tab=data.frame(tab,motif="ACACACRC(RC)")
# colnames(tab)=c("position","freq","motif")
# Freq.b=rbind(Freq.b,tab)
# # 88 3
# 
# Freq.c=table(cgh$start)
# Freq.c=data.frame(Freq.c,motif="CGH")
# colnames(Freq.c)=c("position","freq","motif")
# 
# tab=table(cgh.rc$start)
# tab=data.frame(tab,motif="CGH(RC)")
# colnames(tab)=c("position","freq","motif")
# Freq.c=rbind(Freq.c,tab)
# # 98 3
# 
# Freq.d=table(kaaaaaaa$start)
# Freq.d=data.frame(Freq.d,motif="KAAAAAAA")
# colnames(Freq.d)=c("position","freq","motif")
# 
# tab=table(kaaaaaaa.rc$start)
# tab=data.frame(tab,motif="KAAAAAAA(RC)")
# colnames(tab)=c("position","freq","motif")
# Freq.d=rbind(Freq.d,tab)
# # 88 3
# 
# Freq.e=table(atatatat$start)
# Freq.e=data.frame(Freq.e,motif="ATATATAT")
# colnames(Freq.e)=c("position","freq","motif")
# 44 3

# 5UTR
# H=c("A","C","T")
# # motif 1: CGH - H:A,C,T
# cgh.motif=NULL
# for (i in H){
#   cgh.motif=c(cgh.motif,paste("CG",i,sep=""))
# }
# cgh=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],cgh.motif)
# })
# cgh=lapply(cgh,FUN=function(x) do.call(rbind,x))
# cgh=do.call(rbind,cgh)
# cgh=cgh[!is.na(cgh[,1]),]
# cgh=as.data.frame(cgh)
# # 15531 2
# mean(cgh$start) # 25.03277
# mean(cgh$end) # 27.03277
# 
# # motif 1.rc: H'CG - H':T,G,A
# cgh.motif.rc=lapply(cgh.motif,FUN=function(x) {
#   t=chartr("ATCG","TAGC",x)
#   t=reverse(t)
#   t
# })
# cgh.motif.rc=unlist(cgh.motif.rc)
# cgh.rc=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],cgh.motif.rc)
# })
# cgh.rc=lapply(cgh.rc,FUN=function(x) do.call(rbind,x))
# cgh.rc=do.call(rbind,cgh.rc)
# cgh.rc=cgh.rc[!is.na(cgh.rc[,1]),]
# cgh.rc=as.data.frame(cgh.rc)
# # 15747 2
# mean(cgh.rc$start) # 24.99276
# mean(cgh.rc$end) # 26.99276
# 
# Freq=table(cgh$start)
# Freq=data.frame(Freq,motif="CGH")
# colnames(Freq)=c("position","freq","motif")
# 
# tab=table(cgh.rc$start)
# tab=data.frame(tab,motif="CGH(RC)")
# colnames(tab)=c("position","freq","motif")
# Freq=rbind(Freq,tab)
# 98 3

# promoter
# D=c("A","G","T")
# Y=c("C","T")
# W=c("A","T")
# #motif 1: DCGY - D:A,G,T; Y:C,T
# dcgy.motif=NULL
# for (i in D){
#   for (j in Y){
#     dcgy.motif=c(dcgy.motif,paste(i,"CG",j,sep=""))
#   }
# }
# dcgy=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],dcgy.motif)
# })
# dcgy=lapply(dcgy,FUN=function(x) do.call(rbind,x))
# dcgy=do.call(rbind,dcgy)
# dcgy=dcgy[!is.na(dcgy[,1]),]
# dcgy=as.data.frame(dcgy)
# # 25263 2
# mean(dcgy$start) # 24.58722
# mean(dcgy$end) # 27.58722
# 
# # motif 1.rc: Y'CGD' - Y':G,A; D':T,C,A
# dcgy.motif.rc=lapply(dcgy.motif,FUN=function(x) {
#   t=chartr("ATCG","TAGC",x)
#   t=reverse(t)
#   t
# })
# dcgy.motif.rc=unlist(dcgy.motif.rc)
# dcgy.rc=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],dcgy.motif.rc)
# })
# dcgy.rc=lapply(dcgy.rc,FUN=function(x) do.call(rbind,x))
# dcgy.rc=do.call(rbind,dcgy.rc)
# dcgy.rc=dcgy.rc[!is.na(dcgy.rc[,1]),]
# dcgy.rc=as.data.frame(dcgy.rc)
# # 25184 2
# mean(dcgy.rc$start) # 24.45021
# mean(dcgy.rc$end) # 27.45021
# 
# # motif 2: CGD - D:A,G,T
# cgd.motif=NULL
# for (i in D){
#   cgd.motif=c(cgd.motif,paste("CG",i,sep=""))
# }
# cgd=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],cgd.motif)
# })
# cgd=lapply(cgd,FUN=function(x) do.call(rbind,x))
# cgd=do.call(rbind,cgd)
# cgd=cgd[!is.na(cgd[,1]),]
# cgd=as.data.frame(cgd)
# # 49542 2
# mean(cgd$start) # 25.05074
# mean(cgd$end) # 27.05074
# 
# # motif 2.rc: D'CG - D':T,C,A
# cgd.motif.rc=lapply(cgd.motif,FUN=function(x) {
#   t=chartr("ATCG","TAGC",x)
#   t=reverse(t)
#   t
# })
# cgd.motif.rc=unlist(cgd.motif.rc)
# cgd.rc=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],cgd.motif.rc)
# })
# cgd.rc=lapply(cgd.rc,FUN=function(x) do.call(rbind,x))
# cgd.rc=do.call(rbind,cgd.rc)
# cgd.rc=cgd.rc[!is.na(cgd.rc[,1]),]
# cgd.rc=as.data.frame(cgd.rc)
# # 49223 2
# mean(cgd.rc$start) # 24.94206
# mean(cgd.rc$end) # 26.94206
# 
# # motif 3: AAAAAAWT - W:A,T
# aaaaaawt.motif=NULL
# for (i in W){
#   aaaaaawt.motif=c(aaaaaawt.motif,paste("AAAAAA",i,"T",sep=""))
# }
# aaaaaawt=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],aaaaaawt.motif)
# })
# aaaaaawt=lapply(aaaaaawt,FUN=function(x) do.call(rbind,x))
# aaaaaawt=do.call(rbind,aaaaaawt)
# aaaaaawt=aaaaaawt[!is.na(aaaaaawt[,1]),]
# aaaaaawt=as.data.frame(aaaaaawt)
# # 830 2
# mean(aaaaaawt$start) # 21.36386
# mean(aaaaaawt$end) # 28.36386
# 
# # motif 3.rc: AW'TTTTTT - W':T,A
# aaaaaawt.motif.rc=lapply(aaaaaawt.motif,FUN=function(x) {
#   t=chartr("ATCG","TAGC",x)
#   t=reverse(t)
#   t
# })
# aaaaaawt.motif.rc=unlist(aaaaaawt.motif.rc)
# aaaaaawt.rc=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],aaaaaawt.motif.rc)
# })
# aaaaaawt.rc=lapply(aaaaaawt.rc,FUN=function(x) do.call(rbind,x))
# aaaaaawt.rc=do.call(rbind,aaaaaawt.rc)
# aaaaaawt.rc=aaaaaawt.rc[!is.na(aaaaaawt.rc[,1]),]
# aaaaaawt.rc=as.data.frame(aaaaaawt.rc)
# # 828 2
# mean(aaaaaawt.rc$start) # 22.32367
# mean(aaaaaawt.rc$end) # 29.32367
# 
# # motif 4: AYACACAC - Y:C,T
# ayacacac.motif=NULL
# for (i in Y){
#   ayacacac.motif=c(ayacacac.motif,paste("A",i,"ACACAC",sep=""))
# }
# ayacacac=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],ayacacac.motif)
# })
# ayacacac=lapply(ayacacac,FUN=function(x) do.call(rbind,x))
# ayacacac=do.call(rbind,ayacacac)
# ayacacac=ayacacac[!is.na(ayacacac[,1]),]
# ayacacac=as.data.frame(ayacacac)
# # 500 2
# mean(ayacacac$start) # 21.358
# mean(ayacacac$end) # 28.358
# 
# # motif 4.rc: GTGTGTY'T - Y':G,A
# ayacacac.motif.rc=lapply(ayacacac.motif,FUN=function(x) {
#   t=chartr("ATCG","TAGC",x)
#   t=reverse(t)
#   t
# })
# ayacacac.motif.rc=unlist(ayacacac.motif.rc)
# ayacacac.rc=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],ayacacac.motif.rc)
# })
# ayacacac.rc=lapply(ayacacac.rc,FUN=function(x) do.call(rbind,x))
# ayacacac.rc=do.call(rbind,ayacacac.rc)
# ayacacac.rc=ayacacac.rc[!is.na(ayacacac.rc[,1]),]
# ayacacac.rc=as.data.frame(ayacacac.rc)
# # 531 2
# mean(ayacacac.rc$start) # 21.61394
# mean(ayacacac.rc$end) # 28.61394
# 
# # motif 5: GGGGAAAA
# ggggaaaa.motif=c("GGGGAAAA")
# ggggaaaa=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],ggggaaaa.motif)
# })
# ggggaaaa=lapply(ggggaaaa,FUN=function(x) do.call(rbind,x))
# ggggaaaa=do.call(rbind,ggggaaaa)
# ggggaaaa=ggggaaaa[!is.na(ggggaaaa[,1]),]
# ggggaaaa=as.data.frame(ggggaaaa)
# # 133 2
# mean(ggggaaaa$start) # 20.7218
# mean(ggggaaaa$end) # 27.7218
# 
# # motif 5.rc: TTTTCCCC
# ggggaaaa.motif.rc=lapply(ggggaaaa.motif,FUN=function(x) {
#   t=chartr("ATCG","TAGC",x)
#   t=reverse(t)
#   t
# })
# ggggaaaa.motif.rc=unlist(ggggaaaa.motif.rc)
# ggggaaaa.rc=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],ggggaaaa.motif.rc)
# })
# ggggaaaa.rc=lapply(ggggaaaa.rc,FUN=function(x) do.call(rbind,x))
# ggggaaaa.rc=do.call(rbind,ggggaaaa.rc)
# ggggaaaa.rc=ggggaaaa.rc[!is.na(ggggaaaa.rc[,1]),]
# ggggaaaa.rc=as.data.frame(ggggaaaa.rc)
# # 147 2
# mean(ggggaaaa.rc$start) # 22.63265
# mean(ggggaaaa.rc$end) # 29.63265
# 
# # motif 6: ATATATAT
# atatatat.motif=c("ATATATAT")
# atatatat=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],atatatat.motif)
# })
# atatatat=lapply(atatatat,FUN=function(x) do.call(rbind,x))
# atatatat=do.call(rbind,atatatat)
# atatatat=atatatat[!is.na(atatatat[,1]),]
# atatatat=as.data.frame(atatatat)
# # # 282 2
# # mean(atatatat$start) # 21.95035
# # mean(atatatat$end) # 28.95035
# 
# Freq.a=table(dcgy$start)
# Freq.a=data.frame(Freq.a,motif="DCGY")
# colnames(Freq.a)=c("position","freq","motif")
# 
# tab=table(dcgy.rc$start)
# tab=data.frame(tab,motif="DCGY(RC)")
# colnames(tab)=c("position","freq","motif")
# Freq.a=rbind(Freq.a,tab)
# # 96 3
# 
# Freq.b=table(cgd$start)
# Freq.b=data.frame(Freq.b,motif="CGD")
# colnames(Freq.b)=c("position","freq","motif")
# 
# tab=table(cgd.rc$start)
# tab=data.frame(tab,motif="CGD(RC)")
# colnames(tab)=c("position","freq","motif")
# Freq.b=rbind(Freq.b,tab)
# # 98 3
# 
# Freq.c=table(aaaaaawt$start)
# Freq.c=data.frame(Freq.c,motif="AAAAAAWT")
# colnames(Freq.c)=c("position","freq","motif")
# 
# tab=table(aaaaaawt.rc$start)
# tab=data.frame(tab,motif="AAAAAAWT(RC)")
# colnames(tab)=c("position","freq","motif")
# Freq.c=rbind(Freq.c,tab)
# # 88 3
# 
# Freq.d=table(ayacacac$start)
# Freq.d=data.frame(Freq.d,motif="AYACACAC")
# colnames(Freq.d)=c("position","freq","motif")
# 
# tab=table(ayacacac.rc$start)
# tab=data.frame(tab,motif="AYACACAC(RC)")
# colnames(tab)=c("position","freq","motif")
# Freq.d=rbind(Freq.d,tab)
# # 88 3
# 
# Freq.e=table(ggggaaaa$start)
# Freq.e=data.frame(Freq.e,motif="GGGGAAAA")
# colnames(Freq.e)=c("position","freq","motif")
# 
# tab=table(ggggaaaa.rc$start)
# tab=data.frame(tab,motif="GGGGAAAA(RC)")
# colnames(tab)=c("position","freq","motif")
# Freq.e=rbind(Freq.e,tab)
# # 77 3
# 
# Freq.f=table(atatatat$start)
# Freq.f=data.frame(Freq.f,motif="ATATATAT")
# colnames(Freq.f)=c("position","freq","motif")
# 44 3

# intron
# D=c("A","G","T")
# N=c("A","C","G","T")
# Y=c("C","T")
# V=c("A","C","G")
# R=c("A","G")
# H=c("A","C","T")
# K=c("G","T")
# B=c("C","G","T")
# W=c("A","T")
# # motif 1: DNCGY - D:A,G,T; N:A,C,G,T; Y:C,T
# dncgy.motif=NULL
# for (i in D){
#   for (j in N){
#     for (s in Y){
#       dncgy.motif=c(dncgy.motif,paste(i,j,"CG",s,sep=""))
#     }
#   }
# }
# dncgy=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],dncgy.motif)
# })
# dncgy=lapply(dncgy,FUN=function(x) do.call(rbind,x))
# dncgy=do.call(rbind,dncgy)
# dncgy=dncgy[!is.na(dncgy[,1]),]
# dncgy=as.data.frame(dncgy)
# # 60740 2
# mean(dncgy$start) # 23.72242
# mean(dncgy$end) # 27.72242
# 
# # motif 1.rc: Y'CGN'D' - Y':G,A; N':T,G,C,A; D':T,C,A
# dncgy.motif.rc=lapply(dncgy.motif,FUN=function(x) {
#   t=chartr("ATCG","TAGC",x)
#   t=reverse(t)
#   t
# })
# dncgy.motif.rc=unlist(dncgy.motif.rc)
# dncgy.rc=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],dncgy.motif.rc)
# })
# dncgy.rc=lapply(dncgy.rc,FUN=function(x) do.call(rbind,x))
# dncgy.rc=do.call(rbind,dncgy.rc)
# dncgy.rc=dncgy.rc[!is.na(dncgy.rc[,1]),]
# dncgy.rc=as.data.frame(dncgy.rc)
# # 60217 2
# mean(dncgy.rc$start) # 24.2389
# mean(dncgy.rc$end) # 28.2389
# 
# # motif 2: CGVR - V:A,C,G; R:A,G
# cgvr.motif=NULL
# for (i in V){
#   for (j in R){
#     cgvr.motif=c(cgvr.motif,paste("CG",i,j,sep=""))
#   }
# }
# cgvr=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],cgvr.motif)
# })
# cgvr=lapply(cgvr,FUN=function(x) do.call(rbind,x))
# cgvr=do.call(rbind,cgvr)
# cgvr=cgvr[!is.na(cgvr[,1]),]
# cgvr=as.data.frame(cgvr)
# # 63343 2
# mean(cgvr$start) # 24.58215
# mean(cgvr$end) # 27.58215
# 
# # motif 2.rc: R'V'CG - R':T,C; V':T,G,C
# cgvr.motif.rc=lapply(cgvr.motif,FUN=function(x) {
#   t=chartr("ATCG","TAGC",x)
#   t=reverse(t)
#   t
# })
# cgvr.motif.rc=unlist(cgvr.motif.rc)
# cgvr.rc=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],cgvr.motif.rc)
# })
# cgvr.rc=lapply(cgvr.rc,FUN=function(x) do.call(rbind,x))
# cgvr.rc=do.call(rbind,cgvr.rc)
# cgvr.rc=cgvr.rc[!is.na(cgvr.rc[,1]),]
# cgvr.rc=as.data.frame(cgvr.rc)
# # 64246 2
# mean(cgvr.rc$start) # 24.3304
# mean(cgvr.rc$end) # 27.3304
# 
# # motif 3: AHATATAY - H:A,C,T; Y:C,T
# ahatatay.motif=NULL
# for (i in H){
#   for (j in Y){
#     ahatatay.motif=c(ahatatay.motif,paste("A",i,"ATATA",j,sep=""))
#   }
# }
# ahatatay=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],ahatatay.motif)
# })
# ahatatay=lapply(ahatatay,FUN=function(x) do.call(rbind,x))
# ahatatay=do.call(rbind,ahatatay)
# ahatatay=ahatatay[!is.na(ahatatay[,1]),]
# ahatatay=as.data.frame(ahatatay)
# # 7050
# mean(ahatatay$start) # 22.07546
# mean(ahatatay$end) # 29.07546
# 
# # motif 3.rc: Y'TATATH'T - Y':G,A; H':T,G,A
# ahatatay.motif.rc=lapply(ahatatay.motif,FUN=function(x) {
#   t=chartr("ATCG","TAGC",x)
#   t=reverse(t)
#   t
# })
# ahatatay.motif.rc=unlist(ahatatay.motif.rc)
# ahatatay.rc=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],ahatatay.motif)
# })
# ahatatay.rc=lapply(ahatatay.rc,FUN=function(x) do.call(rbind,x))
# ahatatay.rc=do.call(rbind,ahatatay.rc)
# ahatatay.rc=ahatatay.rc[!is.na(ahatatay.rc[,1]),]
# ahatatay.rc=as.data.frame(ahatatay.rc)
# # 7054 2
# mean(ahatatay.rc$start) # 22.24029
# mean(ahatatay.rc$end) # 29.24029
# 
# # motif 4: KKAAAAAA - K:G,T
# kkaaaaaa.motif=NULL
# for (i in K){
#   for (j in K){
#     kkaaaaaa.motif=c(kkaaaaaa.motif,paste(i,j,"AAAAAA",sep=""))
#   }
# }
# kkaaaaaa=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],kkaaaaaa.motif)
# })
# kkaaaaaa=lapply(kkaaaaaa,FUN=function(x) do.call(rbind,x))
# kkaaaaaa=do.call(rbind,kkaaaaaa)
# kkaaaaaa=kkaaaaaa[!is.na(kkaaaaaa[,1]),]
# kkaaaaaa=as.data.frame(kkaaaaaa)
# # 6381 2
# mean(kkaaaaaa$start) # 22.05626
# mean(kkaaaaaa$end) # 29.05626
# 
# # motif 4.rc: TTTTTTK'K' - K':C,A
# kkaaaaaa.motif.rc=lapply(kkaaaaaa.motif,FUN=function(x) {
#   t=chartr("ATCG","TAGC",x)
#   t=reverse(t)
#   t
# })
# kkaaaaaa.motif.rc=unlist(kkaaaaaa.motif.rc)
# kkaaaaaa.rc=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],kkaaaaaa.motif.rc)
# })
# kkaaaaaa.rc=lapply(kkaaaaaa.rc,FUN=function(x) do.call(rbind,x))
# kkaaaaaa.rc=do.call(rbind,kkaaaaaa.rc)
# kkaaaaaa.rc=kkaaaaaa.rc[!is.na(kkaaaaaa.rc[,1]),]
# kkaaaaaa.rc=as.data.frame(kkaaaaaa.rc)
# # 6300 2
# mean(kkaaaaaa.rc$start) # 22.9127
# mean(kkaaaaaa.rc$end) # 29.9127
# 
# # motif 5: ACACACRC - R:A,G
# acacacrc.motif=NULL
# for (i in R){
#   acacacrc.motif=c(acacacrc.motif,paste("ACACAC",i,"C",sep=""))
# }
# acacacrc=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],acacacrc.motif)
# })
# acacacrc=lapply(acacacrc,FUN=function(x) do.call(rbind,x))
# acacacrc=do.call(rbind,acacacrc)
# acacacrc=acacacrc[!is.na(acacacrc[,1]),]
# acacacrc=as.data.frame(acacacrc)
# # 2434 2
# mean(acacacrc$start) # 21.9129
# mean(acacacrc$end) # 28.9129
# 
# # motif 5.rc: GR'GTGTGT - R':T,C
# acacacrc.motif.rc=lapply(acacacrc.motif,FUN=function(x) {
#   t=chartr("ATCG","TAGC",x)
#   t=reverse(t)
#   t
# })
# acacacrc.motif.rc=unlist(acacacrc.motif.rc)
# acacacrc.rc=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],acacacrc.motif.rc)
# })
# acacacrc.rc=lapply(acacacrc.rc,FUN=function(x) do.call(rbind,x))
# acacacrc.rc=do.call(rbind,acacacrc.rc)
# acacacrc.rc=acacacrc.rc[!is.na(acacacrc.rc[,1]),]
# acacacrc.rc=as.data.frame(acacacrc.rc)
# # 2424 2
# mean(acacacrc.rc$start) # 20.82632
# mean(acacacrc.rc$end) # 27.82632
# 
# # motif 6: AADTTTTK - D:A,G,T; K:G,T
# aadttttk.motif=NULL
# for (i in D){
#   for (j in K){
#     aadttttk.motif=c(aadttttk.motif,paste("AA",i,"TTTT",j,sep=""))
#   }
# }
# aadttttk=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],aadttttk.motif)
# })
# aadttttk=lapply(aadttttk,FUN=function(x) do.call(rbind,x))
# aadttttk=do.call(rbind,aadttttk)
# aadttttk=aadttttk[!is.na(aadttttk[,1]),]
# aadttttk=as.data.frame(aadttttk)
# # 8228 2
# mean(aadttttk$start) # 23.28439
# mean(aadttttk$end) # 30.28439
# 
# # motif 6.rc: K'AAAAD'TT - K':C,A; D':T,C,A
# aadttttk.motif.rc=lapply(aadttttk.motif,FUN=function(x) {
#   t=chartr("ATCG","TAGC",x)
#   t=reverse(t)
#   t
# })
# aadttttk.motif.rc=unlist(aadttttk.motif.rc)
# aadttttk.rc=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],aadttttk.motif.rc)
# })
# aadttttk.rc=lapply(aadttttk.rc,FUN=function(x) do.call(rbind,x))
# aadttttk.rc=do.call(rbind,aadttttk.rc)
# aadttttk.rc=aadttttk.rc[!is.na(aadttttk.rc[,1]),]
# aadttttk.rc=as.data.frame(aadttttk.rc)
# # 8418
# mean(aadttttk.rc$start) # 21.55809
# mean(aadttttk.rc$end) # 28.55809
# 
# # motif 7: CGB - B:C,G,T
# cgb.motif=NULL
# for (i in B){
#   cgb.motif=c(cgb.motif,paste("CG",i,sep=""))
# }
# cgb=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],cgb.motif)
# })
# cgb=lapply(cgb,FUN=function(x) do.call(rbind,x))
# cgb=do.call(rbind,cgb)
# cgb=cgb[!is.na(cgb[,1]),]
# cgb=as.data.frame(cgb)
# # 135928 2
# mean(cgb$start) # 25.08666
# mean(cgb$end) # 27.08666
# 
# # motif 7.rc: B'CG - B':G,C,A
# cgb.motif.rc=lapply(cgb.motif,FUN=function(x) {
#   t=chartr("ATCG","TAGC",x)
#   t=reverse(t)
#   t
# })
# cgb.motif.rc=unlist(cgb.motif.rc)
# cgb.rc=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],cgb.motif.rc)
# })
# cgb.rc=lapply(cgb.rc,FUN=function(x) do.call(rbind,x))
# cgb.rc=do.call(rbind,cgb.rc)
# cgb.rc=cgb.rc[!is.na(cgb.rc[,1]),]
# cgb.rc=as.data.frame(cgb.rc)
# # 136030 2
# mean(cgb.rc$start) # 24.8344
# mean(cgb.rc$end) # 26.8344
# 
# # motif 8: GTGTGY - Y:C,T
# gtgtgy.motif=NULL
# for (i in Y){
#   gtgtgy.motif=c(gtgtgy.motif,paste("GTGTG",i,sep=""))
# }
# gtgtgy=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],gtgtgy.motif)
# })
# gtgtgy=lapply(gtgtgy,FUN=function(x) do.call(rbind,x))
# gtgtgy=do.call(rbind,gtgtgy)
# gtgtgy=gtgtgy[!is.na(gtgtgy[,1]),]
# gtgtgy=as.data.frame(gtgtgy)
# # 9380 2
# mean(gtgtgy$start) # 22.64659
# mean(gtgtgy$end) # 27.64659
# 
# # motif 8.rc: Y'CACAC - Y':G,A
# gtgtgy.motif.rc=lapply(gtgtgy.motif,FUN=function(x) {
#   t=chartr("ATCG","TAGC",x)
#   t=reverse(t)
#   t
# })
# gtgtgy.motif.rc=unlist(gtgtgy.motif.rc)
# gtgtgy.rc=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],gtgtgy.motif.rc)
# })
# gtgtgy.rc=lapply(gtgtgy.rc,FUN=function(x) do.call(rbind,x))
# gtgtgy.rc=do.call(rbind,gtgtgy.rc)
# gtgtgy.rc=gtgtgy.rc[!is.na(gtgtgy.rc[,1]),]
# gtgtgy.rc=as.data.frame(gtgtgy.rc)
# # 9486 2
# mean(gtgtgy.rc$start) # 23.35853
# mean(gtgtgy.rc$end) # 28.35853
# 
# # motif 9: GAGAGAGA
# gagagaga.motif=c("GAGAGAGA")
# gagagaga=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],gagagaga.motif)
# })
# gagagaga=lapply(gagagaga,FUN=function(x) do.call(rbind,x))
# gagagaga=do.call(rbind,gagagaga)
# gagagaga=gagagaga[!is.na(gagagaga[,1]),]
# gagagaga=as.data.frame(gagagaga)
# # 890 2
# mean(gagagaga$start) # 23.14382
# mean(gagagaga$end) # 30.14382
# 
# # motif 9.rc: TCTCTCTC
# gagagaga.motif.rc=lapply(gagagaga.motif,FUN=function(x) {
#   t=chartr("ATCG","TAGC",x)
#   t=reverse(t)
#   t
# })
# gagagaga.motif.rc=unlist(gagagaga.motif.rc)
# gagagaga.rc=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],gagagaga.motif.rc)
# })
# gagagaga.rc=lapply(gagagaga.rc,FUN=function(x) do.call(rbind,x))
# gagagaga.rc=do.call(rbind,gagagaga.rc)
# gagagaga.rc=gagagaga.rc[!is.na(gagagaga.rc[,1]),]
# gagagaga.rc=as.data.frame(gagagaga.rc)
# # 900 2
# mean(gagagaga.rc$start) # 21.12222
# mean(gagagaga.rc$end) # 28.12222
# 
# # motif 10: GAYTACA - Y:C,T
# gaytaca.motif=NULL
# for (i in Y){
#   gaytaca.motif=c(gaytaca.motif,paste("GA",i,"TACA",sep=""))
# }
# gaytaca=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],gaytaca.motif)
# })
# gaytaca=lapply(gaytaca,FUN=function(x) do.call(rbind,x))
# gaytaca=do.call(rbind,gaytaca)
# gaytaca=gaytaca[!is.na(gaytaca[,1]),]
# gaytaca=as.data.frame(gaytaca)
# # 4003 2
# mean(gaytaca$start) # 19.84662
# mean(gaytaca$end) # 25.84662
# 
# # motif 10.rc: TGTAY'TC - Y':G,A
# gaytaca.motif.rc=lapply(gaytaca.motif,FUN=function(x) {
#   t=chartr("ATCG","TAGC",x)
#   t=reverse(t)
#   t
# })
# gaytaca.motif.rc=unlist(gaytaca.motif.rc)
# gaytaca.rc=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],gaytaca.motif.rc)
# })
# gaytaca.rc=lapply(gaytaca.rc,FUN=function(x) do.call(rbind,x))
# gaytaca.rc=do.call(rbind,gaytaca.rc)
# gaytaca.rc=gaytaca.rc[!is.na(gaytaca.rc[,1]),]
# gaytaca.rc=as.data.frame(gaytaca.rc)
# # 3203 2
# mean(gaytaca.rc$start) # 25.69216
# mean(gaytaca.rc$end) # 31.69216
# 
# # motif 11: AAAAAWT - W:A,T
# aaaaawt.motif=NULL
# for (i in W){
#   aaaaawt.motif=c(aaaaawt.motif,paste("AAAAA",i,"T",sep=""))
# }
# aaaaawt=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],aaaaawt.motif)
# })
# aaaaawt=lapply(aaaaawt,FUN=function(x) do.call(rbind,x))
# aaaaawt=do.call(rbind,aaaaawt)
# aaaaawt=aaaaawt[!is.na(aaaaawt[,1]),]
# aaaaawt=as.data.frame(aaaaawt)
# # 14316
# mean(aaaaawt$start) # 22.18371
# mean(aaaaawt$end) # 28.18371
# 
# # motif 11.rc: AW'TTTTT - W':T,A
# aaaaawt.motif.rc=lapply(aaaaawt.motif,FUN=function(x) {
#   t=chartr("ATCG","TAGC",x)
#   t=reverse(t)
#   t
# })
# aaaaawt.motif.rc=unlist(aaaaawt.motif.rc)
# aaaaawt.rc=apply(mut,1,function(x) {
#   stri_locate_all_fixed(x[7],aaaaawt.motif.rc)
# })
# aaaaawt.rc=lapply(aaaaawt.rc,FUN=function(x) do.call(rbind,x))
# aaaaawt.rc=do.call(rbind,aaaaawt.rc)
# aaaaawt.rc=aaaaawt.rc[!is.na(aaaaawt.rc[,1]),]
# aaaaawt.rc=as.data.frame(aaaaawt.rc)
# # 13797 2
# mean(aaaaawt.rc$start) # 2392179
# mean(aaaaawt.rc$end) # 29.92179
# 
# Freq.a=table(dncgy$start)
# Freq.a=data.frame(Freq.a,motif="DNCGY")
# colnames(Freq.a)=c("position","freq","motif")
# 
# tab=table(dncgy.rc$start)
# tab=data.frame(tab,motif="DNCGY(RC)")
# colnames(tab)=c("position","freq","motif")
# Freq.a=rbind(Freq.a,tab)
# # 94 3
# 
# Freq.b=table(cgvr$start)
# Freq.b=data.frame(Freq.b,motif="CGVR")
# colnames(Freq.b)=c("position","freq","motif")
# 
# tab=table(cgvr.rc$start)
# tab=data.frame(tab,motif="CGVR(RC)")
# colnames(tab)=c("position","freq","motif")
# Freq.b=rbind(Freq.b,tab)
# # 96 3
# 
# Freq.c=table(ahatatay$start)
# Freq.c=data.frame(Freq.c,motif="AHATATAY")
# colnames(Freq.c)=c("position","freq","motif")
# 
# tab=table(ahatatay.rc$start)
# tab=data.frame(tab,motif="AHATATAY(RC)")
# colnames(tab)=c("position","freq","motif")
# Freq.c=rbind(Freq.c,tab)
# # 88 3
# 
# Freq.d=table(kkaaaaaa$start)
# Freq.d=data.frame(Freq.d,motif="KKAAAAAA")
# colnames(Freq.d)=c("position","freq","motif")
# 
# tab=table(kkaaaaaa.rc$start)
# tab=data.frame(tab,motif="KKAAAAAA(RC)")
# colnames(tab)=c("position","freq","motif")
# Freq.d=rbind(Freq.d,tab)
# # 88 3
# 
# Freq.e=table(acacacrc$start)
# Freq.e=data.frame(Freq.e,motif="ACACACRC")
# colnames(Freq.e)=c("position","freq","motif")
# 
# tab=table(acacacrc.rc$start)
# tab=data.frame(tab,motif="ACACACRC(RC)")
# colnames(tab)=c("position","freq","motif")
# Freq.e=rbind(Freq.e,tab)
# # 88 3
# 
# Freq.f=table(aadttttk$start)
# Freq.f=data.frame(Freq.f,motif="AADTTTTK")
# colnames(Freq.f)=c("position","freq","motif")
# 
# tab=table(aadttttk.rc$start)
# tab=data.frame(tab,motif="AADTTTTK(RC)")
# colnames(tab)=c("position","freq","motif")
# Freq.f=rbind(Freq.f,tab)
# # 88 3
# 
# Freq.g=table(cgb$start)
# Freq.g=data.frame(Freq.g,motif="CGB")
# colnames(Freq.g)=c("position","freq","motif")
# 
# tab=table(cgb.rc$start)
# tab=data.frame(tab,motif="CGB(RC)")
# colnames(tab)=c("position","freq","motif")
# Freq.g=rbind(Freq.g,tab)
# # 98 3
# 
# Freq.h=table(gtgtgy$start)
# Freq.h=data.frame(Freq.h,motif="GTGTGY")
# colnames(Freq.h)=c("position","freq","motif")
# 
# tab=table(gtgtgy.rc$start)
# tab=data.frame(tab,motif="GTGTGY(RC)")
# colnames(tab)=c("position","freq","motif")
# Freq.h=rbind(Freq.h,tab)
# # 92 3
# 
# Freq.i=table(gagagaga$start)
# Freq.i=data.frame(Freq.i,motif="GAGAGAGA")
# colnames(Freq.i)=c("position","freq","motif")
# 
# tab=table(gagagaga.rc$start)
# tab=data.frame(tab,motif="GAGAGAGA(RC)")
# colnames(tab)=c("position","freq","motif")
# Freq.i=rbind(Freq.i,tab)
# # 88 3
# 
# Freq.j=table(gaytaca$start)
# Freq.j=data.frame(Freq.j,motif="GAYTACA")
# colnames(Freq.j)=c("position","freq","motif")
# 
# tab=table(gaytaca.rc$start)
# tab=data.frame(tab,motif="GAYTACA(RC)")
# colnames(tab)=c("position","freq","motif")
# Freq.j=rbind(Freq.j,tab)
# # 90 3
# 
# Freq.k=table(aaaaawt$start)
# Freq.k=data.frame(Freq.k,motif="AAAAAWT")
# colnames(Freq.k)=c("position","freq","motif")
# 
# tab=table(aaaaawt.rc$start)
# tab=data.frame(tab,motif="AAAAAWT(RC)")
# colnames(tab)=c("position","freq","motif")
# Freq.k=rbind(Freq.k,tab)
# 90 3

# enhancer: Freq.a, Freq.b
# 3utr: Freq.a, Freq.b, Freq.c, Freq.d, Freq.e
# 5utr: Freq
# prom: Freq.a, Freq.b, Freq.c, Freq.d, Freq.e, Freq.f
# intron: Freq.a, Freq.b, Freq.c, Freq.d, Freq.e, Freq.f, Freq.g, Freq.h, Freq.i, Freq.j, Freq.k
Freq.a$position=as.numeric(as.character(Freq.a$position))
Freq.a$pos=Freq.a$position-26
Freq.a$logfreq=log(Freq.a$freq)
Freq.a$f.ratio=apply(Freq.a,1,function(x){
  as.numeric(x[2])/mean(Freq.a[which(Freq.a$motif==x[3]),"freq"])
})

Freq.b$position=as.numeric(as.character(Freq.b$position))
Freq.b$pos=Freq.b$position-26
Freq.b$logfreq=log(Freq.b$freq)
Freq.b$f.ratio=apply(Freq.b,1,function(x){
  as.numeric(x[2])/mean(Freq.b[which(Freq.b$motif==x[3]),"freq"])
})

# Freq.c$position=as.numeric(as.character(Freq.c$position))
# Freq.c$pos=Freq.c$position-26
# Freq.c$logfreq=log(Freq.c$freq)
# Freq.c$f.ratio=apply(Freq.c,1,function(x){
#   as.numeric(x[2])/mean(Freq.c[which(Freq.c$motif==x[3]),"freq"])
# })
# 
# Freq.d$position=as.numeric(as.character(Freq.d$position))
# Freq.d$pos=Freq.d$position-26
# Freq.d$logfreq=log(Freq.d$freq)
# Freq.d$f.ratio=apply(Freq.d,1,function(x){
#   as.numeric(x[2])/mean(Freq.d[which(Freq.d$motif==x[3]),"freq"])
# })
# 
# Freq.e$position=as.numeric(as.character(Freq.e$position))
# Freq.e$pos=Freq.e$position-26
# Freq.e$logfreq=log(Freq.e$freq)
# Freq.e$f.ratio=apply(Freq.e,1,function(x){
#   as.numeric(x[2])/mean(Freq.e[which(Freq.e$motif==x[3]),"freq"])
# })
#
# Freq.f$position=as.numeric(as.character(Freq.f$position))
# Freq.f$pos=Freq.f$position-26
# Freq.f$logfreq=log(Freq.f$freq)
# Freq.f$f.ratio=apply(Freq.f,1,function(x){
#   as.numeric(x[2])/mean(Freq.f[which(Freq.f$motif==x[3]),"freq"])
# })
#
# Freq$position=as.numeric(as.character(Freq$position))
# Freq$pos=Freq$position-26
# Freq$logfreq=log(Freq$freq)
# Freq$f.ratio=apply(Freq,1,function(x){
#   as.numeric(x[2])/mean(Freq[which(Freq$motif==x[3]),"freq"])
# })
# 
# Freq.g$position=as.numeric(as.character(Freq.g$position))
# Freq.g$pos=Freq.g$position-26
# Freq.g$logfreq=log(Freq.g$freq)
# Freq.g$f.ratio=apply(Freq.g,1,function(x){
#   as.numeric(x[2])/mean(Freq.g[which(Freq.g$motif==x[3]),"freq"])
# })
# 
# Freq.h$position=as.numeric(as.character(Freq.h$position))
# Freq.h$pos=Freq.h$position-26
# Freq.h$logfreq=log(Freq.h$freq)
# Freq.h$f.ratio=apply(Freq.h,1,function(x){
#   as.numeric(x[2])/mean(Freq.h[which(Freq.h$motif==x[3]),"freq"])
# })
# 
# Freq.i$position=as.numeric(as.character(Freq.i$position))
# Freq.i$pos=Freq.i$position-26
# Freq.i$logfreq=log(Freq.i$freq)
# Freq.i$f.ratio=apply(Freq.i,1,function(x){
#   as.numeric(x[2])/mean(Freq.i[which(Freq.i$motif==x[3]),"freq"])
# })
# 
# Freq.j$position=as.numeric(as.character(Freq.j$position))
# Freq.j$pos=Freq.j$position-26
# Freq.j$logfreq=log(Freq.j$freq)
# Freq.j$f.ratio=apply(Freq.j,1,function(x){
#   as.numeric(x[2])/mean(Freq.j[which(Freq.j$motif==x[3]),"freq"])
# })
# 
# Freq.k$position=as.numeric(as.character(Freq.k$position))
# Freq.k$pos=Freq.k$position-26
# Freq.k$logfreq=log(Freq.k$freq)
# Freq.k$f.ratio=apply(Freq.k,1,function(x){
#   as.numeric(x[2])/mean(Freq.k[which(Freq.k$motif==x[3]),"freq"])
# })

# positional figures
ggplot(data=Freq.a, aes(x=pos, y=f.ratio, group=motif,color=motif)) +
# ggplot(data=Freq.b, aes(x=pos, y=f.ratio, group=motif,color=motif)) +
# ggplot(data=Freq.c, aes(x=pos, y=f.ratio, group=motif,color=motif)) +
# ggplot(data=Freq.d, aes(x=pos, y=f.ratio, group=motif,color=motif)) +
# ggplot(data=Freq.e, aes(x=pos, y=f.ratio, group=motif,color=motif)) +
# ggplot(data=Freq.f, aes(x=pos, y=f.ratio, group=motif,color=motif)) +
# ggplot(data=Freq.g, aes(x=pos, y=f.ratio, group=motif,color=motif)) +
# ggplot(data=Freq.h, aes(x=pos, y=f.ratio, group=motif,color=motif)) +
# ggplot(data=Freq.i, aes(x=pos, y=f.ratio, group=motif,color=motif)) +
# ggplot(data=Freq.j, aes(x=pos, y=f.ratio, group=motif,color=motif)) +
# ggplot(data=Freq.k, aes(x=pos, y=f.ratio, group=motif,color=motif)) +
# ggplot(data=Freq, aes(x=pos, y=f.raito, group=motif,color=motif)) +
  geom_line() +
  geom_point() +
  ylab(paste("freq/mean freq for that motif in mutated","sequences",sep="\n"))+
  ggtitle(paste("Frequency of motif enriched in mutated sequences","w.r.t. position in mutated sequences",sep="\n")) +
  geom_vline(xintercept=0)+
  theme(text=element_text(size=20))

save.image(file="dreme_final_enh.RData")
# save.image(file="dreme_final_3utr.RData")
# save.image(file="dreme_final_5utr.RData")
# save.image(file="dreme_final_prom.RData")
# save.image(file="dreme_final_intron.RData")

#########################################################################################################################

###
# lasso selection for nucleotide contexts and enriched motifs
###
dnaSeq$dna=as.character(dnaSeq$dna)
dnaSeq$fiveMer=substr(dnaSeq$dna,24,28)
dnaSeq$threeMer=substr(dnaSeq$fiveMer,2,4)
dnaSeq$oneMer=substr(dnaSeq$threeMer,2,2) 

length(unique(dnaSeq$fiveMer)) 
length(unique(dnaSeq$threeMer))
length(unique(dnaSeq$oneMer)) 

# enhancers
###
# get presence of motifs
# include only enriched
# BRCG, BCG
### 
# motifs=c("brcg","bcg")
# motif 1: BRCG - B:C,G,T; R:A,G (window 3)
a=3
dnaSeq$brcg=sapply(dnaSeq$dna,FUN=function(x){
  if(sum(str_detect(substr(x,(26-a),(26+a)),unique(c(brcg.motif,brcg.motif.rc))))>=1)
  {1} else {0}
})

# motif 2: BCG - B:C,G,T (window 2)
a=2
dnaSeq$bcg=sapply(dnaSeq$dna,FUN=function(x){
  if(sum(str_detect(substr(x,(26-a),(26+a)),unique(c(bcg.motif,bcg.motif.rc))))>=1)
  {1} else {0}
})

# 3UTR
###
# get presence of motifs
# include only enriched
# CGBH, ACACACRC, CGH, KAAAAAAA, ATATATAT
###
# motifs=c("cgbh","acacacrc","cgh","kaaaaaaa","atatatat")
# motif 1: CGBH - B:C,G,T; H:A,C,T (window 2)
# a=2
# dnaSeq$cgbh=apply(dnaSeq,1,function(x){
#   if(sum(str_detect(substr(x[7],26-a,26+a),unique(c(cgbh.motif,cgbh.motif.rc))))>=1)
#   {1} else {0}
# })
# 
# # motif 2: ACACACRC - R:A,G (window 25)
# dnaSeq$acacacrc=apply(dnaSeq,1,function(x){
#   if(sum(str_detect(x[7],unique(c(acacacrc.motif,acacacrc.motif.rc))))>=1)
#   {1} else {0}
# })
# 
# # motif 3: CGH - H:A,C,T (window 2)
# a=2
# dnaSeq$cgh=apply(dnaSeq,1,function(x){
#   if(sum(str_detect(substr(x[7],26-a,26+a),unique(c(cgh.motif,cgh.motif.rc))))>=1)
#   {1} else {0}
# })
# 
# # motif 4: KAAAAAAA - K:G,T (window 10)
# a=10
# dnaSeq$kaaaaaaa=apply(dnaSeq,1,function(x){
#   if(sum(str_detect(substr(x[7],26-a,26+a),unique(c(kaaaaaaa.motif,kaaaaaaa.motif.rc))))>=1)
#   {1} else {0}
# })
# 
# # motif 5: ATATATAT (window 25)
# dnaSeq$atatatat=apply(dnaSeq,1,function(x){
#   if(sum(str_detect(x[7],unique(c(atatatat.motif))))>=1)
#   {1} else {0}
# })

# 5UTR
###
# get presence of motifs
# include only enriched
# CGH
###
# motifs=c("cgh")
# motif 1: CGH - H:A,C,T (window 2)
# a=2
# dnaSeq$cgh=apply(dnaSeq,1,function(x){
#   if(sum(str_detect(substr(x[7],26-a,26+a),unique(c(cgh.motif,cgh.motif.rc))))>=1)
#   {1} else {0}
# })

# promoter
###
# get presence of motifs
# include only enriched
# DCGY, CGD, AAAAAAWT, AYACACAC, GGGGAAAA, ATATATAT
###
# motifs=c("dcgy","cgd","aaaaaawt","ayacacac","ggggaaaa","atatatat")
# motif 1: DCGY - D:A,G,T; Y:C,T (window 2)
# a=2
# dnaSeq$dcgy=apply(dnaSeq,1,function(x){
#   if(sum(str_detect(substr(x[7],26-a,26+a),unique(c(dcgy.motif,dcgy.motif.rc))))>=1)
#   {1} else {0}
# })
# 
# # motif 2: CGD - D:A,G,T (window 2)
# a=2
# dnaSeq$cgd=apply(dnaSeq,1,function(x){
#   if(sum(str_detect(substr(x[7],26-a,26+a),unique(c(cgd.motif,cgd.motif.rc))))>=1)
#   {1} else {0}
# })
# 
# # motif 3: AAAAAAWT - W:A,T (window 10)
# a=10
# dnaSeq$aaaaaawt=apply(dnaSeq,1,function(x){
#   if(sum(str_detect(substr(x[7],26-a,26+a),unique(c(aaaaaawt.motif,aaaaaawt.motif.rc))))>=1)
#   {1} else {0}
# })
# 
# # motif 4: AYACACAC - Y:C,T (window 25)
# dnaSeq$ayacacac=apply(dnaSeq,1,function(x){
#   if(sum(str_detect(x[7],unique(c(ayacacac.motif,ayacacac.motif.rc))))>=1)
#   {1} else {0}
# })
# 
# # motif 5: GGGGAAAA (window 5)
# a=5
# dnaSeq$ggggaaaa=apply(dnaSeq,1,function(x){
#   if(sum(str_detect(substr(x[7],26-a,26+a),unique(c(ggggaaaa.motif,ggggaaaa.motif.rc))))>=1)
#   {1} else {0}
# })
# 
# # motif 6: ATATATAT (window 25)
# dnaSeq$atatatat=apply(dnaSeq,1,function(x){
#   if(sum(str_detect(x[7],unique(c(atatatat.motif))))>=1)
#   {1} else {0}
# })

# intron
###
# get presence of motifs
# include both enriched and depleted
# DNCGY, CGVR, AHATATAY, KKAAAAAA, ACACACRC, AADTTTTK, CGB, GTGTGY, GAGAGAGA, GAYTACA, AAAAAWT
###
# motifs=c("dncgy","cgvr","ahatatay","kkaaaaaa","acacacrc","aadttttk","cgb","gtgtgy","gagagaga","gaytaca","aaaaawt")
# motif 1: DNCGY - D:A,G,T; N:A,C,G,T; Y:C,T (window 3)
# a=3
# dnaSeq$dncgy=apply(dnaSeq,1,function(x){
#   if(sum(str_detect(substr(x[7],26-a,26+a),unique(c(dncgy.motif,dncgy.motif.rc))))>=1)
#   {1} else {0}
# })
# 
# # motif 2: CGVR - V:A,C,G; R:A,G (window 3)
# a=3
# dnaSeq$cgvr=apply(dnaSeq,1,function(x){
#   if(sum(str_detect(substr(x[7],26-a,26+a),unique(c(cgvr.motif,cgvr.motif.rc))))>=1)
#   {1} else {0}
# })
# 
# # motif 3: AHATATAY - H:A,C,T; Y:C,T (window 25)
# dnaSeq$ahatatay=apply(dnaSeq,1,function(x){
#   if(sum(str_detect(x[7],unique(c(ahatatay.motif,ahatatay.motif.rc))))>=1)
#   {1} else {0}
# })
# 
# # motif 4: KKAAAAAA - K:G,T (window 10)
# a=10
# dnaSeq$kkaaaaaa=apply(dnaSeq,1,function(x){
#   if(sum(str_detect(substr(x[7],a,b),unique(c(kkaaaaaa.motif,kkaaaaaa.motif.rc))))>=1)
#   {1} else {0}
# })
# 
# # motif 5: ACACACRC - R:A,G (window 25)
# dnaSeq$acacacrc=apply(dnaSeq,1,function(x){
#   if(sum(str_detect(x[7],unique(c(acacacrc.motif,acacacrc.motif.rc))))>=1)
#   {1} else {0}
# })
# 
# # motif 6: AADTTTTK - D:A,G,T; K:G,T (window 10)
# a=10
# dnaSeq$aadttttk=apply(dnaSeq,1,function(x){
#   if(sum(str_detect(substr(x[7],26-a,26+a),unique(c(aadttttk.motif,aadttttk.motif.rc))))>=1)
#   {1} else {0}
# })
# 
# # motif 7: CGB - B:C,G,T (window 2)
# a=2
# dnaSeq$cgb=apply(dnaSeq,1,function(x){
#   if(sum(str_detect(substr(x[7],26-a,26+a),unique(c(cgb.motif,cgb.motif.rc))))>=1)
#   {1} else {0}
# })
# 
# # motif 8: GTGTGY - Y:C,T (window 25)
# dnaSeq$gtgtgy=apply(dnaSeq,1,function(x){
#   if(sum(str_detect(x[7],unique(c(gtgtgy.motif,gtgtgy.motif.rc))))>=1)
#   {1} else {0}
# })
# 
# # motif 9: GAGAGAGA (window 25)
# dnaSeq$gagagaga=apply(dnaSeq,1,function(x){
#   if(sum(str_detect(x[7],unique(c(gagagaga.motif,gagagaga.motif.rc))))>=1)
#   {1} else {0}
# })
# 
# # motif 10: GAYTACA - Y:C,T (window 25)
# dnaSeq$gaytaca=apply(dnaSeq,1,function(x){
#   if(sum(str_detect(x[7],unique(c(gaytaca.motif,gaytaca.motif.rc))))>=1)
#   {1} else {0}
# })
# 
# # motif 11: AAAAAWT - W:A,T (window 10)
# a=10
# dnaSeq$aaaaawt=apply(dnaSeq,1,function(x){
#   if(sum(str_detect(substr(x[7],26-a,26+a),unique(c(aaaaawt.motif,aaaaawt.motif.rc))))>=1)
#   {1} else {0}
# })

# pair up with rev comp
# rev comp for 5mer
dnaSeq$five.rc=sapply(dnaSeq$fiveMer,FUN=function(x) {
  t=reverse(chartr("ATCG","TAGC",x))
}) 
# rev comp for 3mer
dnaSeq$three.rc=sapply(dnaSeq$threeMer,FUN=function(x) {
  t=reverse(chartr("ATCG","TAGC",x))
}) 
# rev comp for 1mer
dnaSeq$one.rc=sapply(dnaSeq$oneMer,FUN=function(x){
  t=reverse(chartr("ATCG","TAGC",x))
}) 

# 5mer and 5mer.rc pair
# 512 unique five.pair
dnaSeq$five.pair=ifelse(substr(dnaSeq$fiveMer,3,3) %in% c("A","G"),dnaSeq$fiveMer,dnaSeq$five.rc)
# 3mer and 3mer.rc pair
# 32 unique three.pair
dnaSeq$three.pair=ifelse(substr(dnaSeq$threeMer,2,2) %in% c("A","G"),dnaSeq$threeMer,dnaSeq$three.rc)
# 1mer and 1mer.rc pair
# 2 unique one.pair
dnaSeq$one.pair=ifelse(dnaSeq$oneMer %in% c("A","G"),dnaSeq$oneMer,dnaSeq$one.rc)

# extract flanks for 3mer and 5mer
# right/3' flank of 3mer
# 8 unique right flank of 3mer
dnaSeq$three.right=sapply(dnaSeq$three.pair,FUN=function(x) {
  substr(x,2,3)
}) 
# left/5' flank of 3mer
# 8 unique left flank of 3mer
dnaSeq$three.left=sapply(dnaSeq$three.pair,FUN=function(x) {
  substr(x,1,2)
}) 
# right/3' flank of 5mer
# 32 unique right flank of 5mer
dnaSeq$five.right=sapply(dnaSeq$five.pair,FUN=function(x) {
  substr(x,3,5)
})
# left/5' flank of 5mer
# 32 unique left flank of 5mer
dnaSeq$five.left=sapply(dnaSeq$five.pair,FUN=function(x) {
  substr(x,1,3)
})

tab2=aggregate(mut~one.pair+three.pair+three.right+three.left+five.pair+five.right+five.left+brcg+bcg,dnaSeq,FUN=function(x){y=sum(x==1);z=sum(x==0);return(cbind(y,z))}) # enhancer
# tab2=aggregate(mut~one.pair+three.pair+three.right+three.left+five.pair+five.right+five.left+cgbh+acacacrc+cgh+kaaaaaaa+atatatat,dnaSeq,FUN=function(x){y=sum(x==1);z=sum(x==0);return(cbind(y,z))}) # 3UTR
# tab2=aggregate(mut~one.pair+three.pair+three.right+three.left+five.pair+five.right+five.left+cgh,dnaSeq,FUN=function(x){y=sum(x==1);z=sum(x==0);return(cbind(y,z))})  # 5UTR
# tab2=aggregate(mut~one.pair+three.pair+three.right+three.left+five.pair+five.right+five.left+dcgy+cgd+aaaaaawt+ayacacac+ggggaaaa+atatatat,dnaSeq,FUN=function(x){y=sum(x==1);z=sum(x==0);return(cbind(y,z))}) # promoter
# tab2=aggregate(mut~one.pair+three.pair+three.right+three.left+five.pair+five.right+five.left+dncgy+cgvr+ahatatay+kkaaaaaa+acacacrc+aadttttk+cgb+gtgtgy+gagagaga+gaytaca+aaaaawt,dnaSeq,FUN=function(x){y=sum(x==1);z=sum(x==0);return(cbind(y,z))}) # intron
nrow(unique(dnaSeq[,c(motifs,"one.pair","three.pair","three.right","three.left","five.pair","five.right","five.left")])) 
dim(tab2) 
df=tab2$mut # success(mutated)-1, failure(non-mutated)-0
colnames(df)=c("Success","Failure") 
df=data.frame(df,tab2[,1:(ncol(tab2)-1)]) 
nrow(unique(df[,3:9])) # 512
dim(df)
df=df[,c(3:ncol(df),1:2)]

# lasso logistic 
one=as.factor(df$one.pair) # 2
three=as.factor(df$three.pair) # 32
three.right=as.factor(df$three.right) # 8
three.left=as.factor(df$three.left) # 8
five=as.factor(df$five.pair) # 512
five.right=as.factor(df$five.right) # 32
five.left=as.factor(df$five.left) # 32

# convert to dummy variables with all factor levels
xfactors=cbind(model.matrix(~one-1),model.matrix(~three-1),model.matrix(~three.right-1),model.matrix(~three.left-1),model.matrix(~five-1),model.matrix(~five.right-1),model.matrix(~five.left-1),df[,motifs])
x <- as.matrix(data.frame(xfactors))

y_data=cbind(df$Failure,df$Success)
colnames(y_data)=c("Failure","Success")

### run logistic lasso with intercept, all levels (code5.R, stabs2)
# 4.08pm-4.08pm enhancers
# 5.48pm-5.48pm 3UTR
# 5.11pm-5.11pm 5UTR
# 1.58pm-1.58pm prom
# 3.19pm-3.22pm intron

length(stabs2$selected.min)
coef.min <- stabs2$coef.min[,stabs2$selected.min] 
coef.min.sign=data.frame(feature = colnames(coef.min),sign=numeric(ncol(coef.min)))
# t.test  to determine sign of features selected
for (i in 1:ncol(coef.min))
{
  if (t.test(coef.min[,i],alternative = "less")$p.value <= 0.05)
  {
    coef.min.sign$sign[i] = -1
  } else {coef.min.sign$sign[i] = 1}
  print(colnames(coef.min)[i])
}
feature=names(stabs2$freq.min)
f=stabs2$freq.min
freq=data.frame(feature,f)
sel=merge(coef.min.sign,freq,by="feature")
sel[order(sel$f,decreasing=TRUE),]
# enhancer: 26
# feature sign    f
# 14      threeAAG    1 1.00
# 20  three.leftCG    1 1.00
# 22 three.rightAA   -1 1.00
# 10 five.rightGAG   -1 0.99
# 24 three.rightGC    1 0.99
# 19  three.leftCA   -1 0.97
# 5   five.leftTCA    1 0.96
# 7  five.rightAAG   -1 0.96
# 3   five.leftGTA    1 0.94
# 11 five.rightGCA    1 0.94
# 25      threeTAG   -1 0.94
# 15      threeAGA    1 0.92
# 23 three.rightAC    1 0.92
# 13 five.rightGGC   -1 0.91
# 21  three.leftGA   -1 0.91
# 18      threeGGC    1 0.90
# 4   five.leftTAA    1 0.89
# 16      threeAGG    1 0.88
# 9  five.rightAGT    1 0.87
# 6   five.leftTGG    1 0.86
# 8  five.rightACC    1 0.86
# 17      threeCGT    1 0.85
# 2   five.leftCAG    1 0.81
# 1   five.leftAGA   -1 0.80
# 26      threeTGG   -1 0.80
# 12 five.rightGCC    1 0.76

# 3UTR: 30
# feature sign    f
# 15      kaaaaaaa    1 1.00
# 22  three.leftCA   -1 1.00
# 23  three.leftCG    1 1.00
# 27 three.rightAC    1 1.00
# 28 three.rightGC    1 1.00
# 17      threeAAG    1 0.99
# 26 three.rightAA   -1 0.99
# 1       acacacrc    1 0.98
# 3   five.leftAGA   -1 0.98
# 12 five.rightAGT    1 0.98
# 19      threeAGG    1 0.97
# 11 five.rightACC    1 0.96
# 13 five.rightGCA    1 0.96
# 6   five.leftTCA    1 0.92
# 20      threeAGT   -1 0.92
# 5   five.leftGTA    1 0.91
# 7   five.leftTGG    1 0.91
# 21      threeCGT    1 0.91
# 2       atatatat    1 0.90
# 4   five.leftGAA   -1 0.89
# 9  five.rightAAG   -1 0.89
# 16      threeAAA   -1 0.89
# 30      threeTGG   -1 0.84
# 24  three.leftGA   -1 0.83
# 25  three.leftTG   -1 0.82
# 29      threeTAC    1 0.81
# 14     fiveTTAAA    1 0.78
# 8  five.rightAAC   -1 0.77
# 10 five.rightAAT    1 0.76
# 18      threeAGA    1 0.76

# 5UTR: 17 
# feature sign    f
# 13  three.leftCG    1 1.00
# 4   five.leftTGG    1 0.99
# 12      threeCGT    1 0.99
# 17 three.rightGC    1 0.99
# 10      threeAAG    1 0.98
# 15 three.rightAA   -1 0.98
# 7  five.rightGCA    1 0.95
# 11      threeAGG    1 0.95
# 16 three.rightAC    1 0.94
# 6  five.rightAGT    1 0.93
# 14  three.leftGA   -1 0.90
# 9  five.rightGGC   -1 0.87
# 8  five.rightGCC    1 0.85
# 1   five.leftACG    1 0.83
# 3   five.leftTCA    1 0.81
# 2   five.leftAGA   -1 0.80
# 5  five.rightAAG   -1 0.78

# promoter: 33
# feature sign    f
# 1       aaaaaawt    1 1.00
# 2       ayacacac    1 1.00
# 8   five.leftTCA    1 1.00
# 9   five.leftTGG    1 1.00
# 14 five.rightGCA    1 1.00
# 20      threeAAG    1 1.00
# 22      threeAGG    1 1.00
# 23      threeCGT    1 1.00
# 25  three.leftCA   -1 1.00
# 26  three.leftCG    1 1.00
# 28 three.rightAA   -1 1.00
# 29 three.rightAC    1 1.00
# 12 five.rightAGT    1 0.99
# 18 five.rightGGT    1 0.99
# 30 three.rightGC    1 0.99
# 5   five.leftAGA   -1 0.98
# 21      threeAGA    1 0.97
# 33      threeTGG   -1 0.97
# 6   five.leftGTA    1 0.96
# 13 five.rightGAG   -1 0.96
# 24      threeGGC    1 0.96
# 11 five.rightAAG   -1 0.95
# 15 five.rightGCC    1 0.95
# 27  three.leftGA   -1 0.94
# 10 five.rightAAC   -1 0.92
# 16 five.rightGCG   -1 0.91
# 31      threeTAC    1 0.84
# 7   five.leftTAA    1 0.82
# 17 five.rightGGC   -1 0.81
# 32      threeTAG   -1 0.80
# 4   five.leftACG    1 0.79
# 3           dcgy    1 0.77
# 19     fiveTTAAA    1 0.76

# intron: 64
# feature sign    f
# 1        aaaaawt    1 1.00
# 3       acacacrc    1 1.00
# 4       ahatatay    1 1.00
# 5          dncgy    1 1.00
# 9   five.leftAGA   -1 1.00
# 18  five.leftGTA    1 1.00
# 21  five.leftTCA    1 1.00
# 23  five.leftTGG    1 1.00
# 26 five.rightAAG   -1 1.00
# 29 five.rightAGT    1 1.00
# 30 five.rightATA   -1 1.00
# 31 five.rightGAG   -1 1.00
# 32 five.rightGCA    1 1.00
# 33 five.rightGCC    1 1.00
# 41      kkaaaaaa    1 1.00
# 45      threeAAG    1 1.00
# 46      threeAGA    1 1.00
# 47      threeAGG    1 1.00
# 48      threeAGT   -1 1.00
# 51      threeCGT    1 1.00
# 54  three.leftCA   -1 1.00
# 55  three.leftCG    1 1.00
# 57 three.rightAA   -1 1.00
# 58 three.rightAC    1 1.00
# 59 three.rightGC    1 1.00
# 64      threeTGG   -1 1.00
# 2       aadttttk    1 0.99
# 28 five.rightACC    1 0.99
# 40        gtgtgy    1 0.99
# 11  five.leftATA    1 0.98
# 35 five.rightGGT    1 0.98
# 37     fiveTTAAA    1 0.98
# 39       gaytaca   -1 0.98
# 61      threeTAC    1 0.98
# 53      threeGGC    1 0.97
# 62      threeTAG   -1 0.97
# 38      gagagaga    1 0.96
# 63      threeTGA    1 0.96
# 10  five.leftAGG   -1 0.95
# 14  five.leftCTG   -1 0.95
# 19  five.leftTAA    1 0.95
# 42          oneA   -1 0.95
# 13  five.leftCAA    1 0.94
# 34 five.rightGGC   -1 0.94
# 7      fiveGGAAA    1 0.93
# 25 five.rightAAC   -1 0.93
# 8   five.leftACA   -1 0.92
# 22  five.leftTGA    1 0.92
# 16  five.leftGAG   -1 0.89
# 56  three.leftGA   -1 0.88
# 43          oneG    1 0.87
# 15  five.leftGAA   -1 0.86
# 27 five.rightAAT    1 0.85
# 49      threeCAC   -1 0.85
# 6      fiveAAATT    1 0.83
# 24  five.leftTTA   -1 0.83
# 36 five.rightGTT   -1 0.82
# 44      threeAAA   -1 0.81
# 20  five.leftTAG    1 0.79
# 12  five.leftATG    1 0.78
# 17  five.leftGCA   -1 0.78
# 50      threeCGA   -1 0.76
# 52      threeGAG   -1 0.76
# 60      threeTAA    1 0.75

length(stabs2$selected.1se) 
coef.1se <- stabs2$coef.1se[,stabs2$selected.1se] 
coef.1se.sign=data.frame(feature = colnames(coef.1se),sign=numeric(ncol(coef.1se)))
# t.test  to determine sign of features selected
for (i in 1:ncol(coef.1se))
{
  if (t.test(coef.1se[,i],alternative = "less")$p.value <= 0.05)
  {
    coef.1se.sign$sign[i] = -1
  } else {coef.1se.sign$sign[i] = 1}
  print(colnames(coef.1se)[i])
}
feature=names(stabs2$freq.1se)
f=stabs2$freq.1se
freq=data.frame(feature,f)
sel=merge(coef.1se.sign,freq,by="feature") 
sel[order(sel$f,decreasing=TRUE),]
# enhancer: 8
# feature sign    f
# 4  three.leftCG    1 1.00
# 6 three.rightAA   -1 1.00
# 8 three.rightGC    1 1.00
# 2      threeAAG    1 0.98
# 7 three.rightAC    1 0.94
# 3  three.leftCA   -1 0.93
# 5  three.leftGA   -1 0.92
# 1 five.rightGCA    1 0.76

# 3UTR: 7
# feature sign    f
# 4  three.leftCG    1 1.00
# 7 three.rightGC    1 0.98
# 5 three.rightAA   -1 0.97
# 3  three.leftCA   -1 0.96
# 6 three.rightAC    1 0.88
# 1      kaaaaaaa    1 0.86
# 2      threeAAG    1 0.85

# 5UTR: 8
# feature sign    f
# 4  three.leftCG    1 1.00
# 6 three.rightAA   -1 1.00
# 8 three.rightGC    1 1.00
# 5  three.leftGA   -1 0.93
# 7 three.rightAC    1 0.85
# 2      threeAAG    1 0.83
# 3      threeCGT    1 0.78
# 1 five.rightGCA    1 0.75

# promoter: 20
# feature sign    f
# 9  five.rightGCA    1 1.00
# 10      threeAAG    1 1.00
# 11      threeAGG    1 1.00
# 14  three.leftCA   -1 1.00
# 15  three.leftCG    1 1.00
# 18 three.rightAC    1 1.00
# 19 three.rightGC    1 1.00
# 12      threeCGT    1 0.99
# 16  three.leftGA   -1 0.97
# 17 three.rightAA   -1 0.97
# 20      threeTGG   -1 0.96
# 4   five.leftAGA   -1 0.95
# 1       aaaaaawt    1 0.91
# 8  five.rightAGT    1 0.91
# 6   five.leftTGG    1 0.89
# 7  five.rightAAG   -1 0.87
# 5   five.leftTCA    1 0.83
# 2       ayacacac    1 0.81
# 3           dcgy    1 0.76
# 13      threeGGC    1 0.75

# intron: 35
# feature sign    f
# 1        aaaaawt    1 1.00
# 4       ahatatay    1 1.00
# 5          dncgy    1 1.00
# 15 five.rightGCA    1 1.00
# 22      threeAAG    1 1.00
# 23      threeAGG    1 1.00
# 27  three.leftCA   -1 1.00
# 28  three.leftCG    1 1.00
# 30 three.rightAA   -1 1.00
# 31 three.rightAC    1 1.00
# 32 three.rightGC    1 1.00
# 35      threeTGG   -1 1.00
# 6   five.leftAGA   -1 0.99
# 18      kkaaaaaa    1 0.99
# 17        gtgtgy    1 0.98
# 9   five.leftTCA    1 0.97
# 11 five.rightAAG   -1 0.97
# 24      threeAGT   -1 0.96
# 29  three.leftGA   -1 0.96
# 2       aadttttk    1 0.94
# 21      threeAAA   -1 0.94
# 13 five.rightAGT    1 0.93
# 10  five.leftTGG    1 0.92
# 3       acacacrc    1 0.87
# 25      threeCGT    1 0.86
# 8   five.leftGTA    1 0.85
# 12 five.rightACC    1 0.85
# 16 five.rightGCC    1 0.85
# 14 five.rightATA   -1 0.83
# 26      threeGGC    1 0.81
# 19          oneA   -1 0.79
# 33      threeTAC    1 0.79
# 7   five.leftAGG   -1 0.77
# 34      threeTAG   -1 0.77
# 20          oneG    1 0.76

length(stabs2$selected.1st) 
coef.1st <- stabs2$coef.1st[,stabs2$selected.1st] 
coef.1st.sign=data.frame(feature = colnames(coef.1st),sign=numeric(ncol(coef.1st)))
# t.test  to determine sign of features selected
for (i in 1:ncol(coef.1st))
{
  if (t.test(coef.1st[,i],alternative = "less")$p.value <= 0.05)
  {
    coef.1st.sign$sign[i] = -1
  } else {coef.1st.sign$sign[i] = 1}
  print(colnames(coef.1st)[i])
}
feature=names(stabs2$freq.1st)
f=stabs2$freq.1st
freq=data.frame(feature,f)
sel=merge(coef.1st.sign,freq,by="feature")
sel[order(sel$f,decreasing=TRUE),]
# enhancer: 7
# feature sign    f
# 3  three.leftCG    1 1.00
# 5 three.rightAA   -1 1.00
# 7 three.rightGC    1 1.00
# 1      threeAAG    1 0.96
# 6 three.rightAC    1 0.96
# 4  three.leftGA   -1 0.93
# 2  three.leftCA   -1 0.92

# 3UTR: 10
# feature sign    f
# 5   three.leftCG    1 1.00
# 9  three.rightAC    1 1.00
# 10 three.rightGC    1 1.00
# 4   three.leftCA   -1 0.99
# 8  three.rightAA   -1 0.99
# 3       threeAAG    1 0.96
# 2       kaaaaaaa    1 0.85
# 1  five.rightAGT    1 0.82
# 7   three.leftTG   -1 0.82
# 6   three.leftGA   -1 0.77

# 5UTR: 10
# feature sign    f
# 6   three.leftCG    1 1.00
# 8  three.rightAA   -1 1.00
# 10 three.rightGC    1 1.00
# 7   three.leftGA   -1 0.98
# 9  three.rightAC    1 0.90
# 4       threeAAG    1 0.81
# 1   five.leftTGG    1 0.80
# 3  five.rightGCA    1 0.80
# 2  five.rightAGT    1 0.76
# 5       threeCGT    1 0.76

# promoter: 10
# feature sign    f
# 6   three.leftCG    1 1.00
# 10 three.rightGC    1 1.00
# 5   three.leftCA   -1 0.99
# 9  three.rightAC    1 0.99
# 7   three.leftGA   -1 0.98
# 8  three.rightAA   -1 0.98
# 1  five.rightGCA    1 0.97
# 2       threeAAG    1 0.97
# 3       threeAGG    1 0.85
# 4       threeCGT    1 0.77

# intron: 9
# feature sign    f
# 4  three.leftCA   -1 1.00
# 5  three.leftCG    1 1.00
# 7 three.rightAA   -1 1.00
# 8 three.rightAC    1 1.00
# 9 three.rightGC    1 1.00
# 3      threeAAG    1 0.98
# 6  three.leftGA   -1 0.97
# 1       aaaaawt    1 0.93
# 2      kkaaaaaa    1 0.81

###
# fit logistic glm
###
# using features selected by lasso.1se
sel=names(stabs2$selected.1se) 
sel=x[,sel]
logr=glm(y_data[,c("Success","Failure")]~sel,family=binomial(link="logit"))
# enhancer
# Call:
#   glm(formula = y_data[, c("Success", "Failure")] ~ sel, family = binomial(link = "logit"))
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -3.6430  -1.0642  -0.2032   0.7113   6.1346  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)      -0.57098    0.01914 -29.832  < 2e-16 ***
#   selthreeAAG       0.65799    0.05897  11.157  < 2e-16 ***
#   selthree.rightAA -0.31223    0.04168  -7.491 6.82e-14 ***
#   selthree.rightAC  0.29875    0.04270   6.996 2.63e-12 ***
#   selthree.rightGC  0.28858    0.04595   6.280 3.38e-10 ***
#   selthree.leftCA  -0.19124    0.04056  -4.715 2.42e-06 ***
#   selthree.leftCG   1.92508    0.05155  37.341  < 2e-16 ***
#   selthree.leftGA  -0.23536    0.04457  -5.281 1.29e-07 ***
#   selfive.rightGCA  0.29745    0.07321   4.063 4.84e-05 ***
#   ---
#   Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 3871.3  on 721  degrees of freedom
# Residual deviance: 1285.7  on 713  degrees of freedom
# AIC: 3601.9
# 
# Number of Fisher Scoring iterations: 4

# 3UTR
# Call:
#   glm(formula = y_data[, c("Success", "Failure")] ~ sel, family = binomial(link = "logit"))
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -7.2954  -1.1960  -0.1894   1.0560  13.3818  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)      -0.56764    0.01100 -51.593   <2e-16 ***
#   selthreeAAG       0.59798    0.03439  17.386   <2e-16 ***
#   selthree.rightAA -0.27638    0.02282 -12.109   <2e-16 ***
#   selthree.rightAC  0.34629    0.02545  13.605   <2e-16 ***
#   selthree.rightGC  0.38610    0.02511  15.377   <2e-16 ***
#   selthree.leftCA  -0.23117    0.02448  -9.442   <2e-16 ***
#   selthree.leftCG   2.20284    0.03498  62.967   <2e-16 ***
#   selkaaaaaaa       1.13654    0.05785  19.646   <2e-16 ***
#   ---
#   Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 10982.3  on 1211  degrees of freedom
# Residual deviance:  3396.4  on 1204  degrees of freedom
# AIC: 6336.6
# 
# Number of Fisher Scoring iterations: 4

# 5UTR
# Call:
#   glm(formula = y_data[, c("Success", "Failure")] ~ sel, family = binomial(link = "logit"))
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -4.0212  -0.9127  -0.0429   0.7684   4.1366  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)      -0.64501    0.02131 -30.265  < 2e-16 ***
#   selthreeAAG       0.71188    0.07484   9.512  < 2e-16 ***
#   selthreeCGT       0.55390    0.09745   5.684 1.31e-08 ***
#   selthree.rightAA -0.40733    0.05731  -7.108 1.18e-12 ***
#   selthree.rightAC  0.28109    0.05321   5.282 1.27e-07 ***
#   selthree.rightGC  0.24674    0.04170   5.917 3.28e-09 ***
#   selthree.leftCG   1.24856    0.04103  30.432  < 2e-16 ***
#   selthree.leftGA  -0.30828    0.05301  -5.815 6.05e-09 ***
#   selfive.rightGCA  0.35176    0.07567   4.649 3.34e-06 ***
#   ---
#   Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 2923.61  on 511  degrees of freedom
# Residual deviance:  944.67  on 503  degrees of freedom
# AIC: 2883.4
# 
# Number of Fisher Scoring iterations: 4

# promoter
# Call:
#   glm(formula = y_data[, c("Success", "Failure")] ~ sel, family = binomial(link = "logit"))
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -6.2972  -1.2085  -0.2970   0.9789  11.4394  
# 
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)      -0.65573    0.01180 -55.565  < 2e-16 ***
#   selthreeAAG       0.50908    0.03536  14.396  < 2e-16 ***
#   selthreeAGG       0.40706    0.03183  12.788  < 2e-16 ***
#   selthreeCGT       0.53515    0.05345  10.012  < 2e-16 ***
#   selthreeGGC       0.10306    0.03766   2.737  0.00620 ** 
#   selthreeTGG      -0.32025    0.03993  -8.021 1.05e-15 ***
#   selthree.rightAA -0.11123    0.02519  -4.415 1.01e-05 ***
#   selthree.rightAC  0.36895    0.02418  15.260  < 2e-16 ***
#   selthree.rightGC  0.19123    0.02599   7.356 1.89e-13 ***
#   selthree.leftCA  -0.36240    0.02652 -13.663  < 2e-16 ***
#   selthree.leftCG   1.35426    0.02710  49.969  < 2e-16 ***
#   selthree.leftGA  -0.18397    0.02817  -6.530 6.57e-11 ***
#   selfive.rightAAG -0.34014    0.05115  -6.650 2.93e-11 ***
#   selfive.rightAGT  0.48865    0.04027  12.135  < 2e-16 ***
#   selfive.rightGCA  0.42117    0.03687  11.424  < 2e-16 ***
#   selfive.leftAGA  -0.37984    0.04983  -7.623 2.48e-14 ***
#   selfive.leftTCA   0.42857    0.04394   9.753  < 2e-16 ***
#   selfive.leftTGG   0.27345    0.03389   8.068 7.15e-16 ***
#   seldcgy           0.11848    0.04025   2.943  0.00325 ** 
#   selaaaaaawt       1.14238    0.06891  16.579  < 2e-16 ***
#   selayacacac       0.55623    0.05992   9.283  < 2e-16 ***
#   ---
#   Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 12295.6  on 1173  degrees of freedom
# Residual deviance:  2994.7  on 1153  degrees of freedom
# AIC: 6184.3
# 
# Number of Fisher Scoring iterations: 4

# intron
# Call:
#   glm(formula = y_data[, c("Success", "Failure")] ~ sel, family = binomial(link = "logit"))
# 
# Deviance Residuals: 
#   Min        1Q    Median        3Q       Max  
# -12.8587   -1.1907   -0.2836    0.8176   21.3793  
# 
# Coefficients: (1 not defined because of singularities)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)      -0.542554   0.007223 -75.114  < 2e-16 ***
#   seloneA          -0.020667   0.012184  -1.696  0.08985 .  
# seloneG                 NA         NA      NA       NA    
# selthreeAAA      -0.298946   0.016838 -17.754  < 2e-16 ***
#   selthreeAAG       0.489973   0.016463  29.763  < 2e-16 ***
#   selthreeAGG       0.293803   0.015395  19.084  < 2e-16 ***
#   selthreeAGT      -0.265269   0.019072 -13.909  < 2e-16 ***
#   selthreeCGT       0.318885   0.030201  10.559  < 2e-16 ***
#   selthreeGGC       0.231140   0.020144  11.475  < 2e-16 ***
#   selthreeTAC       0.116814   0.021056   5.548 2.89e-08 ***
#   selthreeTAG      -0.391405   0.022719 -17.228  < 2e-16 ***
#   selthreeTGG      -0.351192   0.018184 -19.313  < 2e-16 ***
#   selthree.rightAA -0.140696   0.012826 -10.969  < 2e-16 ***
#   selthree.rightAC  0.247276   0.013578  18.211  < 2e-16 ***
#   selthree.rightGC  0.027269   0.017198   1.586  0.11283    
# selthree.leftCA  -0.468734   0.013096 -35.791  < 2e-16 ***
#   selthree.leftCG   1.977038   0.018852 104.871  < 2e-16 ***
#   selthree.leftGA  -0.168310   0.013798 -12.198  < 2e-16 ***
#   selfive.rightAAG -0.311660   0.020686 -15.066  < 2e-16 ***
#   selfive.rightACC  0.258660   0.019710  13.123  < 2e-16 ***
#   selfive.rightAGT  0.482488   0.016408  29.407  < 2e-16 ***
#   selfive.rightATA -0.223318   0.018611 -11.999  < 2e-16 ***
#   selfive.rightGCA  0.440229   0.020282  21.705  < 2e-16 ***
#   selfive.rightGCC  0.301975   0.020826  14.500  < 2e-16 ***
#   selfive.leftAGA  -0.346831   0.020261 -17.118  < 2e-16 ***
#   selfive.leftAGG  -0.192611   0.017135 -11.241  < 2e-16 ***
#   selfive.leftGTA   0.265934   0.020160  13.191  < 2e-16 ***
#   selfive.leftTCA   0.465319   0.018136  25.657  < 2e-16 ***
#   selfive.leftTGG   0.157380   0.015898   9.899  < 2e-16 ***
#   seldncgy          0.150321   0.020725   7.253 4.07e-13 ***
#   selahatatay       0.327456   0.018847  17.374  < 2e-16 ***
#   selkkaaaaaa       0.783532   0.023276  33.663  < 2e-16 ***
#   selacacacrc       0.576368   0.031428  18.339  < 2e-16 ***
#   selaadttttk       0.435432   0.026782  16.258  < 2e-16 ***
#   selgtgtgy         0.040163   0.014596   2.752  0.00593 ** 
#   selaaaaawt        0.309662   0.021601  14.335  < 2e-16 ***
#   ---
#   Signif. codes:  0 ?***? 0.001 ?**? 0.01 ?*? 0.05 ?.? 0.1 ? ? 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 78482  on 7262  degrees of freedom
# Residual deviance: 18238  on 7228  degrees of freedom
# AIC: 31565
# 
# Number of Fisher Scoring iterations: 4

###
# barplot of coefficients against motifs/sequence contexts
###
n=ncol(sel)+1
df=data.frame(coef=coef(logr),feat=names(coef(logr)),type=numeric(n))
df$type=sapply(df$feat,FUN=function(x) {
  if (grepl("one",x)){
    "1mer"
  } else if (grepl("three",x) & grepl("right",x)){
    "3mer.rightflank"
  } else if (grepl("three",x) & grepl("left",x)){
    "3mer.leftflank"
  } else if (grepl("three",x)){
    "3mer"
  } else if (grepl("five",x) & grepl("right",x)){
    "5mer.rightflank"
  } else if (grepl("five",x) & grepl("left",x)){
    "5mer.leftflank"
  } else if (grepl("sel",x)){
    "motif"
  } else {"intercept"}
})
df=df[order(df$coef,decreasing=TRUE),]
df$feat=as.character(df$feat)
df$feat=ifelse(grepl("sel",df$feat),substr(df$feat,4,nchar(df$feat)),df$feat)
df$type=factor(df$type,levels=c("intercept","1mer","3mer.rightflank","3mer.leftflank","3mer","5mer.rightflank","5mer.leftflank","motif"))
df$feature=apply(df,1,function(x) {
  if(x[3]=="1mer"){
    paste(x[3],substr(x[2],4,4),sep="_")
  } else if (x[3]=="3mer.rightflank"){
    paste(x[3],substr(x[2],12,13),sep="_")
  } else if (x[3]=="3mer.leftflank"){
    paste(x[3],substr(x[2],11,12),sep="_")
  } else if (x[3]=="3mer"){
    paste(x[3],substr(x[2],6,8),sep="_")
  } else if (x[3]=="5mer.rightflank"){
    paste(x[3],substr(x[2],11,13),sep="_")
  } else if (x[3]=="5mer.leftflank"){
    paste(x[3],substr(x[2],10,12),sep="_")
  } else if (x[3]=="motif"){
    paste(x[3],x[2],sep="_")
  } else (x[3])
})
df=df[-which(df$type=="intercept"),] # 8

df$colour=as.character(df$type)
df$colour=ifelse(df$type %in% c("3mer.rightflank","3mer.leftflank"),"3mer.flank",df$colour)
df$colour=ifelse(df$colour %in% c("5mer.rightflank","5mer.leftflank"),"5mer.flank",df$colour)
df$colour=factor(df$colour,levels=c("1mer","3mer.flank","3mer","5mer.flank","motif"))
dff=df[order(df$colour,df$coef),]
dff$feature=factor(dff$feature,levels=dff$feature)
ggplot(dff,aes(x=feature,y=coef,fill=colour))+geom_bar(stat="identity")+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab("motif/sequence context")+ylab("regression coefficient")

save.image(file="dreme_final_enh.RData")
# save.image(file="dreme_final_3utr.RData")
# save.image(file="dreme_final_5utr.RData")
# save.image(file="dreme_final_prom.RData")
# save.image(file="dreme_final_intron.RData")
