# Rdissertation
Codes in R used for Dissertation

<h5>Plot % of uniquely mapped NG reads to PM genome</h5>
```
data1<-read.csv("forR.csv")
library(ggplot2)
colnames(data1)
qplot(Percentage.uniquely.mapped.to.Pmgenome*100, data = data1, geom = "histogram", xlim=c(0,50), main="Histogram of mapped reads", xlab="Percentage of good reads uniquely mapped to PM genome",ylab="Number of samples")
```
Plot with error bars
```
g<-ggplot(data1, aes(trt,bio_yield))+
  stat_summary(fun.y=mean, geom="bar",fill="grey50")+
  stat_summary(fun.data=mean_cl_normal, geom="errorbar")
```
Plot with boot strap bars
```
ggplot(data1, aes(trt,bio_yield))+
  stat_summary(fun.y=mean, geom="bar",fill="grey40")+
  stat_summary(fun.data=mean_cl_boot, geom="pointrange")
```
Combining all results of gbs snps fro summarizing
```
setwd("C:/Users/Dev Paudel/Dropbox/1_WangLab/GBS/snp_summary")
library(plyr)
library(VariantAnnotation)
library(compare)
library(sqldf)

#tassel
tassel<-read.table("tassel43_allcustomSNPlog.txt", header=T)
l_tassel<-tassel[,c("Chr","SNPPosition","Alleles","nTagsAtLocus","nReads","nTaxaCovered")]
l_tassel$Locus=paste(l_tassel$Chr,l_tassel$SNPPosition, sep="_")
colnames(l_tassel)<-c("tas_Chr","tas_SNPPosition","tas_Alleles","tas_nTagsAtLocus","tas_nReads","tas_nTaxaCovered","Locus")

#GBS-SNP-CROP
snpcrop<-read.table("snpcrop.Tassel_in.hmp.txt", header=T)
l_snpcrop <- snpcrop[,c("rs","chrom","pos","N190","N122")]
colnames(l_snpcrop)<-c("Locus","crop_chr","crop_pos","crop_N190","crop_N122")

#GATK
l_gatk<-read.table("locus_gatkq30dp11.txt", header=T)
l_gatk$Locus=paste(l_gatk$CHROM,l_gatk$POS, sep="_")
colnames(l_gatk)<-c("gatk_chr","gatk_pos","gatk_ref","gatk_alt","Locus")

#SAMTOOLS
l_sam <-read.table("locus_samtoolsq30dp11.txt", header=T)
l_sam$Locus=paste(l_sam$CHROM,l_sam$POS, sep="_")
colnames(l_sam)<-c("sam_chr","sam_pos","sam_ref","sam_alt","Locus")

#Freebayes
l_free <-read.table("locus_freebayes5covq30dp11.txt", header=T)
l_free$Locus=paste(l_free$CHROM,l_free$POS, sep="_")
colnames(l_free)<-c("free_chr","free_pos","free_ref","free_alt","Locus")

#Extract common snps
View(merge(l_free, l_sam, by="Locus"))
```
