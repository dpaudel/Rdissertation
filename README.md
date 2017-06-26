# Rdissertation

Codes in R used for Dissertation
<h5>Plots & summary for 2016 phenotyping</h5>

```
setwd("C:/Users/dpaudel/Dropbox/2_UFlabAltpeter/2016 phenotyping")
flower16<-read.csv("floweringdate_2016.csv", header=T)
library(agricolae)
library(plyr)
head(flower16)
str(flower16)
plot(flower16$DAP,flower16$Line)
myaov1<-lm(flower16$DAP~flower16$Line)
summary(myaov1)
#plot(myaov1)
cdata <- ddply(flower16, c("Line"), summarise,mean = mean(DAP),sd=sd(DAP),se= sd/sqrt(3))
cdata
plot(hist(cdata$mean), col="blue", main="Histogram for first flowering date 2016 mapping population", xlab="Days after planting")
abline(v=373.33, lwd=2,col="red" ); text(371, 30, "N122")
abline(v=420, lwd=2,col="green" ); text(418, 30, "N190")
#write.csv(cdata, file="summary2016.csv")

#For germplasm collection
germ16<-read.csv("flowering_germplasm2016.csv",header=T)
myaov2<-lm(germ16$DAP~germ16$Genotype)
summary(myaov2)
plot(myaov2)
#Remove NA for getting true means
germ16_na<-na.omit(germ16)
gdata <- ddply(germ16_na, c("Genotype"), summarise,mean = mean(DAP),sd=sd(DAP),se= sd/sqrt(3))
head(gdata)
#View(gdata)
plot(hist(gdata$mean), col="blue", main="Histgram of first flowering date 2016 - Germplasm", xlab="Days after planting")
abline(v=464, lwd=2,col="red" ); text(462,20, "N122", col="red")
abline(v=486, lwd=2,col="green" ); text(488, 20, "N190", col="green")
#write.csv(gdata, file="summary2016_germplasm.csv")
#using aggregate to group data
aggregate(germ16_na$DAP, by=list(germ16_na$Genotype), FUN=function(x)c(mean=mean(x),sd=sd(x)))

```

<h5>Merge all data for linkage map</h5>

```
setwd("C:/Users/Dev Paudel/Desktop/vcf/filter_48")
library(plyr)
dtassel<-read.table("dtassel_final.txt", header=T)
dcrop<-read.table("dcrop_final.txt",header=T)
crop<-read.table("crop_final.txt", header=T)
samtools<-read.table("samtools_final.txt", header=T)
gatk<-read.table("gatk_final.txt", header=T)
tassel<-read.table("tassel_final.txt", header=T)
combined<-rbind.fill(dtassel,dcrop,crop,samtools,gatk,tassel)
dim(combined)
head(combined)
write.csv(combined, file="combined_file.csv")
```

<h5>Plots of various generations with plant heigth (for 2nd committee meeting) </h5>

```
setwd("C:/Users/Dev Paudel/Dropbox/2_UFlabAltpeter/phenoytping_data")
raw<-read.csv("data_all.csv")
library(dplyr)
dim(raw)
colnames(raw)
str(raw)
#Extract only relevant data
data<-raw[,8:18]
theme_set(theme_gray(base_size = 18))

#Subset for 787 and plot
sub787<- subset(data,Trt=="787" | Trt =="787_S2" | Trt=="787_S3"| Trt=="787_S4")
p<-ggplot(data=sub787, aes(x=Trt, y=Height_cm))
p+geom_boxplot()+ geom_jitter(width=0.2)+
  theme_grey(base_size = 18)+
  xlab("Generation of 787")+
  ylab("Plant Height (cm)")


#Subset for 603 and plot
sub603<- subset(data, Trt=="603"| Trt=="603_S2"|Trt=="603_S3"|Trt=="603_S4")
q<-ggplot(data=sub603, aes(x=Trt, y=Height_cm))
q+geom_boxplot()+ geom_jitter(width=0.2)+
  theme_grey(base_size = 18)+
  xlab("Generation of 603")+
  ylab("Plant height (cm)")

#Subset for MS787 and plot
sub_ms787 <- subset(data, Trt=="MS787_S1"|Trt=="MS787_S2"|Trt=="MS787_S3")
r<-ggplot(data=sub_ms787, aes(x=Trt, y=Height_cm))
r+geom_boxplot()+geom_jitter(width=0.2)+
  theme_grey(base_size = 18)+
  xlab("Generation of back-crosses")+
  ylab("Plant height (cm)")

#Subset for Schank
sub_schank<-subset(data, Trt=="Schank"| Trt=="Schank_S1")
s<-ggplot(data=sub_schank,aes(x=Trt, y=Height_cm))
s+geom_boxplot()+geom_jitter(width=0.2)+
  theme_grey(base_size=18)+
  xlab("Generation of Schank")+
  ylab("Plant height (cm)")
```

<h5>Plot % of uniquely mapped NG reads to PM genome</h5>

```
setwd("C:/Users/Dev Paudel/Desktop/stacks")
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

<h5>Summarise data of 2012 and 2013 phenotypes</h5>

```

#For 2012 data
flower12<-read.csv("flowering_2012_14.csv", header=T)
head(flower12)
str(flower12)
myaov12<-lm(flower12$DAP2012~flower12$Line)
summary(myaov12)
colnames(flower12)
summary12<-aggregate(flower12$DAP2012, by=list(flower12$Line), FUN=function(x)c(mean=mean(x),sd=sd(x)))
#write.csv(summary12, file="summary2012.csv")

#For 2013 data

myaov13<-lm(flower12$DAP2013~flower12$Line)
summary(myaov13)
summary13<-aggregate(flower12$DAP2013, by=list(flower12$Line), FUN=function(x)c(mean=mean(x),sd=sd(x)))
write.csv(summary13, file="summary2013.csv")
```

