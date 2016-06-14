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
