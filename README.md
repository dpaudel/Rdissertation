# Rdissertation
Codes in R used for Dissertation

<h5>Plot % of uniquely mapped NG reads to PM genome</h5>
```
data1<-read.csv("forR.csv")
library(ggplot2)
colnames(data1)
qplot(Percentage.uniquely.mapped.to.Pmgenome*100, data = data1, geom = "histogram", xlim=c(0,50), main="Histogram of mapped reads", xlab="Percentage of good reads uniquely mapped to PM genome",ylab="Number of samples")
```
