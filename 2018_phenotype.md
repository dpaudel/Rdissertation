Phenotypic data analysis 2018

```
setwd("~/Dropbox/2_UFlabAltpeter/2018_planting")
data2018 <- readxl::read_xlsx("phenotypicData/2018-08-08-03-15-51_fieldbook_map_database.xlsx")
```
#### assign cross type
```
data2018$cross_type <- substr(data2018$seed_name,1,1)
```
#### Remove missing values
```
data2018 <- data2018[data2018$cross_type !="*",]
```

```
data2018$cross_type <- as.factor(data2018$cross_type)
height <- data2018[data2018$trait=="height",]
tiller <- data2018[data2018$trait == "TillerNum",]
p_tiller <- ggplot(data=tiller, aes(x=seed_name, y=value))
p_tiller + geom_point(aes(colour = factor(cross_type)))+
  ylab("Number of tillers")+
  xlab("Line name")

p_height <- ggplot(data=height, aes(x=seed_name, y=value))
p_height + geom_point(aes(colour = factor(cross_type)))+
  ylab("Height (inch.)")+
  xlab("Line name")
```
