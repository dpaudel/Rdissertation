# 2018 phenotypic data analysis

setwd("~/Dropbox/2_UFlabAltpeter/2018_planting")
library(tidyverse)
library(agricolae)
library(dpaudel)
library(lawstat)
library(extrafont)
library(car)
loadfonts(device = "win")
windowsFonts(Arial=windowsFont("TT Arial"))

data2018 <- read.csv("phenotypicData/2018-10-06_map_database_complete.csv")
# assign cross type
#data2018$cross_type <- substr(data2018$seed_name,1,1)

# Remove missing values
data2018 <- data2018[data2018$cross_type !="*",]
data2018$cross_type <- as.factor(data2018$cross_type)
data2018$block <- as.factor(data2018$block)
data2018$seed_name <- as.factor(data2018$seed_name)
# Summarize data; CAUTION: leaf lenght and height is in inches and biomass is in lbs
summary_data18 <- data2018 %>% group_by(trait, cross_type) %>% summarize(mean=mean(value), min=min(value), max=max(value))
copy2clipboard(summary_data18)
# Subset data
height <- data2018[data2018$trait=="height",]
tiller <- data2018[data2018$trait == "TillerNum",]
stemdiam <- data2018[data2018$trait=="StemDiam",]
LeafLength <- data2018[data2018$trait=="LeafLength",]
LeafWidth <-  data2018[data2018$trait=="LeafWidth",]
biomass <- data2018[data2018$trait=="BioMass",]
sample_wt <- data2018[data2018$trait=="Sample_wt_gm",]
sample_drywt <- data2018[data2018$trait=="Dry_wt_gm",]
conversionRate <- data2018[data2018$trait=="ConversionRate",]

#### Plant height ####

height$value_cm <- height$value*2.54 # convert inch to cm

# Levene test - median
leveneTest(value_cm~cross_type, data=height) # 2.702e-15 ***
leveneTest(value_cm~seed_name, data=height[height$cross_type=="A",]) # 0.07781
leveneTest(value_cm~seed_name, data=height[height$cross_type=="B",]) # 0.9584
leveneTest(value_cm~seed_name, data=height[height$cross_type=="C",]) # 0.8109
leveneTest(value_cm~seed_name, data=height[height$cross_type=="D",]) # 0.6924

# scatterplot
p_height <- ggplot(data=height, aes(x=seed_name, y=value_cm))
p_height + geom_point(aes(colour = factor(cross_type)))+
  ylab("Plant height (cm)")+ xlab("PMN hybrid lines")+ theme_uf() + theme(axis.ticks = element_blank(), axis.text.x = element_blank())

# Boxplot
ggplot(height, aes(x=cross_type, y=value_cm)) +
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(position=position_jitter(width=.1, height=0), aes(colour=cross_type))+ theme_uf()+
  ylab("Plant height (cm)")+xlab("Cross Type")+
  theme(legend.position="none") # remove legend

# Facet wrap
ggplot(height, aes(value_cm))+ geom_histogram(binwidth = 10, color="black", fill="grey40") +
  theme_uf()+ylab("")+xlab("Plant height (cm)")+
  facet_wrap(~cross_type)


#### Tiller Number ####

# Levene test - median
leveneTest(value~cross_type, data=tiller) # < 2.2e-16 ***
leveneTest(value~seed_name, data=tiller[tiller$cross_type=="A",]) # 0.07578
leveneTest(value~seed_name, data=tiller[tiller$cross_type=="B",]) # 0.6613
leveneTest(value~seed_name, data=tiller[tiller$cross_type=="C",]) # 0.1015
leveneTest(value~seed_name, data=tiller[tiller$cross_type=="D",]) # 0.6265

# scatterplot
p_tiller <- ggplot(data=tiller, aes(x=seed_name, y=value))
p_tiller + geom_point(aes(colour = factor(cross_type)))+
  ylab("Number of tillers")+ xlab("PMN hybrid lines")+ theme_uf() + theme(axis.ticks = element_blank(), axis.text.x = element_blank())

# Boxplot
ggplot(tiller, aes(x=cross_type, y=value)) +
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(position=position_jitter(width=.1, height=0), aes(colour=cross_type))+ theme_uf()+
  ylab("Number of tillers")+xlab("Cross Type")+
  theme(legend.position="none") # remove legend

# Facet wrap
ggplot(tiller, aes(value))+ geom_histogram(binwidth = 2, color="black", fill="grey40") +
  theme_uf()+ylab("")+xlab("Number of tillers")+
  facet_wrap(~cross_type)




#### Stem Diameter ####

# Levene test - median
leveneTest(value~cross_type, data=stemdiam) # 1.194e-05 ***
leveneTest(value~seed_name, data=stemdiam[stemdiam$cross_type=="A",]) # 0.3174
leveneTest(value~seed_name, data=stemdiam[stemdiam$cross_type=="B",]) # 0.9829
leveneTest(value~seed_name, data=stemdiam[stemdiam$cross_type=="C",]) # 0.4913
leveneTest(value~seed_name, data=stemdiam[stemdiam$cross_type=="D",]) # 0.8523

# scatterplot
p_stemdiam <- ggplot(data=stemdiam, aes(x=seed_name, y=value))
p_stemdiam + geom_point(aes(colour = factor(cross_type)))+
  ylab("Stem diameter (mm)")+ xlab("PMN hybrid lines")+ theme_uf() + theme(axis.ticks = element_blank(), axis.text.x = element_blank())

# Boxplot
ggplot(stemdiam, aes(x=cross_type, y=value)) +
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(position=position_jitter(width=.1, height=0), aes(colour=cross_type))+ theme_uf()+
  ylab("Stem diameter (mm)")+xlab("Cross Type")+
  theme(legend.position="none") # remove legend

# Facet wrap
ggplot(stemdiam, aes(value))+ geom_histogram(binwidth = 1, color="black", fill="grey40") +
  theme_uf()+ylab("")+xlab("Stem diameter (mm)")+
  facet_wrap(~cross_type)


#### Leaf length ####
LeafLength$value_cm <- LeafLength$value*2.54 # convert inch to cm

# Levene test - median
leveneTest(value_cm~cross_type, data=LeafLength) # < 2.2e-16 ***
leveneTest(value_cm~seed_name, data=LeafLength[LeafLength$cross_type=="A",]) # 0.2382
leveneTest(value_cm~seed_name, data=LeafLength[LeafLength$cross_type=="B",]) # 0.5612
leveneTest(value_cm~seed_name, data=LeafLength[LeafLength$cross_type=="C",]) # 0.6267
leveneTest(value_cm~seed_name, data=LeafLength[LeafLength$cross_type=="D",]) # 0.1045

# scatterplot
p_LeafLength <- ggplot(data=LeafLength, aes(x=seed_name, y=value_cm))
p_LeafLength + geom_point(aes(colour = factor(cross_type)))+
  ylab("Leaf length (cm)")+ xlab("PMN hybrid lines")+ theme_uf() + theme(axis.ticks = element_blank(), axis.text.x = element_blank())

# Boxplot
ggplot(LeafLength, aes(x=cross_type, y=value_cm)) +
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(position=position_jitter(width=.1, height=0), aes(colour=cross_type))+ theme_uf()+
  ylab("Leaf length (cm)")+xlab("Cross Type")+
  theme(legend.position = "none")

# Facet wrap
ggplot(LeafLength, aes(value_cm))+ geom_histogram(binwidth = 3, color="black", fill="grey40") +
  theme_uf()+ylab("")+xlab("Leaf length (cm)")+
  facet_wrap(~cross_type)

#### Leaf Width ####

# Levene test - median
leveneTest(value~cross_type, data=LeafWidth) # 0.05293 .
leveneTest(value~seed_name, data=LeafWidth[LeafWidth$cross_type=="A",]) # 0.4596
leveneTest(value~seed_name, data=LeafWidth[LeafWidth$cross_type=="B",]) # 0.8497
leveneTest(value~seed_name, data=LeafWidth[LeafWidth$cross_type=="C",]) # 0.7394
leveneTest(value~seed_name, data=LeafWidth[LeafWidth$cross_type=="D",]) # 0.1875

# scatterplot
p_LeafWidth <- ggplot(data=LeafWidth, aes(x=seed_name, y=value))
p_LeafWidth + geom_point(aes(colour = factor(cross_type)))+
  ylab("Leaf width (mm)")+ xlab("PMN hybrid lines")+ theme_uf() + theme(axis.ticks = element_blank(), axis.text.x = element_blank())

# Boxplot
ggplot(LeafWidth, aes(x=cross_type, y=value)) +
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(position=position_jitter(width=.1, height=0), aes(colour=cross_type))+ theme_uf()+
  ylab("Leaf width (mm)")+xlab("Cross Type")+
  theme(legend.position="none") # remove legend

# Facet wrap
ggplot(LeafWidth, aes(value))+ geom_histogram(binwidth = 2, color="black", fill="grey40") +
  theme_uf()+ylab("")+xlab("Leaf width (mm)")+
  facet_wrap(~cross_type)

#### Biomass - sample not used <use total_biomass below> ####
biomass$value_kg <- biomass$value*0.453592 # Convert lbs to kgs
# Levene test - median
leveneTest(value_kg~cross_type, data=biomass) # < 2.2e-16 ***
leveneTest(value_kg~seed_name, data=biomass[biomass$cross_type=="A",]) # 0.079 .
leveneTest(value_kg~seed_name, data=biomass[biomass$cross_type=="B",]) # 0.7688
leveneTest(value_kg~seed_name, data=biomass[biomass$cross_type=="C",]) # 0.7688
leveneTest(value_kg~seed_name, data=biomass[biomass$cross_type=="D",]) # 0.02005 *

# scatterplot
p_biomass <- ggplot(data=biomass, aes(x=seed_name, y=value_kg))
p_biomass + geom_point(aes(colour = factor(cross_type)))+
  ylab("Biomass")+ xlab("PMN hybrid lines")+ theme_uf() + theme(axis.ticks = element_blank(), axis.text.x = element_blank())

# Boxplot
ggplot(biomass, aes(x=cross_type, y=value_kg)) +
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(position=position_jitter(width=.1, height=0), aes(colour=cross_type))+ theme_uf()+
  ylab("Fresh Biomass (kg/plant)")+xlab("Cross Type")+
  theme(legend.position="none") # remove legend

# Facet wrap
ggplot(biomass, aes(value_kg))+ geom_histogram(binwidth = 1, color="black", fill="grey40") +
  theme_uf()+ylab("")+xlab("Fresh Biomass (kg/plant)")+
  facet_wrap(~cross_type)

#### Biomass - added samples ####
biomass$value_kg <- biomass$value*0.453592 # Convert lbs to kgs
sample_wt$value_kg <- sample_wt$value/1000 # convert gm to kgs
sample_wt$value_kg[is.na(sample_wt$value_kg)] <- 0
total_biomass <- cbind(biomass[,],sample_wt$value_kg)
#View(total_biomass)
total_biomass$total_biomass_kg <- total_biomass$value_kg+total_biomass$`sample_wt$value_kg`
# the final biomass to use in total_biomass$total_biomass_kg
# Levene test - median
leveneTest(total_biomass_kg~cross_type, data=total_biomass) # < 2.2e-16 ***
leveneTest(total_biomass_kg~seed_name, data=total_biomass[total_biomass$cross_type=="A",]) # 0.02145 *
leveneTest(total_biomass_kg~seed_name, data=total_biomass[total_biomass$cross_type=="B",]) # 0.6081
leveneTest(total_biomass_kg~seed_name, data=total_biomass[total_biomass$cross_type=="C",]) # 0.07944 .
leveneTest(total_biomass_kg~seed_name, data=total_biomass[total_biomass$cross_type=="D",]) # 0.7515

# scatterplot
p_biomass2 <- ggplot(data=total_biomass, aes(x=seed_name, y=total_biomass_kg))
p_biomass2 + geom_point(aes(colour = factor(cross_type)))
p_biomass2 + geom_point(aes(colour = factor(cross_type)))+
  ylab("Biomass")+ xlab("PMN hybrid lines")+ theme_uf() + theme(axis.ticks = element_blank(), axis.text.x = element_blank())

# Boxplot
ggplot(total_biomass, aes(x=cross_type, y=total_biomass_kg)) +
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(position=position_jitter(width=.1, height=0), aes(colour=cross_type))+ theme_uf()+
  ylab("Fresh Biomass (kg/plant)")+xlab("Cross Type")+
  theme(legend.position="none") # remove legend

# Facet wrap
ggplot(total_biomass, aes(total_biomass_kg))+ geom_histogram(binwidth = 1, color="black", fill="grey40") +
  theme_uf()+ylab("")+xlab("Fresh Biomass (kg/plant)")+
  facet_wrap(~cross_type)


#### Biomass tons per ha ####
# 1 plant area = 12 sq. ft = 0.000111484 ha.

total_biomass$tons_per_ha <- (total_biomass$total_biomass_kg/0.000111484)/1000
hist(total_biomass$tons_per_ha)

#### drybiomass tons per ha ####
total_biomass_dm <- cbind(total_biomass[,],conversionRate$value)
#View(total_biomass_dm)
total_biomass_dm$dm_tpha <- (total_biomass_dm$total_biomass_kg*total_biomass_dm$`conversionRate$value`/0.000111484)/1000
summary(total_biomass_dm$dm_tpha)
hist(total_biomass_dm$dm_tpha)

# Levene test - median
leveneTest(dm_tpha~cross_type, data=total_biomass_dm) # < 2.2e-16 ***
leveneTest(dm_tpha~seed_name, data=total_biomass_dm[total_biomass_dm$cross_type=="A",]) # 0.02077 *
leveneTest(dm_tpha~seed_name, data=total_biomass_dm[total_biomass_dm$cross_type=="B",]) # 0.6462
leveneTest(dm_tpha~seed_name, data=total_biomass_dm[total_biomass_dm$cross_type=="C",]) # 0.1246
leveneTest(dm_tpha~seed_name, data=total_biomass_dm[total_biomass_dm$cross_type=="D",]) # 0.9255

# scatterplot
p_dm_tpha <- ggplot(data=total_biomass_dm, aes(x=seed_name, y=dm_tpha))
p_dm_tpha + geom_point(aes(colour = factor(cross_type)))
p_dm_tpha + geom_point(aes(colour = factor(cross_type)))+
  ylab("Dry Biomass (tons/ha.")+ xlab("PMN hybrid lines")+ theme_uf() + theme(axis.ticks = element_blank(), axis.text.x = element_blank())

# Boxplot
ggplot(total_biomass_dm, aes(x=cross_type, y=dm_tpha)) +
  geom_boxplot(outlier.shape=NA) + #avoid plotting outliers twice
  geom_jitter(position=position_jitter(width=.1, height=0), aes(colour=cross_type))+ theme_uf()+
  ylab("Dry Biomass (tons/ha.)")+xlab("Cross Type")+
  theme(legend.position="none") # remove legend

# Facet wrap
ggplot(total_biomass_dm, aes(dm_tpha))+ geom_histogram(binwidth = 1, color="black", fill="grey40") +
  theme_uf()+ylab("")+xlab("Dry Biomass (tons/ha)")+
  facet_wrap(~cross_type)



#############################################################
############ AOV ####
stemdiam1 <- aov(value~cross_type, data=stemdiam)
summary(stemdiam1)
stemdiam1
diagnostics(residuals(stemdiam1),fitted(stemdiam1))
hsd_stemdiam1 <- HSD.test(stemdiam1, "cross_type", group=T)
hsd_stemdiam1
stemdiamA <- stemdiam[stemdiam$cross_type=="A",]
(stemdiamA1 <- aov(value~seed_name, data=stemdiamA))
summary(stemdiamA1)
diagnostics(residuals(stemdiamA1), fitted(stemdiamA1))
hsd_stemdiamA1 <- HSD.test(stemdiamA1, "seed_name", group=T)
hsd_stemdiamA1
cv.model(stemdiam1)
cv.model(stemdiamA1)

#### CV ####
# Extract mean, sd, and cv; cv is same for cm vs inch or kg vs lbs
meansd_2018 <- data2018 %>% group_by(trait, cross_type) %>% summarise(mean=mean(value), sd=sd(value), cv=sd(value)*100/mean(value), min=min(value), max=max(value))

ggplot(meansd_2018, aes(y=cv, x=cross_type,  fill=cross_type))+  geom_bar(stat='identity')+
  facet_wrap(~trait)+ theme_uf()+
  xlab("Cross type") + ylab("Coefficient of variation (%)")+
  theme(legend.position="none")


# Summarise values , done individually because some labels have been transformed to cm and kg
height %>% group_by(cross_type) %>% summarise(mean=mean(value_cm), sd=sd(value_cm), cv=sd(value_cm)*100/mean(value_cm), min=min(value_cm), max=max(value_cm)) %>% View()
LeafLength %>%  group_by(cross_type) %>% summarise(mean=mean(value_cm), sd=sd(value_cm), cv=sd(value_cm)*100/mean(value_cm), min=min(value_cm), max=max(value_cm)) %>% View()
total_biomass %>% group_by(cross_type) %>% summarise(mean=mean(total_biomass_kg), sd=sd(value_kg), cv=sd(value_kg)*100/mean(value_kg), min=min(value_kg), max=max(value_kg)) %>% View()
# After adding sample biomass to first 3 replicates
total_biomass %>% group_by(cross_type) %>% summarise(mean=mean(total_biomass_kg), sd=sd(total_biomass_kg), cv=sd(total_biomass_kg)*100/mean(total_biomass_kg), min=min(total_biomass_kg), max=max(total_biomass_kg)) %>% View()
tiller %>% group_by(cross_type) %>% summarise(mean=mean(value), sd=sd(value), cv=sd(value)*100/mean(value), min=min(value), max=max(value)) %>% View()
stemdiam %>% group_by(cross_type) %>% summarise(mean=mean(value), sd=sd(value), cv=sd(value)*100/mean(value), min=min(value), max=max(value)) %>% View()
LeafWidth %>% group_by(cross_type) %>% summarise(mean=mean(value), sd=sd(value), cv=sd(value)*100/mean(value), min=min(value), max=max(value)) %>% View()
total_biomass_dm %>% group_by(cross_type) %>% summarise(mean=mean(dm_tpha), sd=sd(dm_tpha), cv=sd(dm_tpha)*100/mean(dm_tpha), min=min(dm_tpha), max=max(dm_tpha)) %>% View()

#### AOV ####
## height
height1 <- aov(value_cm~cross_type*block, data=height)
diagnostics(resid=height1$residuals, fitvalue = height1$fitted.values)
summary(height1)
(HSD.test(height1, "cross_type"))$groups
(HSD.test(height1, "block"))$groups

## tiller number
tiller1 <- aov(value~cross_type*block, data=tiller)
diagnostics(resid=tiller1$residuals, fitvalue = tiller1$fitted.values)
summary(tiller1)
(HSD.test(tiller1, "cross_type"))$groups


## Stem diameter
stemdiam1 <- aov(value~cross_type*block, data=stemdiam)
diagnostics(resid=stemdiam1$residuals, fitvalue = stemdiam1$fitted.values)
summary(stemdiam1)
(HSD.test(stemdiam1, "cross_type"))$groups


## Leaf length
leaflength1 <- aov(value_cm~cross_type*block, data=LeafLength)
diagnostics(resid=leaflength1$residuals, fitvalue = leaflength1$fitted.values)
summary(leaflength1)
(HSD.test(leaflength1, "cross_type"))$groups


## Leaf Width
leafwidth1 <- aov(value~cross_type*block, data=LeafWidth)
diagnostics(resid=leafwidth1$residuals, fitvalue = leafwidth1$fitted.values)
summary(leafwidth1)
(HSD.test(leafwidth1, "cross_type"))$groups


## Fresh biomass per plant
biomass1 <- aov(value_kg~cross_type*block, data=biomass)
diagnostics(resid=biomass1$residuals, fitvalue = biomass1$fitted.values)
summary(biomass1)
(HSD.test(biomass1, "cross_type"))$groups

## Dry biomass tons per ha
dm1 <- aov(dm_tpha~cross_type*block, data=total_biomass_dm)
diagnostics(resid=dm1$residuals, fitvalue = dm1$fitted.values)
summary(dm1) # p = 3.03e-10 ***
(HSD.test(dm1, "cross_type"))$groups
###################################################

sumheight <- summarySE(height, measurevar="value_cm", groupvars=c("cross_type","block"))
pd <- position_dodge(0.1)
ggplot(sumheight, aes(x=cross_type, y=value_cm, colour=block)) +
  geom_errorbar(aes(ymin=sumheight$value_cm - sumheight$ci, ymax= sumheight$value_cm + sumheight$ci), width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd)


ggplot(sumheight, aes(x=cross_type, y=value_cm, fill=block)) +
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=value_cm-ci, ymax=value_cm+ci),
                width=.2,                    # Width of the error bars
                position=position_dodge(.9))+
  theme_uf()

hsdsummary <- HSD.test(height1, "cross_type")
hsdsummary$cross_type <- hsdsummary$groups
##### To add
ggplot(height, aes(x=hsdsummary$cross_type, y=hsdsummary$value_cm, fill=hsdsummary$cross_type)) +
  geom_bar(position=position_dodge(), stat="identity",
           colour="black", # Use black outlines,
           size=.3) + scale_fill_brewer() +      # Thinner lines
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL),
                size=.5,    # Thinner lines
                width=.2,
                position=position_dodge(.9)) +
  ylab("Plant height 'cm'") +
  ylim(0,0.5)+
  theme_bw()+
  theme(axis.title=element_text( size = rel(1.8) , color="black")) +
  theme(axis.text = element_text(colour = "black", size="18")) +
  theme(legend.title = element_text(size = "18"))+
  theme(legend.text = element_text(size = "18")) +
  theme(plot.title = element_text(size = rel(3)))+
  theme(axis.title.x = element_blank())+
  geom_text(data = hsdsummary, mapping = aes(label = M, y = 0.3),angle=0, hjust = 0)

