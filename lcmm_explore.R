setwd("~/R_directory/MAPS/6_14_16_data/")
library(gridExtra)
library(grid)
library(ggplot2)
library(lcmm)
library(plyr)
library(dplyr)
rm(list=ls())

###----------------- read in data -----------------------
### combinning the data
map <- read.csv(file="map_6_14_16.csv",sep=",",header=T)
nseq <- read.csv(file="otu_table_stats.csv",sep=",", header=T)
map2 <- inner_join(map,nseq)
stool <- read.csv(file="Updated_5_23_16 MAPS_stool_MCB_processing set up.csv", sep=",", header=T)
stool <- stool[,(1:5)]
map <- inner_join(map2, stool)
setdiff(map2$Fecal_Sample_ID_OLD, map$Fecal_Sample_ID_OLD)
map2 <- map %>% group_by(Daily_Code) %>% top_n(1, nsequence) %>% filter(nsequence>=10000) #pick the highest number of sequencing for each daily code
map <- select(map2, Sample, Subject_ID, Fecal_Sample_ID, PNA)
nbaby <- count(map, Subject_ID) # number of samples per baby

## read in the a-diversity
a.div <- read.csv(file="alpha_even_10000.csv",sep=",", header=T)
a.div <- inner_join(map, a.div)
a.div$PNA <- as.character(a.div$PNA)
a.div$Subject_ID <- as.character(a.div$Subject_ID)
a.div.30 <- subset(a.div, as.numeric(PNA)<31)

## demo, nnns
demo <- read.csv(file="MAPS demographic data_3_11_16.csv")
demo <- demo[,c("Subject_ID", "BIrth_GA", "Gender", "Baby_Race", "Delivery", "PROM", "Twins", "Birth_weight", "SNAPEII")]
demo$Subject_ID <- as.character(demo$Subject_ID)
nnns <- read.csv(file="NNNS_final_8_25_15.csv")
nnns$Subject_ID <- as.character(nnns$Subject_ID)
demonnns <- inner_join(demo, nnns)
a.div.30 <- left_join(a.div.30, demonnns)

## feeding
feeding <- read.csv(file="feeding_first_50.csv", sep=",", header=T)
feeding$Subject_ID <- as.character(feeding$Subject_ID)
feeding$PNA <- as.character(feeding$PNA)
count <- feeding %>% filter(as.numeric(PNA)<31) %>% group_by(Subject_ID,PNA)%>%summarise(summbm=sum(MBM, na.rm=TRUE), sumdbm=sum(DBM, na.rm=TRUE), sumformula=sum(Formula, na.rm=TRUE))
count <- count %>% mutate(pmbm=summbm/(summbm+sumdbm+sumformula),pdbm=sumdbm/(summbm+sumdbm+sumformula),pform=sumformula/(summbm+sumdbm+sumformula))
a.div.30 <- inner_join(a.div.30, count[,c("Subject_ID", "PNA", "pmbm","pdbm","pform")])

count2 <- feeding %>% filter(as.numeric(PNA)<11) %>%group_by(Subject_ID)%>%summarise(summbm=sum(MBM, na.rm=TRUE), sumdbm=sum(DBM, na.rm=TRUE), sumformula=sum(Formula, na.rm=TRUE))
count2 <- count2 %>% mutate(pmbm.10=summbm/(summbm+sumdbm+sumformula),pdbm.10=sumdbm/(summbm+sumdbm+sumformula),pform.10=sumformula/(summbm+sumdbm+sumformula)) 

count2$FT.10[count2$pmbm.10>=count2$pdbm.10 & count2$pmbm.10>=count2$pform.10] <- "mbm"
count2$FT.10[count2$pdbm.10>count2$pmbm.10 & count2$pdbm.10>count2$pform.10] <- "dbm"
count2$FT.10[count2$pform.10>count2$pmbm.10 & count2$pform.10 >count2$pdbm.10] <- "formula"
count2 <- select(count2, Subject_ID, pmbm.10, pdbm.10, pform.10,FT.10)
a.div.30 <- left_join(a.div.30, count2)

## antibiotics use, the first 10 days
anti <- read.csv(file="Daily_antibiotic.40.50.csv")
anti$Subject_ID <- as.character(anti$Subject_ID)
anti$PNA <- as.character(anti$PNA)
anti.10 <- anti%>% filter(as.numeric(PNA)<11) %>%count(Subject_ID)%>%rename(days.of.antibiotics.f10=n)
a.div.30 <- left_join(a.div.30, anti.10)
a.div.30$days.of.antibiotics.f10[is.na(a.div.30$days.of.antibiotics.f10)] <- 0

## niss
niss <- read.csv(file="MAPS_NISS_contact_daily_45.csv")
niss$Subject_ID <- as.character(niss$Subject_ID)
niss$PNA <- as.character(niss$PNA)
niss.daily  <- niss%>% mutate(freq.acute.daily=Acute5+Acute4+Acute3+Acute2, freq.chronic.daily=Chronic4+Chronic3+Chronic2)%>%select(Subject_ID,PNA, weighted.acute.daily, weighted.chronic.daily, freq.acute.daily, freq.chronic.daily, Kangaroo_Care_Daily, Breastfeeding_Daily, Contact_Total_Daily)
niss.10 <- niss %>% filter(PNA<11) %>% group_by (Subject_ID) %>% summarise(sum.weighted.acute.10=sum(weighted.acute.daily, na.rm=TRUE), sum.weighted.chronic.10=sum(weighted.chronic.daily, na.rm=TRUE), sum.kangroo.10=sum(Kangaroo_Care_Daily, na.rm=TRUE), sum.breast.10=sum(Breastfeeding_Daily, na.rm=TRUE))
a.div.30 <- left_join(a.div.30,niss.daily)

a.div.30$PNA <- as.numeric(a.div.30$PNA)
####---------------------LCMM  modeling -----------------------------
### PNA Feeding,first 10 days percentage

m1mom <- hlme(simpson~PNA+pmbm,random=~PNA,subject='Subject_ID',ng=1,data=a.div.30, cor="AR"(PNA), maxiter=5e3)
summary(m1mom)

m2mom<-hlme(simpson~PNA+pmbm,random=~PNA,mixture=~PNA,classmb=~pmbm.10, cor="AR"(PNA),subject='Subject_ID',ng=2,data=a.div.30,B=m1mom)
summary(m2mom)

people1 <- as.data.frame(m2mom$pprob[,1:2])
fd1 <- left_join(a.div.30, people1)
p1 <- ggplot(fd1, aes(PNA, simpson, group=Subject_ID, colour=as.character(fd1$class))) + geom_line() + geom_smooth(aes(group=as.character(fd1$class)), method="loess", size=2, se=F)  + scale_y_continuous(limits = c(0,1)) + labs(x="PNA",y="simpson",colour="Latent Class");p1
p2 <- ggplot(fd1, aes(PNA, simpson, group=Subject_ID, colour=as.character(fd1$class))) + geom_smooth(aes(group=Subject_ID, colour=as.character(fd1$class)),size=0.5, se=F) + geom_smooth(aes(group=as.character(fd1$class)), method="loess", size=2, se=T)  + scale_y_continuous(limits = c(0,1)) + labs(x="PNA",y="simpson",colour="Latent Class");p2

