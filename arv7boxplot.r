######################2019-04-15###################################################################
##Group 1: AR-V7 negative                                                                       ###
##Group 2: AR-V7 positive with expression in the bottom half of positive samples (ie AR-V7 low) ###
##Group 3: AR-V7 positive with expression in the top half of positive samples (ie AR-V7 high    ###
###################################################################################################  

library(tidyverse)
library(data.table)
library(ggpubr)
##AR_V7 negtive, positive (high), and positive (low)##
arv7<-read.delim("AR_V7_V9.stats",stringsAsFactors = FALSE)
arv7.neg<-arv7[arv7$AR_V7==0,]$sample_id #37
arv7.pos<-arv7[arv7$AR_V7>0,]
arv7.pos.hi<-arv7.pos[arv7.pos$AR_V7>=median(arv7.pos$AR_V7),]$sample_id #27
arv7.pos.low<-arv7.pos[arv7.pos$AR_V7<median(arv7.pos$AR_V7),]$sample_id #26

##TACSTD2 expressoin in above three groups##
dat<-read.delim("V1.TACSTD2.rpkm.count.txt",stringsAsFactors = FALSE)
datT<-transpose(dat[,]) #function from data.table
colnames(datT) <- rownames(dat)
rownames(datT) <- colnames(dat)
##or##
test<-setnames(datT, rownames(dat))
##remove 1row##
datT$sample_id<-row.names(datT)
datT$rpkm<-datT$`1`
datT$`1`<-NULL
datT<-datT[-1,]

##assign groups##
for (i in 1:nrow(datT)){
  datT$group[i]<-ifelse(datT$sample_id[i]%in%arv7.neg,"AR_V7 negtive",ifelse(datT$sample_id[i]%in%arv7.pos.hi,"AR_V7 positive (high)","AR_V7 positive (low)"))
}
datT$group<-as.factor(datT$group)
datT$rpkm<-as.numeric(datT$rpkm)
##plot##
p <- ggboxplot(datT, x = "group", y = "rpkm",
               color = "group", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
               add = "jitter", shape = "group",
               xlab=FALSE,
               ylab="TACSTD2 expressoin (RPKM)")
p
# Add p-values comparing groups
# Specify the comparisons you want
my_comparisons <- list( c("AR_V7 negtive", "AR_V7 positive (high)"), c("AR_V7 positive (high)", "AR_V7 positive (low)"), c("AR_V7 negtive", "AR_V7 positive (low)") )
p1<-p + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 300)    #global pvalue position





