###############2019-04-12#######################################
##TACSTD2 (RPKM) V1 vs V2 ######################################
##PROMOTE on cal42n#############################################
##/home/tywang/Projects/lncRNAs/PROMOTE/known_counts/results####
################################################################
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

samples<-read.delim("matched.samples",header=FALSE)
v1<-read.delim("V1.TACSTD2.rpkm.count.txt",header=TRUE)
v2<-read.delim("V2.TACSTD2.rpkm.count.txt",header=TRUE)

##reshape v1##
v1_t<-as.data.frame(t(v1)[-1,])
v1_t$rpkm<-as.numeric(as.character(v1_t$`t(v1)[-1, ]`))
v1_t$`t(v1)[-1, ]`<-NULL
v1_t$group<-"Baseline"
v1_t$sampleid<-row.names(v1_t)

##reshapeV2###
v2_t<-as.data.frame(t(v2)[-1,])
v2_t$rpkm<-as.numeric(as.character(v2_t$`t(v2)[-1, ]`))
v2_t$`t(v2)[-1, ]`<-NULL
v2_t$group<-"Week 12"
v2_t$sampleid<-row.names(v2_t)


##wide##
v1_t.m<-full_join(v1_t,samples,by=c("sampleid"="V1"))
v1v2_t.m<-full_join(v1_t.m,v2_t,by=c("V2"="sampleid"))
v1v2_t.m<-na.omit(v1v2_t.m)
v1v2_t.m$Baseline<-v1v2_t.m$rpkm.x
v1v2_t.m$rpkm.x<-NULL
v1v2_t.m$`Week 12`<-v1v2_t.m$rpkm.y
v1v2_t.m$rpkm.y<-NULL
v1v2_t.m$group.x<-v1v2_t.m$group.y<-NULL
v1v2_t.m$V2<-NULL

##long##
v1v2_t.m_long<-gather(v1v2_t.m,condition,RPKM,Baseline:`Week 12`,factor_key=TRUE)

for (i in 1:nrow(v1v2_t.m_long)){
  v1v2_t.m_long$treatment[i]<-strsplit(v1v2_t.m_long$sampleid[i],"_")[[1]][3]
}
v1v2_t.m_long$treatment<-as.factor(v1v2_t.m_long$treatment)

##plot##
mytheme <- theme_classic() %+replace% 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(face="bold",angle=90))  
  
p1 <- ggplot(data = v1v2_t.m_long, aes(x = condition, y = RPKM, group = sampleid, colour = sampleid)) +
  mytheme +
  #coord_trans(y="log10", limy=c(1000,6000)) +
  labs(y = "TACSTD2 (RPKM)") + 
  geom_line(size=1) + theme(legend.position="none")

##without color
p2 <- ggplot(data = v1v2_t.m_long, aes(x = condition, y = RPKM, group = sampleid)) +
  mytheme +
  #coord_trans(y="log10", limy=c(1000,6000)) +
  labs(y = "TACSTD2 (RPKM)") + 
  geom_line(size=0.5,alpha=0.5) + theme(legend.position="none")
  

########2019-05-10################
##boxplot with paired line########
##facet_wrap by responder or not##

input<-"tier2.info.txt"
respond<-read.delim(input,header = TRUE,check.names = FALSE,stringsAsFactors = FALSE)
respond<-respond[,c(7,16)]
v1v2_t.m_long.join<-full_join(v1v2_t.m_long,respond,by=c("sampleid"="sample_id"))
v1v2_t.m_long.join<-na.omit(v1v2_t.m_long.join)
nrow(v1v2_t.m_long.join[v1v2_t.m_long.join$condition=="Week 12",]) #63
nrow(v1v2_t.m_long.join[v1v2_t.m_long.join$composite_progression==0,]) #74
nrow(v1v2_t.m_long.join[v1v2_t.m_long.join$composite_progression==1,]) #52

for (i in 1:nrow(v1v2_t.m_long.join)){
  v1v2_t.m_long.join$type[i]<-ifelse(v1v2_t.m_long.join$composite_progression[i]==0,"Responder (N=37)","Non-responder (N=26)")
}

p<-ggpaired(v1v2_t.m_long.join, x = "condition", y = "RPKM",
         color = "condition", line.color ="gray", line.size = 0.4,
         palette = "npg")+
         stat_compare_means(paired = TRUE,label.y = 210)+
         facet_wrap(~type)+
         labs(y = "TACSTD2 (RPKM)")

# Use only p.format as label. Remove method name.
p + stat_compare_means(label = "p.format", paired = TRUE)


##without boxplot##
p2 <- ggplot(data = v1v2_t.m_long.join, aes(x = condition, y = RPKM, group = sampleid)) +
  mytheme +
  #coord_trans(y="log10", limy=c(1000,6000)) +
  labs(y = "TACSTD2 (RPKM)") + 
  geom_line(size=0.5,alpha=0.5) + theme(legend.position="none")+
  facet_wrap(~type)



