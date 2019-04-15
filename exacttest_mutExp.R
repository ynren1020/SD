################2019-04-12##################################
##AR,TP53,PTEN,RB1 mutation vs TACSTD2(high and low)########
##scott paper for revision##################################
##/home/yren/project/scott/paper/mutation###################
############################################################
library(dplyr)
library(tidyr)

v1<-read.delim("V1.TACSTD2.rpkm.count.txt",header=TRUE)
##reshape v1##
v1_t<-as.data.frame(t(v1)[-1,])
v1_t$rpkm<-as.numeric(as.character(v1_t$`t(v1)[-1, ]`))
v1_t$`t(v1)[-1, ]`<-NULL
v1_t$group<-"Baseline"
v1_t$sampleid<-row.names(v1_t)
##sample mutation info##
tier<-read.delim("tier2.info.txt")
tier$EXN<-toupper(tier$exn)
tiersub<-tier[,c("sample_id","EXN")]

##join two datasets##
v1_mu<-full_join(tiersub,v1_t,by=c("sample_id"="sampleid"))
v1_mu<-na.omit(v1_mu)

##subset v1_mu by TACSTD1 median rpkm  44.81262##
v1_mu_high<-filter(v1_mu,rpkm>= 44.81262)#38
v1_mu_low<-filter(v1_mu,rpkm< 44.81262)#38

exn_hi<-as.data.frame(v1_mu_high$EXN)
exn_low<-as.data.frame(v1_mu_low$EXN)

write.table(exn_hi,"exn_hi.txt",quote=FALSE,col.names = FALSE,row.names = FALSE,sep="\t")
write.table(exn_low,"exn_low.txt",quote=FALSE,col.names = FALSE,row.names = FALSE,sep="\t")

##test##
##AR##
genemat <- matrix(c(20, 21,17, 16), nrow = 2,
                  dimnames =
                    list(c("TACSTD1.low", "TACSTD1.high"), 
                         c("AR.mut", "AR.nomut")))
AR<-as.data.frame(genemat)
write.table(AR,"AR.fisher.txt",sep="\t",quote=FALSE,col.names = TRUE,row.names = TRUE)
##matrix look like this:
#           sv.mut sv.nomut
#var.low       20       17
#var.high      21       16

p.val<-fisher.test(genemat, alternative = "two.sided")$p.value #1

##TP53##
genemat <- matrix(c(12, 13,25, 24), nrow = 2,
                  dimnames =
                    list(c("TACSTD1.low", "TACSTD1.high"),
                         c("TP53.mut", "TP53.nomut")))
TP53<-as.data.frame(genemat)
write.table(TP53,"TP53.fisher.txt",sep="\t",quote=FALSE,col.names = TRUE,row.names = TRUE)

p.val<-fisher.test(genemat, alternative = "two.sided")$p.value #1

##PTEN##
genemat <- matrix(c(5, 7,32, 30), nrow = 2,
                  dimnames =
                    list(c("TACSTD1.low", "TACSTD1.high"),
                         c("PTEN.mut", "PTEN.nomut")))
PTEN<-as.data.frame(genemat)
write.table(PTEN,"PTEN.fisher.txt",sep="\t",quote=FALSE,col.names = TRUE,row.names = TRUE)


p.val<-fisher.test(genemat, alternative = "two.sided")$p.value #0.75

##RB1##
genemat <- matrix(c(4, 7,33, 30), nrow = 2,
                  dimnames =
                    list(c("TACSTD1.low", "TACSTD1.high"),
                         c("RB1.mut", "RB1.nomut")))
RB1<-as.data.frame(genemat)
write.table(RB1,"RB1.fisher.txt",sep="\t",quote=FALSE,col.names = TRUE,row.names = TRUE)

p.val<-fisher.test(genemat, alternative = "two.sided")$p.value #0.51



