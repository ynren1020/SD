####################################2019-04-23########################################
##junctions number in each sample#####################################################
######################################################################################

library(dplyr)
library(tidyr)

#junc<-read.delim("linear_ARV_known.txt",header=TRUE,stringsAsFactors = FALSE)


args <- commandArgs(TRUE)
dtb<-read.delim(args[1],header=FALSE,stringsAsFactors = FALSE)
dtbsub<-dtb[,c(4,5)]
dtbsub$sample<-strsplit(args[1],'[.]')[[1]][1]

write.table(dtbsub,paste0(strsplit(args[1],'[.]')[[1]][1],"juncN.txt"),quote = FALSE,sep="\t",col.names = FALSE,row.names = FALSE)

##reshape long to wide##
samples<-read.delim("samples.juncN.txt",header=FALSE,stringsAsFactors = FALSE)
samples<-samples%>%rename("ARv_Name"="V1","junction"="V2")%>%spread(V3,junction)

write.table(samples,"ARv_junc.txt",sep="\t",quote=FALSE,col.names = TRUE,row.names = FALSE)


