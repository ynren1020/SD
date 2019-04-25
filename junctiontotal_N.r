####################################2019-04-23########################################
##total junctions number in each sample for normalization#############################
######################################################################################

#library(dplyr)
#library(tidyr)

#junc<-read.delim("DTB-266-BL-T-RNA.janno",header=TRUE,stringsAsFactors = FALSE)
#junctotal<-as_tibble(sum(junc[,5]))
#colnames(junctotal)<-strsplit("DTB-266-BL-T-RNA.janno",'[.]')[[1]][1]

#write.table(junctotal,paste0(strsplit("DTB-266-BL-T-RNA.janno",'[.]')[[1]][1],"total.txt"),quote = FALSE,sep="\t",col.names = TRUE,row.names = FALSE)

args <- commandArgs(TRUE)
dtb<-read.delim(args[1],header=TRUE,stringsAsFactors = FALSE)
dtbsub<-as.data.frame(sum(dtb[,5]))
colnames(dtbsub)<-strsplit(args[1],'[.]')[[1]][1]

write.table(dtbsub,paste0(strsplit(args[1],'[.]')[[1]][1],".total.txt"),quote = FALSE,sep="\t",col.names = TRUE,row.names = FALSE)

##normalize ARv_junc.txt with sample.total.final.txt##
ar<-read.delim("ARv_junc.txt",header=TRUE,stringsAsFactors = FALSE)
total<-read.delim("sample.total.final.txt",header = TRUE,stringsAsFactors = FALSE)

totalsub<-total[,colnames(total)%in%colnames(ar)]

##combine##
ar_totalsub<-rbind(ar,totalsub)
test<-ar_totalsub[,2:ncol(ar_totalsub)]
rownames(test)<-ar_totalsub$ARv_Name
##normlized by total in row10##
for (i in 1:ncol(test)){
  test[,i]<-(test[,i]/test[10,i])*(10^6)
  
}

##remove 10throw##
test<-test[-10,]
##add column##
test$ARv_Name<-rownames(test)
##move ARv_Name column to first column##
test<-test%>%select("ARv_Name",everything())

write.table(test,"ARv_junc_normalize.txt",col.names = TRUE,row.names = FALSE,quote=FALSE,sep="\t")








