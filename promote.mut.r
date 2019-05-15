#########################2019-05-15#################
##promote data mutation info########################
####################################################

library(tidyverse)

input2<-"tier2.info.txt"
tier2<-read.delim(input2,header = TRUE,stringsAsFactors = FALSE)
tier2<-tier2[,c("exn","sample_id")]
tier2$exn<-toupper(tier2$exn)

input<-"promote.XRCC2.mutation.txt"
atm<-read.delim(input,header = FALSE,stringsAsFactors = FALSE)
for (i in 1:nrow(atm)){
  atm$exn[i]<-strsplit(strsplit(atm$V1[i],"[.]")[[1]][1],"_")[[1]][2]
}

atm.tier2<-full_join(atm,tier2,by="exn")
atm.tier2<-atm.tier2%>%filter(!is.na(atm.tier2$V1))

mutation<-data.frame(Sample=unique(atm.tier2$sample_id),Gene="Mutation",Alteration="HOMDEL",Type="CNA")

write.table(mutation,toupper(input),sep = "\t",quote = FALSE,col.names = TRUE,row.names = FALSE)
                     