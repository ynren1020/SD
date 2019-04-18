#########################2019-04-18################################
##KLF5 expression boxplots (Abi/enza naïve | Abi and/or treated)###
###################################################################
library(tidyverse)
library(data.table)
library(ggpubr)
library(stringr)

treat<-read.csv("WC-SU2C-second_line_therapy.csv",sep=",",header=TRUE,stringsAsFactors = FALSE)
klf5<-read.delim("2018_04_15_matrix_rna_tpm.txt",stringsAsFactors = FALSE,header=TRUE)
klf5<-klf5[klf5$IDENTIFIER=="KLF5",]
klf5T<-transpose(klf5)
rownames(klf5T) <- colnames(klf5)
colnames(klf5T) <- rownames(klf5)
klf5T$sample_id<-rownames(klf5T)
klf5T<-klf5T[-1,]%>%rename("RPKM"="7137") #is it RPKM

for (i in 1:nrow(klf5T)){
  klf5T$sample[i]<-str_replace_all(klf5T$sample_id[i],'[.]',"-")
}
klf5T$sample_id<-NULL
klf5T<-klf5T%>%rename("sample_id"="sample")

##join##
klf5T_tr<-full_join(klf5T,treat,by="sample_id")
klf5T_tr<-na.omit(klf5T_tr)

##create status based on both naive (Naive) or not (Treated)##
for (i in 1:nrow(klf5T_tr)){
  klf5T_tr$status[i]<-ifelse(klf5T_tr$enza[i]=="naive"&klf5T_tr$abi[i]=="naive","Naive","Treated")
}

klf5T_tr$status<-as.factor(klf5T_tr$status)
klf5T_tr$RPKM<-as.numeric(klf5T_tr$RPKM)
nrow(klf5T_tr[klf5T_tr$status=="Naive",])#35
nrow(klf5T_tr[klf5T_tr$status=="Treated",])#63

##plot##
p <- ggboxplot(klf5T_tr, x = "status", y = "RPKM",
               color = "status", palette =c("#00AFBB", "#E7B800"),
               add = "jitter", shape = "status",
               xlab=FALSE,
               ylab="KLF5 expressoin (RPKM)")
p


# Add p-values comparing groups
p1<-p + stat_compare_means(label.y = 200)+ # Add p-value
  ##change tick mark labels##
p2<-p1+scale_x_discrete(labels=c("Naive" = "Naive\n(N=35)", "Treated" = "Treated\n(N=63)"))
p2

ggsave("klf5boxplot.pdf",width = 5,height = 5)
