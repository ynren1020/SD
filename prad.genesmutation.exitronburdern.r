##########################2019-05-16###########################
##prad.maf filtered by 20genes and separate prad samples into##
##two groups: mutations or not and then compare exitron #######
##burden between two groups####################################
###############################################################

library(tidyverse)

input<-"prad.genes.mutation.txt"
mut<-read.delim(input,header=FALSE,stringsAsFactors = FALSE)
tumor_id<-unique(mut$V16) #(matched tumor barcode) #42

##/home/tywang/Projects/Exitron/TCGA/TCGA.info.txt (fileid, tumor_id)
##contains RNAseq and WGS data
input2<-"TCGA.info.txt"
info<-read.delim(input2,header=TRUE,stringsAsFactors = FALSE)
info.prad<-info[info$PROJECT=="TCGA-PRAD",]
##20genes mutation##
info.prad.mut<-info.prad[info.prad$TUMOR%in%tumor_id,]
info.prad.nomut<-setdiff(info.prad,info.prad.mut)

##fileid groups#
mutfile<-as.data.frame(paste0(info.prad.mut$FILE_ID,".exitron"))
nomutfile<-as.data.frame(paste0(info.prad.nomut$FILE_ID,".exitron"))

write.table(mutfile,"prad.mut.fileid.txt",quote=FALSE,col.names = FALSE,row.names = FALSE,sep = "\t")
write.table(nomutfile,"prad.nomut.fileid.txt",quote=FALSE,col.names = FALSE,row.names = FALSE,sep = "\t")


##compute two proportions' z test
res <- prop.test(x = c(3, 25), n = c(42, 450))
# Printing the results
res 
#2-sample test for equality of proportions with continuity correction

#data:  c(3, 25) out of c(42, 450)
#X-squared = 0.0058427, df = 1, p-value = 0.9391
#alternative hypothesis: two.sided
#95 percent confidence interval:
#  -0.07785438  0.10960041
#sample estimates:
#  prop 1     prop 2 
#0.07142857 0.05555556 
##fisher exact test##
exitron.mut<-matrix(c(3,39,25,425),nrow=2,dimnames=list(c("exitron","no-exitron"),
                                                        c("mutation","no-mutation")))

fisher.test(exitron.mut, alternative = "greater") #pval=0.4347 OR=1.306903

fisher.test(exitron.mut, alternative = "less") #p-value = 0.7912,

fisher.test(exitron.mut, alternative = "two.sided") #p-value = 0.7231 

##patients separate into two groups based on 20 genes mut or not##
##use all exitron numbers##
allexitron.mut<-read.delim("prad.mut.allexitron.txt",header=FALSE,stringsAsFactors = FALSE)
allexitron.nomut<-read.delim("prad.nomut.allexitron.txt",header=FALSE,stringsAsFactors = FALSE)

t.test(allexitron.mut$V1,allexitron.nomut$V1)

#Welch Two Sample t-test

#data:  allexitron.mut$V1 and allexitron.nomut$V1
#t = 2.2609, df = 45.875, p-value = 0.02856
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  2.437356 42.033120
#sample estimates:
#  mean of x mean of y 
#172.0952  149.8600 

#nonparametric#
wilcox.test(allexitron.mut$V1,allexitron.nomut$V1,alternative = "two.sided")

#Wilcoxon rank sum test with continuity correction

#data:  allexitron.mut$V1 and allexitron.nomut$V1
#W = 11556, p-value = 0.01687
#alternative hypothesis: true location shift is not equal to 0

##boxplot##
library(ggpubr)
allexitron.mut$group<-"DNA repair genes mutated (N=42)"
allexitron.nomut$group<-"DNA repair genes not mutated (N=450)"
df<-rbind(allexitron.mut,allexitron.nomut)


p <- ggboxplot(df, x = "group", y = "V1",
               color = "group", palette =c("#00AFBB", "#E7B800"),
               add = "jitter", shape = "group",xlab=FALSE,ylab="Exitron Load")+
               stat_compare_means(method = "t.test",label.y = 400,label = "p.format")
  
p


##############2019-05-29#################################################
##Look for unique exitrons (genes) in 20genes mutation or not group######
##find out the shared exitrons (genes) in 20genes mutation or not group##
#########################################################################
prad.mut<-read.delim("prad.mut.allsamples.exitron.0529.txt",header = FALSE,stringsAsFactors = FALSE)
prad.nomut<-read.delim("prad.nomut.allsamples.exitron.0529.txt",header = FALSE,stringsAsFactors = FALSE)

prad.mut<-unite(prad.mut,"coord",c("V1","V2","V3","V6"),sep=":")
prad.nomut<-unite(prad.nomut,"coord",c("V1","V2","V3","V6"),sep=":")

##unique corrd##
prad.mut.exitronU<-unique(prad.mut$coord) #3227
prad.nomut.exitronU<-unique(prad.nomut$coord) #14413
 
shared.exitron<-intersect(prad.mut.exitronU,prad.nomut.exitronU) #2167
prad.mut.exitronUfinal<-setdiff(prad.mut.exitronU,shared.exitron) #1060
prad.nomut.exitronUfinal<-setdiff(prad.nomut.exitronU,shared.exitron) #12246
write.table(shared.exitron,"prad.mutNomut.sharedexitron.txt",quote = FALSE,col.names = FALSE,row.names = FALSE,sep="\t")
write.table(prad.mut.exitronUfinal,"prad.mut.exitronUnique.txt",quote = FALSE,col.names = FALSE,row.names = FALSE,sep="\t")
write.table(prad.nomut.exitronUfinal,"prad.nomut.exitronUnique.txt",quote = FALSE,col.names = FALSE,row.names = FALSE,sep="\t")


##unique genes##
prad.mut.genes<-unique(prad.mut$V7) #1760
prad.nomut.genes<-unique(prad.nomut$V7) #4806

shared.genes<-intersect(prad.mut.genes,prad.nomut.genes) #1545
prad.mut.genesU<-setdiff(prad.mut.genes,shared.genes) #215
prad.nomut.genesU<-setdiff(prad.nomut.genes,shared.genes) #3261


write.table(shared.genes,"prad.mutNomut.sharedgenes.txt",quote = FALSE,col.names = FALSE,row.names = FALSE,sep="\t")
write.table(prad.mut.genesU,"prad.mut.genesUnique.txt",quote = FALSE,col.names = FALSE,row.names = FALSE,sep="\t")
write.table(prad.nomut.genesU,"prad.nomut.genesUnique.txt",quote = FALSE,col.names = FALSE,row.names = FALSE,sep="\t")


