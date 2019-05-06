##########################2019-05-03############################
##Normalized AR and AR-V7 expression boxplot####################

#Samples with Gene Gain + Enhancer Gain (n=48)
#DTB-024-PRO     DTB-035-BL         DTB-037-BL         DTB-059-BL         DTB-061-BL         DTB-064-BL         DTB-077-PRO     DTB-085-BL                DTB-097-PRO     DTB-100-BL         DTB-111-PRO     DTB-112-BL         DTB-119-PRO     DTB-126-BL         DTB-137-PRO                DTB-140-BL         DTB-146-BL         DTB-156-BL         DTB-167-PRO     DTB-170-BL         DTB-172-BL         DTB-173-BL                DTB-190-BL         DTB-202-BL         DTB-216-PRO     DTB-220-BL         DTB-232-PRO     DTB-265-PRO     DTB-060-BL                DTB-071-BL         DTB-080-BL         DTB-091-BL         DTB-128-BL         DTB-149-BL         DTB-151-BL         DTB-183-BL                DTB-205-BL         DTB-206-BL         DTB-214-BL         DTB-223-BL         DTB-252-BL         DTB-255-BL         DTB-258-BL                DTB-266-BL       DTB-019-PRO     DTB-034-BL         DTB-138-BL         DTB-213-BL
#Samples with Gene Gain + Enhancer Gain + AR-SV (n=20)
#DTB-009-BL DTB-023-BL DTB-069-BL DTB-089-BL  DTB-090-PRO DTB-092-BL  DTB-094-BL DTB-101-BL  DTB-165-PRO  DTB-175-BL  DTB-186-BL   DTB-201-PRO  DTB-005-BL  DTB-018-BL   DTB-063-BL    DTB-141-BL    DTB-143-BL    DTB-176-BL    DTB-102-PRO     DTB-104-BL
#Samples with no AR Alteration (n= 13)
#DTB-036-BL DTB-042-BL DTB-053-BL  DTB-055-PRO DTB-067-PRO DTB-074-BL  DTB-098-PRO2   DTB-135-PRO  DTB-194-PRO     DTB-222-BL   DTB-040-BL   DTB-159-BL   DTB-204-BL



library(tidyverse)
library(data.table)
library(ggpubr)
library(stringr)

dat<-read.delim("2018_04_15_matrix_rna_tpm.txt",stringsAsFactors = FALSE,header=TRUE)
AR.dat<-dat[dat$IDENTIFIER=="AR",]

AR.datT<-transpose(AR.dat)
rownames(AR.datT) <- colnames(AR.dat)
colnames(AR.datT) <- rownames(AR.dat)
AR.datT$sample_id<-rownames(AR.datT)
AR.datT<-AR.datT[-1,]%>%rename("RPKM"="813") #is it RPKM? it is TPM

for (i in 1:nrow(AR.datT)){
  AR.datT$sample[i]<-str_replace_all(AR.datT$sample_id[i],'[.]',"-")
}
AR.datT$sample_id<-NULL
AR.datT<-AR.datT%>%rename("sample_id"="sample")

##groups##
##Samples with Gene Gain + Enhancer Gain:48
gene.enhancer<-c('DTB-024-PRO', 'DTB-035-BL', 'DTB-037-BL', 'DTB-059-BL', 'DTB-061-BL', 'DTB-064-BL',  'DTB-077-PRO', 'DTB-085-BL', 'DTB-097-PRO','DTB-100-BL',  'DTB-111-PRO', 'DTB-112-BL', 'DTB-119-PRO',  'DTB-126-BL', 'DTB-137-PRO','DTB-140-BL', 'DTB-146-BL', 'DTB-156-BL','DTB-167-PRO',     'DTB-170-BL', 'DTB-172-BL', ' DTB-173-BL', 'DTB-190-BL','DTB-202-BL','DTB-216-PRO','DTB-220-BL', 'DTB-232-PRO',' DTB-265-PRO', 'DTB-060-BL', 'DTB-071-BL','DTB-080-BL', 'DTB-091-BL', 'DTB-128-BL', 'DTB-149-BL', 'DTB-151-BL', 'DTB-183-BL', 'DTB-205-BL', 'DTB-206-BL',  'DTB-214-BL',  'DTB-223-BL',  'DTB-252-BL',  'DTB-255-BL', 'DTB-258-BL', 'DTB-266-BL', 'DTB-019-PRO','DTB-034-BL', 'DTB-138-BL', 'DTB-213-BL')
#Samples with Gene Gain + Enhancer Gain + AR-SV:20
gene.enhancer.arsv<-c('DTB-009-BL', 'DTB-023-BL', 'DTB-069-BL', 'DTB-089-BL',  'DTB-090-PRO', 'DTB-092-BL',  'DTB-094-BL', 'DTB-101-BL',  'DTB-165-PRO',  'DTB-175-BL',  'DTB-186-BL',   'DTB-201-PRO',  'DTB-005-BL',  'DTB-018-BL',   'DTB-063-BL',    'DTB-141-BL',    'DTB-143-BL',    'DTB-176-BL',    'DTB-102-PRO',     'DTB-104-BL')
#Samples with no AR Alteration:13
noarsv<-c('DTB-036-BL', 'DTB-042-BL', 'DTB-053-BL',  'DTB-055-PRO', 'DTB-067-PRO', 'DTB-074-BL',  'DTB-098-PRO2',   'DTB-135-PRO',  'DTB-194-PRO',     'DTB-222-BL',   'DTB-040-BL',   'DTB-159-BL',   'DTB-204-BL')

samples<-c(gene.enhancer,gene.enhancer.arsv,noarsv)
samples<-as.data.frame(samples)

for (i in 1:nrow(samples)){
  samples$group[i]<-ifelse(samples$samples[i]%in%gene.enhancer,"Gene Gain + Enhancer Gain",ifelse(samples$samples[i]%in%gene.enhancer.arsv,"Gene Gain + Enhancer Gain + AR-SV",ifelse(samples$samples[i]%in%noarsv,"no AR Alteration",NA)))
}

AR.datT<-full_join(AR.datT,samples,by=c("sample_id"="samples"))

write.table(AR.datT,"AR.datT.group.txt",quote=FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")

##calculate DTB-173-BL,DTB-265-PRO,DTB-063-BL samples'TPM##
##1)Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK)
##2)Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
##3)Divide the RPK values by the “per million” scaling factor. This gives you TPM.

##read in data##
#ar.count<-read.delim("AR.header.txt",header = TRUE,stringsAsFactors = FALSE)
dtb<-read.delim("DTB_name_counts.txt",header = TRUE,stringsAsFactors = FALSE)
dtb$DTB.173.BL.rpk<-with(dtb,X.home.ryang.ddn.ucsf_crpc.AppResults.DTB.173.BL.T.RNA.Files.alignments.DTB.173.BL.T.RNA.alignments.bam/Length*1000)
dtb$DTB.265.Pro.rpk<-with(dtb,X.home.ryang.ddn.ucsf_crpc.AppResults.DTB.265.Pro.T.RNA.Files.alignments.DTB.265.Pro.T.RNA.alignments.bam/Length*1000)
dtb$DTB.063.BL.rpk<-with(dtb,X.home.ryang.ddn.ucsf_crpc.AppResults.DTB.063.BL.NS.T.RNA.Files.alignments.DTB.063.BL.NS.T.RNA.alignments.bam/Length*1000)

dtb$DTB.173.BL.tpm<-dtb$DTB.173.BL.rpk/sum(dtb$DTB.173.BL.rpk)*10^6
dtb$DTB.265.Pro.tpm<-dtb$DTB.265.Pro.rpk/sum(dtb$DTB.265.Pro.rpk)*10^6
dtb$DTB.063.BL.tpm<-dtb$DTB.063.BL.rpk/sum(dtb$DTB.063.BL.rpk)*10^6

#AR tpm#
ar.tpm<-dtb[dtb$Geneid=="AR",c(1,13:15)]
#Geneid DTB.173.BL.tpm DTB.265.Pro.tpm DTB.063.BL.tpm
#55216     AR       1040.978        2183.542       1841.138

AR.datT$RPKM[100]<-1040.978
AR.datT$RPKM[101]<-2183.542
AR.datT$RPKM[102]<-1841.138

AR.datT<-na.omit(AR.datT)
nrow(AR.datT[AR.datT$group=="Gene Gain + Enhancer Gain",]) #48
nrow(AR.datT[AR.datT$group=="Gene Gain + Enhancer Gain + AR-SV",]) #20
nrow(AR.datT[AR.datT$group=="no AR Alteration",]) #13


AR.datT$group<-as.factor(AR.datT$group)
AR.datT$RPKM<-as.numeric(AR.datT$RPKM)
AR.datT$logRPKM<-log10(AR.datT$RPKM)


##plot##
p <- ggboxplot(AR.datT, x = "group", y = "logRPKM",
               color = "group", palette =c("#00AFBB", "#E7B800","#FC4E07"),
               add = "jitter", shape = "group",
               xlab=FALSE,
               ylab="AR expressoin (logTPM)")
p


# Add p-values comparing groups
my_comparisons <- list( c("Gene Gain + Enhancer Gain", "Gene Gain + Enhancer Gain + AR-SV"), c("Gene Gain + Enhancer Gain + AR-SV", "no AR Alteration"), c("no AR Alteration", "Gene Gain + Enhancer Gain") )
p1<-p + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 6) +
  rotate_x_text(90)

##change tick mark labels##
p2<-p1+scale_x_discrete(labels=c("Gene Gain + Enhancer Gain" = "Gene Gain\nEnhancer Gain\n(N=48)", "Gene Gain + Enhancer Gain + AR-SV" = "Gene Gain\nEnhancer Gain\nAR-SV\n(N=20)",
                                 "no AR Alteration" = "no AR Alteration\n(N=13)"))
p2

ggsave("ARboxplot.pval.pdf",width = 8,height = 8)



###############AR-V7 normalized junctions for boxplot#######################
arv7<-read.delim("AR-V7.juncN.txt",header = FALSE,stringsAsFactors = FALSE)
arv7<-rename(arv7,"sample_id"="V2","AR-V7"="V1") #109

total<-read.delim("sample.total.txt",header = TRUE,stringsAsFactors = FALSE)
totalT<-transpose(total)
rownames(totalT) <- colnames(total)
colnames(totalT) <- rownames(total)
totalT$sample_id<-rownames(totalT)
totalT<-totalT[-1,]%>%rename("total"="1") 

for (i in 1:nrow(totalT)){
  totalT$sample[i]<-str_replace_all(totalT$sample_id[i],'[.]',"-")
}

totalT$sample_id<-NULL
totalT<-totalT%>%rename("sample_id"="sample") #168
totalT$sample_id<-toupper(totalT$sample_id)

##join totalT and arv7##

arv7.total<-full_join(arv7,totalT,by="sample_id")
arv7.total<-na.omit(arv7.total) #109
#write.table(arv7.total,"arv7.total.txt",quote = FALSE,col.names = FALSE,row.names = FALSE,sep = "\t") #for match with arv7_sample.txt, get arv7.total.match.txt for boxplot

#arv7.match<-read.delim("arv7.total.match.txt",header = FALSE,stringsAsFactors = FALSE)
#arv7.match<-rename(arv7.match,"sample_id"="V2","AR-V7"="V1","total"="V3")

arv7.total$sample<-str_replace(arv7.total$sample_id,"-T-RNA", "")
arv7.total$sample<-str_replace(arv7.total$sample,"-NS", "")
arv7.total$sample<-str_replace(arv7.total$sample,"-D", "")
arv7.total$sample_id<-NULL
arv7.total<-arv7.total%>%rename("sample_id"="sample")
arv7.total<-arv7.total[!duplicated(arv7.total$sample_id), ]
samples$samples<-as.character(samples$samples)
arv7.total<-full_join(arv7.total,samples,by=c("sample_id"="samples"))
arv7.total[arv7.total$sample_id=="DTB-173-BL","group"]<-"Gene Gain + Enhancer Gain"
arv7.total[arv7.total$sample_id=="DTB-265-PRO","group"]<-"Gene Gain + Enhancer Gain"
#DTB-173-BL 28
#DTB-265-Pro 42

#for (i in 1:nrow(arv7.total)){
#  arv7.total$group[i]<-ifelse(arv7.total$sample_id[i]%in%gene.enhancer,"Gene Gain + Enhancer Gain",ifelse(arv7.total$sample_id[i]%in%gene.enhancer.arsv,"Gene Gain + Enhancer Gain + AR-SV",ifelse(arv7.total$sample_id[i]%in%noarsv,"no AR Alteration",NA)))
#}

arv7.total<-na.omit(arv7.total)


nrow(arv7.total[arv7.total$group=="Gene Gain + Enhancer Gain",]) #42
nrow(arv7.total[arv7.total$group=="Gene Gain + Enhancer Gain + AR-SV",]) #16
nrow(arv7.total[arv7.total$group=="no AR Alteration",]) #4

##normalize junction##
arv7.total$arv7_N<-with(arv7.total,`AR-V7`/total*10^6)
arv7.total$group<-as.factor(arv7.total$group)

#################for id match####
arv7_sample<-c(gene.enhancer,gene.enhancer.arsv,noarsv)
arv7_sample<-as.data.frame(arv7_sample)
write.table(arv7_sample,"arv7_sample.txt",quote = FALSE,col.names = FALSE,row.names = FALSE,sep = "\t")
arv7$sample_id<-toupper(arv7$sample_id)
write.table(arv7,"AR-V7.juncN.txt",quote = FALSE,col.names = FALSE,row.names = FALSE,sep = "\t")
#grep -f arv7_sample.txt AR-V7.juncN.txt >AR-V7.juncN.match.txt


###boxplot#######

p <- ggboxplot(arv7.total, x = "group", y = "arv7_N",
               color = "group", palette =c("#00AFBB", "#E7B800","#FC4E07"),
               add = "jitter", shape = "group",
               xlab=FALSE,
               ylab="AR-V7 Normalized Junction Numbers")
p


# Add p-values comparing groups
my_comparisons <- list( c("Gene Gain + Enhancer Gain", "Gene Gain + Enhancer Gain + AR-SV"), c("Gene Gain + Enhancer Gain + AR-SV", "no AR Alteration"), c("no AR Alteration", "Gene Gain + Enhancer Gain") )
p1<-p + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 22) +
  rotate_x_text(90)

##change tick mark labels##
p2<-p1+scale_x_discrete(labels=c("Gene Gain + Enhancer Gain" = "Gene Gain\nEnhancer Gain\n(N=42)", "Gene Gain + Enhancer Gain + AR-SV" = "Gene Gain\nEnhancer Gain\nAR-SV\n(N=16)",
                                 "no AR Alteration" = "no AR Alteration\n(N=4)"))
p2

ggsave("AR-V7boxplot.pval.pdf",width = 8,height = 8)



