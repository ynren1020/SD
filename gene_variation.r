########2019-02-26########
##copy number,SNP, SV ####
##frequency###############

library(dplyr)
##load data##
cn<-read.delim("2018_04_15_matrix_CN_integer_symbol_copycat.txt",stringsAsFactors = FALSE)
snp<-read.delim("DriverGenes-1-s2.0-S0092867418308420-mmc5 (1).txt",stringsAsFactors = FALSE)
snp1<-snp[1:187,]
snp2<-snp[188:260,]
snp2<-dplyr::rename(snp2,sample_id=gene,gene=sample_id)
snp2<-snp2[,c(2,1,3:13)]
snp<-rbind(snp1,snp2)
sv<-read.delim("supplementary_table_S3_CN+SV_drivergenes.txt",stringsAsFactors = FALSE)

##interested gene##
genes<-c("ERG","ETV1","ETV4","ETV5","BRAF","HRAS","MYC","CHD1","SPOP","IDH1","AR","AR enhancer","FOXA1","NCOR1","NCOR2","ASXL2","ZBTB16","TP53","PTEN","AKT1","PIK3CA","RB1","CDKN2A","CCND1","APC",
         "CTNNB1","ZNRF3","KMT2C","KMT2D","KDM6A","HDAC4","NKX3-1","AXL","MED12","ZFHX3","GNAS","BRCA2","BRCA1","ATM","CDK12","ERCC2","PRKDC","MLH1","MSH2","MSH6")




##subset dataset based on interested genes##
cnsub<-cn[cn$IDENTIFIER%in%genes,]
cnsub$gene<-cnsub$IDENTIFIER
##chr info for cnv filter##
snpinfo<-snp[,c("gene","chromosome")]
snpinfo<-snpinfo[!duplicated(snpinfo),]

##NA leads to trouble, so make chromosome NA free###
cnsubnew<-left_join(cnsub,snpinfo,by="gene") ##include chromosome info for threthold filter cnv
cnsubnew[cnsubnew$IDENTIFIER=="AXL","chromosome"]<-"chr19"
cnsubnew[cnsubnew$IDENTIFIER=="CCND1","chromosome"]<-"chr11"
cnsubnew[cnsubnew$IDENTIFIER=="CHD1","chromosome"]<-"chr5"
cnsubnew[cnsubnew$IDENTIFIER=="HDAC4","chromosome"]<-"chr2"
cnsubnew[cnsubnew$IDENTIFIER=="MYC","chromosome"]<-"chr8"
cnsubnew[cnsubnew$IDENTIFIER=="NKX3-1","chromosome"]<-"chr8"
cnsubnew[cnsubnew$IDENTIFIER=="PIK3CA","chromosome"]<-"chr3"
cnsubnew[cnsubnew$IDENTIFIER=="ZNRF3","chromosome"]<-"chr22"
cnsubnew[cnsubnew$IDENTIFIER=="ETV1","chromosome"]<-"chr7"
cnsubnew[cnsubnew$IDENTIFIER=="ETV4","chromosome"]<-"chr17"


##snp sv subset##
snpsub<-snp[snp$gene%in%genes,]
svsub<-sv[sv$gene.1%in%genes|sv$gene.2%in%genes,]

##how many unique samples##
length(unique(colnames(cn))) #102 cn
length(unique(snp$sample_id)) #102 snp
length(unique(sv$sample))     #89  sv

samples<-unique(colnames(cn))[-1] ##this is the unique sample## 
#samples2<-as.character(unique(snp$sample_id))
#samples3<-as.character(unique(sv$sample))
#samples_uniq<-unique(c(samples,samples2,samples3)) #200 ##samples were used different sep, (. and -)




##################################################################################################################
cnAR<-cn[cn$IDENTIFIER=="AR",]
cnARcngain<-cnAR[,which(cnAR[,]>1.5)]
ncol(cnARcngain) #64

cnMYC<-cn[cn$IDENTIFIER=="MYC",]
cnMYCcngain<-cnMYC[,which(cnMYC[,]>3.10)]
ncol(cnMYCcngain) #34

cnETV5<-cn[cn$IDENTIFIER=="ETV5",]
cnETV5cngain<-cnETV5[,which(cnETV5[,]>3.8)]
ncol(cnETV5cngain) #2

cnETV1<-cn[cn$IDENTIFIER=="ETV1",]
cnETV1cngain<-cnETV1[,which(cnETV1[,]>3.8)]
ncol(cnETV1cngain) #1

cnCCND1<-cn[cn$IDENTIFIER=="CCND1",]
cnCCND1cngain<-cnCCND1[,which(cnCCND1[,]>3.10)]
ncol(cnCCND1cngain) #8

cnFOXA1<-cn[cn$IDENTIFIER=="FOXA1",]
cnFOXA1cngain<-cnFOXA1[,which(cnFOXA1[,]>3.0)]
ncol(cnFOXA1cngain) #8

cnCTNNB1<-cn[cn$IDENTIFIER=="CTNNB1",]
cnCTNNB1cngain<-cnCTNNB1[,which(cnCTNNB1[,]>3.0)]
ncol(cnCTNNB1cngain) #3

cnBRAF2<-cn[cn$IDENTIFIER=="BRAF",]
cnBRAF2cngain<-cnBRAF2[,which(cnBRAF2[,]>3.3)]
ncol(cnBRAF2cngain) #2



##loss##
cn[cn$IDENTIFIER=="TP53",which(cn[cn$IDENTIFIER=="TP53",]<1.45)] #43
cn[cn$IDENTIFIER=="PTEN",which(cn[cn$IDENTIFIER=="PTEN",]<1.45)] #42
cn[cn$IDENTIFIER=="RB1",which(cn[cn$IDENTIFIER=="RB1",]<1.0)] #16
cn[cn$IDENTIFIER=="BRCA2",which(cn[cn$IDENTIFIER=="BRCA2",]<1.0)] #8
cn[cn$IDENTIFIER=="ATM",which(cn[cn$IDENTIFIER=="ATM",]<1.21)] #3
cn[cn$IDENTIFIER=="HDAC4",which(cn[cn$IDENTIFIER=="HDAC4",]<1.1)] #4
cn[cn$IDENTIFIER=="APC",which(cn[cn$IDENTIFIER=="APC",]<1.2)] #5
cn[cn$IDENTIFIER=="CHD1",which(cn[cn$IDENTIFIER=="CHD1",]<1.12)] #7 #<1.0:5; 
cn[cn$IDENTIFIER=="NCOR1",which(cn[cn$IDENTIFIER=="NCOR1",]<1.12)]
##Test############################################################################################################
cn[cn$IDENTIFIER=="KDM6A",which(cn[cn$IDENTIFIER=="KDM6A",]>1.4)]
cn[cn$IDENTIFIER=="MED12",which(cn[cn$IDENTIFIER=="MED12",]<0.6)]
cn[cn$IDENTIFIER=="AXL",which(cn[cn$IDENTIFIER=="AXL",]<1.65)]
cn[cn$IDENTIFIER=="PRKDC",which(cn[cn$IDENTIFIER=="PRKDC",]>3.0)]

############################################################################################################################

#####################
##USE THE THRETHOLD FOR COPY NUMBER GAIN and loss##
GAIN_NONSEX = 3
LOSS_SINGLE_NONSEX = 1.65
LOSS_DOUBLE_NONSEX = 0.6
GAIN_SEX = 1.4
LOSS_SEX = 0.6
##gain##
cnsubnewsex<-cnsubnew[cnsubnew$chromosome=="chrX",]
cnsubnewnonsex<-cnsubnew[cnsubnew$chromosome!="chrX",]



#cnsubstat<-cnsub[,1:102]
#rownames(cnsubstat)<-cnsubstat$IDENTIFIER
#cnsubstat$IDENTIFIER<-NULL
##how many cn snp and sv##
##copy number##
#cnsubwoAR<-cnsub[cnsub$IDENTIFIER!="AR",]
##USE THIS FOR LOOP for cn_freq statistics#################################
for (i in 1:nrow(cnsubnew)){
  cnsubnew_num<-cnsubnew[,2:102]
 ifelse(cnsubnew$chromosome[i]=="chrX",cnsubnew$sum[i]<-length(cnsubnew_num[i,which(cnsubnew_num[i,]>1.4)]),cnsubnew$sum[i]<-length(cnsubnew_num[i,which(cnsubnew_num[i,]>3.0)]))
  cnsubnew$homo[i]<-length(cnsubnew_num[i,which(cnsubnew_num[i,]<0.6)])
 ifelse(cnsubnew$chromosome[i]!="chrX",cnsubnew$heter[i]<-length(cnsubnew_num[i,which(0.60<cnsubnew_num[i,]&cnsubnew_num[i,]<1.65)]),cnsubnew$heter[i]<-NA)
  }

cn_freq<-cnsubnew[,c("IDENTIFIER","sum","homo","heter")]%>%
  rename(gene=IDENTIFIER)%>%
  rename(cngain=sum)
 

##snp##
snp_freq<-snpsub %>%
  group_by(gene) %>%
  summarise(snpcount=n_distinct(sample_id)) #unique sample ?#  not uniuqe use =n()

##sv##
##two genes are both exist##
svsub1<-na.omit(svsub)
svsub1_gene1<-svsub1%>%
  group_by(gene.1)%>%
  summarise(svcount=n_distinct(sample))%>%
  rename(gene=gene.1)
  
svsub1_gene2<-svsub1%>%
  group_by(gene.2)%>%
  summarise(svcount=n_distinct(sample))%>%
  rename(gene=gene.2)
svsub1_freq<-rbind(svsub1_gene1,svsub1_gene2)
svsub1_freq<-unique(svsub1_freq)
svsub1_freq[svsub1_freq$gene=="FOXA1","svcount"]<-2

##only gene.2 exist##
svsub2<-svsub[is.na(svsub$gene.1),]
svsub2_freq<-svsub2 %>%
  group_by(gene.2) %>%
  summarise(svcount=n_distinct(sample)) %>%
  rename(gene=gene.2)
##sv frequency altogher##
sv_freq<-full_join(svsub1_freq,svsub2_freq,by="gene")
sv_freq[is.na(sv_freq)]<-0
sv_freq$svcount<-sv_freq$svcount.x+sv_freq$svcount.y
sv_freq$svcount.x<-sv_freq$svcount.y<-NULL

##gene variation frequency altogether##
gene_freq<-full_join(cn_freq,snp_freq,by="gene")
gene_freq_all<-full_join(gene_freq,sv_freq,by="gene")
gene_freq_all[is.na(gene_freq_all)]<-0

write.csv(gene_freq_all,"gene_freq_all_variant.csv",sep=",",row.names = FALSE,col.names = TRUE,quote=FALSE)
