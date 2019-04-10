#################2019-03-08######
##fisher-exact test##############
#################################
AR<-read.delim("DNA_sv_AR.txt")
samples<-as.character(unique(colnames(cn))[-1])
##replace delimiter##
samples_all<-NULL
for (i in 1:length(samples)){
  samples_all[i]<-paste0(strsplit(samples[i],"[.]")[[1]],collapse="-")
}

#negative and positive samples##
samples_pos<-as.character(AR$sample_id) #unique id 23
samples_neg<-dplyr::setdiff(samples_all,samples_pos) #78

##FUNCTION FOR ALL VARIATION FISHER EXACT TEST#################################
exact.test<-function(gene.name,gain,loss) {
cngene<-cn[cn$IDENTIFIER==gene.name,]
cngene<-cngene[,-1]
cngene_gainloss<-cngene[,which(cngene[,]>gain |cngene[,]<loss)] #39
for(i in 1:length(colnames(cngene_gainloss))){
  colnames(cngene_gainloss)[i]<-paste0(strsplit(colnames(cngene_gainloss)[i],"[.]")[[1]],collapse="-")
}

gene_cn<-colnames(cngene_gainloss)
gene_snp<-snp[snp$gene==gene.name,"sample_id"]
gene_sv<-na.omit(sv[sv$gene.1==gene.name|sv$gene.2==gene.name,"sample"])

gene_all<-unique(c(gene_cn,gene_snp,gene_sv)) #39
gene_no<-dplyr::setdiff(samples_all,gene_all) #62

var.pos<-sum(gene_all%in%samples_pos) #9
var.neg<-sum(gene_all%in%samples_neg) #30
novar.pos<-sum(gene_no%in%samples_pos)#14
novar.neg<-sum(gene_no%in%samples_neg)#48

##test##
genemat <- matrix(c(var.pos, novar.pos,var.neg, novar.neg), nrow = 2,
                  dimnames =
                    list(c("var.yes", "var.no"),
                         c("sv.pos", "sv.neg")))

p.val<-fisher.test(genemat, alternative = "two.sided")$p.value
OR<-fisher.test(genemat, alternative = "two.sided")$estimate
conf<-fisher.test(genemat, alternative = "two.sided")$conf.int
genestat<-dplyr::data_frame(gene_name=gene.name,p.value=p.val,odds.ratio=OR,CI_L=conf[1],CI_H=conf[2])

}

genes_nonsex<-list("ERG","ETV1","ETV4","ETV5","BRAF","HRAS","MYC","CHD1","SPOP","IDH1","FOXA1","NCOR1","NCOR2","ASXL2","ZBTB16","TP53","PTEN","AKT1","PIK3CA","RB1","CDKN2A","CCND1","APC",
         "CTNNB1","ZNRF3","KMT2C","KMT2D","HDAC4","NKX3-1","AXL","ZFHX3","GNAS","BRCA2","BRCA1","ATM","CDK12","ERCC2","PRKDC","MLH1","MSH2","MSH6")
genes_sex<-list("MED12","KDM6A","AR")

O_nonsex<-as.data.frame(t(mapply(exact.test,genes_nonsex,3.0,1.65)))
O_sex<-as.data.frame(t(mapply(exact.test,genes_sex,1.4,0.6)))
O_all<-as.data.frame(rbind(O_nonsex,O_sex))
df <- apply(O_all,2,as.character)
write.table(df,"var.fisher.txt",quote=FALSE,col.names = TRUE,row.names = FALSE,sep="\t")
#############################################################################################################################
##2019-04-10##################
##pvalue to FDR###############
##BH adjust or fdr############
##############################
fdr<-read.csv("var.fisher.csv",header=TRUE,sep=",")
fdr$fdr<-p.adjust(fdr$p.value,method="fdr",n=length(fdr$p.value))
##PTEN is not significant any more after adjustment.
##adj.p=1 for all genes.


##2019-03-08##
##pathway exact test#####################################################
#ETS alterations: ERG+ETV1+ETV4+ETV5
#RAS pathway: BRAF+HRAS
#PI3K pathway: PTEN+AKT1+PIK3CA
#Rb pathway: RB1+CDKN2A+CCND1
#WNT pathway: APC+CTNNB1+ZNRF3
#epigenetic regulators: KMT2C+KMT2D+KDM6A+HDAC4
#DNA repair pathway: BRCA2+BRCA1+ATM+CDK12+ERCC2+PRKDC+MLH1+MSH2+MSH6

##pathway gene statistics function##
gene.stat<-function(gene.name){
cngene<-cn[cn$IDENTIFIER==gene.name,]
cngene<-cngene[,-1]
cngene_gainloss<-cngene[,which(cngene[,]>3.0 |cngene[,]<1.65)] #39
for(i in 1:length(colnames(cngene_gainloss))){
  colnames(cngene_gainloss)[i]<-paste0(strsplit(colnames(cngene_gainloss)[i],"[.]")[[1]],collapse="-")
}

gene_cn<-colnames(cngene_gainloss)
gene_snp<-snp[snp$gene==gene.name,"sample_id"]
gene_sv<-na.omit(sv[sv$gene.1==gene.name|sv$gene.2==gene.name,"sample"])

gene_all<-paste0(unique(c(gene_cn,gene_snp,gene_sv)),collapse = ";") #39
gene_no<-paste0(dplyr::setdiff(samples_all,gene_all),collapse=";") #62
gene_info<-dplyr::data_frame(gene_name=gene.name,gene_all=gene_all,gene_no=gene_no)
}

##pathway statistics##
ETS<-list("ERG","ETV1","ETV4","ETV5")
RAS<-list("BRAF","HRAS")
PI3K<-list("PTEN","AKT1","PIK3CA")
Rb<-list("RB1","CDKN2A","CCND1")
WNT<-list("APC","CTNNB1","ZNRF3")
epigenetic<-list("KMT2C","KMT2D","KDM6A","HDAC4")
repair<-list("BRCA2","BRCA1","ATM","CDK12","ERCC2","PRKDC","MLH1","MSH2","MSH6")



###test it if it works!!!!!###############
pathway.test<-function(pathway){
O<-do.call(rbind,lapply(pathway,gene.stat))
geneall<-unique(strsplit(paste(unlist(t(O[,"gene_all"])), collapse=";"),";")[[1]]) #75
geneno<-dplyr::setdiff(samples_all,geneall) #26

var.pos<-sum(geneall%in%samples_pos) #17
var.neg<-sum(geneall%in%samples_neg) #58
novar.pos<-sum(geneno%in%samples_pos)#6
novar.neg<-sum(geneno%in%samples_neg)#20

##test##
genemat <- matrix(c(var.pos, novar.pos,var.neg, novar.neg), nrow = 2,
                  dimnames =
                    list(c("var.yes", "var.no"),
                         c("sv.pos", "sv.neg")))

p.val<-fisher.test(genemat, alternative = "two.sided")$p.value
OR<-fisher.test(genemat, alternative = "two.sided")$estimate
conf<-fisher.test(genemat, alternative = "two.sided")$conf.int
pathway.stat<-dplyr::data_frame(pathway_name=toString(pathway),p.value=p.val,odds.ratio=OR,CI_L=conf[1],CI_H=conf[2])

}
test<-pathway.test(ETS)

#######################################ETS pathway####################################
ERG<-gene.stat("ERG")
ETV1<-gene.stat("ETV1")
ETV4<-gene.stat("ETV4")
ETV5<-gene.stat("ETV5")
ETSpath<-rbind(ERG,ETV1,ETV4,ETV5)
ETSpath.gene_all<-unique(c(strsplit(ETSpath$gene_all[1],";")[[1]],strsplit(ETSpath$gene_all[2],";")[[1]],strsplit(ETSpath$gene_all[3],";")[[1]],strsplit(ETSpath$gene_all[4],";")[[1]]))
ETSpath.gene_no<-dplyr::setdiff(samples_all,ETSpath.gene_all) #26

var.pos<-sum(ETSpath.gene_all%in%samples_pos) #17
var.neg<-sum(ETSpath.gene_all%in%samples_neg) #58
novar.pos<-sum(ETSpath.gene_no%in%samples_pos)#6
novar.neg<-sum(ETSpath.gene_no%in%samples_neg)#20

##test##
genemat <- matrix(c(var.pos, novar.pos,var.neg, novar.neg), nrow = 2,
                  dimnames =
                    list(c("var.yes", "var.no"),
                         c("sv.pos", "sv.neg")))

p.val<-fisher.test(genemat, alternative = "two.sided")$p.value
OR<-fisher.test(genemat, alternative = "two.sided")$estimate
conf<-fisher.test(genemat, alternative = "two.sided")$conf.int
pathway.stat<-dplyr::data_frame(pathway_name="ETS",p.value=p.val,odds.ratio=OR,CI_L=conf[1],CI_H=conf[2])

####################RAS BRAF+HRAS
ERG<-gene.stat("BRAF")
ETV1<-gene.stat("HRAS")
ETSpath<-rbind(ERG,ETV1)
ETSpath.gene_all<-unique(c(strsplit(ETSpath$gene_all[1],";")[[1]],strsplit(ETSpath$gene_all[2],";")[[1]])) #25
ETSpath.gene_no<-dplyr::setdiff(samples_all,ETSpath.gene_all) #76

var.pos<-sum(ETSpath.gene_all%in%samples_pos) #17
var.neg<-sum(ETSpath.gene_all%in%samples_neg) #58
novar.pos<-sum(ETSpath.gene_no%in%samples_pos)#6
novar.neg<-sum(ETSpath.gene_no%in%samples_neg)#20

##test##
genemat <- matrix(c(var.pos, novar.pos,var.neg, novar.neg), nrow = 2,
                  dimnames =
                    list(c("var.yes", "var.no"),
                         c("sv.pos", "sv.neg")))

p.val<-fisher.test(genemat, alternative = "two.sided")$p.value
OR<-fisher.test(genemat, alternative = "two.sided")$estimate
conf<-fisher.test(genemat, alternative = "two.sided")$conf.int
RAS.stat<-dplyr::data_frame(pathway_name="RAS",p.value=p.val,odds.ratio=OR,CI_L=conf[1],CI_H=conf[2])

##################PI3K:PTEN+AKT1+PIK3CA
ERG<-gene.stat("PTEN")
ETV1<-gene.stat("AKT1")
ETV4<-gene.stat("PIK3CA")

ETSpath<-rbind(ERG,ETV1,ETV4)
ETSpath.gene_all<-unique(c(strsplit(ETSpath$gene_all[1],";")[[1]],strsplit(ETSpath$gene_all[2],";")[[1]],strsplit(ETSpath$gene_all[3],";")[[1]]))
ETSpath.gene_no<-dplyr::setdiff(samples_all,ETSpath.gene_all) #26

var.pos<-sum(ETSpath.gene_all%in%samples_pos) #17
var.neg<-sum(ETSpath.gene_all%in%samples_neg) #58
novar.pos<-sum(ETSpath.gene_no%in%samples_pos)#6
novar.neg<-sum(ETSpath.gene_no%in%samples_neg)#20

##test##
genemat <- matrix(c(var.pos, novar.pos,var.neg, novar.neg), nrow = 2,
                  dimnames =
                    list(c("var.yes", "var.no"),
                         c("sv.pos", "sv.neg")))

p.val<-fisher.test(genemat, alternative = "two.sided")$p.value
OR<-fisher.test(genemat, alternative = "two.sided")$estimate
conf<-fisher.test(genemat, alternative = "two.sided")$conf.int
PI3K.stat<-dplyr::data_frame(pathway_name="PI3K",p.value=p.val,odds.ratio=OR,CI_L=conf[1],CI_H=conf[2])

######################################Rb pathway: RB1+CDKN2A+CCND1

ERG<-gene.stat("RB1")
ETV1<-gene.stat("CDKN2A")
ETV4<-gene.stat("CCND1")

ETSpath<-rbind(ERG,ETV1,ETV4)
ETSpath.gene_all<-unique(c(strsplit(ETSpath$gene_all[1],";")[[1]],strsplit(ETSpath$gene_all[2],";")[[1]],strsplit(ETSpath$gene_all[3],";")[[1]]))
ETSpath.gene_no<-dplyr::setdiff(samples_all,ETSpath.gene_all) #26

var.pos<-sum(ETSpath.gene_all%in%samples_pos) #17
var.neg<-sum(ETSpath.gene_all%in%samples_neg) #58
novar.pos<-sum(ETSpath.gene_no%in%samples_pos)#6
novar.neg<-sum(ETSpath.gene_no%in%samples_neg)#20

##test##
genemat <- matrix(c(var.pos, novar.pos,var.neg, novar.neg), nrow = 2,
                  dimnames =
                    list(c("var.yes", "var.no"),
                         c("sv.pos", "sv.neg")))

p.val<-fisher.test(genemat, alternative = "two.sided")$p.value
OR<-fisher.test(genemat, alternative = "two.sided")$estimate
conf<-fisher.test(genemat, alternative = "two.sided")$conf.int
Rb.stat<-dplyr::data_frame(pathway_name="Rb",p.value=p.val,odds.ratio=OR,CI_L=conf[1],CI_H=conf[2])
###############WNT pathway: APC+CTNNB1+ZNRF3
ERG<-gene.stat("APC")
ETV1<-gene.stat("CTNNB1")
ETV4<-gene.stat("ZNRF3")

ETSpath<-rbind(ERG,ETV1,ETV4)
ETSpath.gene_all<-unique(c(strsplit(ETSpath$gene_all[1],";")[[1]],strsplit(ETSpath$gene_all[2],";")[[1]],strsplit(ETSpath$gene_all[3],";")[[1]]))
ETSpath.gene_no<-dplyr::setdiff(samples_all,ETSpath.gene_all) #26

var.pos<-sum(ETSpath.gene_all%in%samples_pos) #17
var.neg<-sum(ETSpath.gene_all%in%samples_neg) #58
novar.pos<-sum(ETSpath.gene_no%in%samples_pos)#6
novar.neg<-sum(ETSpath.gene_no%in%samples_neg)#20

##test##
genemat <- matrix(c(var.pos, novar.pos,var.neg, novar.neg), nrow = 2,
                  dimnames =
                    list(c("var.yes", "var.no"),
                         c("sv.pos", "sv.neg")))

p.val<-fisher.test(genemat, alternative = "two.sided")$p.value
OR<-fisher.test(genemat, alternative = "two.sided")$estimate
conf<-fisher.test(genemat, alternative = "two.sided")$conf.int
WNT.stat<-dplyr::data_frame(pathway_name="WNT",p.value=p.val,odds.ratio=OR,CI_L=conf[1],CI_H=conf[2])

##################epigenetic regulators: KMT2C+KMT2D+KDM6A+HDAC4
ERG<-gene.stat("KMT2C")
ETV1<-gene.stat("KMT2D")
ETV4<-gene.stat("KDM6A")
ETV5<-gene.stat("HDAC4")
ETSpath<-rbind(ERG,ETV1,ETV4,ETV5)
ETSpath.gene_all<-unique(c(strsplit(ETSpath$gene_all[1],";")[[1]],strsplit(ETSpath$gene_all[2],";")[[1]],strsplit(ETSpath$gene_all[3],";")[[1]],strsplit(ETSpath$gene_all[4],";")[[1]]))
ETSpath.gene_no<-dplyr::setdiff(samples_all,ETSpath.gene_all) #26

var.pos<-sum(ETSpath.gene_all%in%samples_pos) #17
var.neg<-sum(ETSpath.gene_all%in%samples_neg) #58
novar.pos<-sum(ETSpath.gene_no%in%samples_pos)#6
novar.neg<-sum(ETSpath.gene_no%in%samples_neg)#20

##test##
genemat <- matrix(c(var.pos, novar.pos,var.neg, novar.neg), nrow = 2,
                  dimnames =
                    list(c("var.yes", "var.no"),
                         c("sv.pos", "sv.neg")))

p.val<-fisher.test(genemat, alternative = "two.sided")$p.value
OR<-fisher.test(genemat, alternative = "two.sided")$estimate
conf<-fisher.test(genemat, alternative = "two.sided")$conf.int
epigenetic.stat<-dplyr::data_frame(pathway_name="epigenetic",p.value=p.val,odds.ratio=OR,CI_L=conf[1],CI_H=conf[2])

#####################DNA repair pathway: BRCA2+BRCA1+ATM+CDK12+ERCC2+PRKDC+MLH1+MSH2+MSH6
ERG<-gene.stat("BRCA2")
ETV1<-gene.stat("BRCA1")
ETV4<-gene.stat("ATM")
ETV5<-gene.stat("CDK12")
gene5<-gene.stat("ERCC2")
gene6<-gene.stat("PRKDC")
gene7<-gene.stat("MLH1")
gene8<-gene.stat("MSH2")
gene9<-gene.stat("MSH6")

ETSpath<-rbind(ERG,ETV1,ETV4,ETV5,gene5,gene6,gene7,gene8,gene9)
ETSpath.gene_all<-unique(c(strsplit(ETSpath$gene_all[1],";")[[1]],strsplit(ETSpath$gene_all[2],";")[[1]],strsplit(ETSpath$gene_all[3],";")[[1]],strsplit(ETSpath$gene_all[4],";")[[1]],strsplit(ETSpath$gene_all[5],";")[[1]],strsplit(ETSpath$gene_all[6],";")[[1]],strsplit(ETSpath$gene_all[7],";")[[1]],strsplit(ETSpath$gene_all[8],";")[[1]],strsplit(ETSpath$gene_all[9],";")[[1]]))
ETSpath.gene_no<-dplyr::setdiff(samples_all,ETSpath.gene_all) #26

var.pos<-sum(ETSpath.gene_all%in%samples_pos) #17
var.neg<-sum(ETSpath.gene_all%in%samples_neg) #58
novar.pos<-sum(ETSpath.gene_no%in%samples_pos)#6
novar.neg<-sum(ETSpath.gene_no%in%samples_neg)#20

##test##
genemat <- matrix(c(var.pos, novar.pos,var.neg, novar.neg), nrow = 2,
                  dimnames =
                    list(c("var.yes", "var.no"),
                         c("sv.pos", "sv.neg")))

p.val<-fisher.test(genemat, alternative = "two.sided")$p.value
OR<-fisher.test(genemat, alternative = "two.sided")$estimate
conf<-fisher.test(genemat, alternative = "two.sided")$conf.int
repair.stat<-dplyr::data_frame(pathway_name="repair",p.value=p.val,odds.ratio=OR,CI_L=conf[1],CI_H=conf[2])

df_pathway<-rbind(pathway.stat,RAS.stat,PI3K.stat,Rb.stat,WNT.stat,epigenetic.stat,repair.stat)
write.csv(df_pathway,"pathway.fisher.csv",sep=",",quote=FALSE,col.names = TRUE,row.names = TRUE)
#########################################################################################################################


##MYC gain fisher exact test##
cngene<-cn[cn$IDENTIFIER=="MYC",]
cngene<-cngene[,-1]
cngenegain<-cngene[,which(cngene[,]>3.0)] #38
for(i in 1:length(colnames(cngenegain))){
colnames(cngenegain)[i]<-paste0(strsplit(colnames(cngenegain)[i],"[.]")[[1]],collapse="-")
}

gainpos<-sum(colnames(cngenegain)%in%samples_pos) #9
gainneg<-length(colnames(cngenegain))-gainpos #29

cngenenogain<-cngene[,which(cngene[,]<=3.0)] #63
for(i in 1:length(colnames(cngenenogain))){
  colnames(cngenenogain)[i]<-paste0(strsplit(colnames(cngenenogain)[i],"[.]")[[1]],collapse="-")
}

nogainpos<-sum(colnames(cngenenogain)%in%samples_pos) #14
nogainneg<-sum(colnames(cngenenogain)%in%samples_neg) #49

genemat <- matrix(c(gainpos, nogainpos,gainneg, nogainneg), nrow = 2,
                      dimnames =
                        list(c("gainyes", "gainno"),
                             c("svpos", "svneg")))

p.val<-fisher.test(genemat, alternative = "two.sided")$p.value
OR<-fisher.test(genemat, alternative = "two.sided")$estimate
conf<-fisher.test(genemat, alternative = "two.sided")$conf.int
#names(fisher.test(genemat, alternative = "two.sided"))
#[1] "p.value"     "conf.int"    "estimate"    "null.value"  "alternative" "method"      "data.name"
genestat<-dplyr::data_frame(gene_name="gene",p.value=p.val,odds.ratio=OR,CI_L=conf[1],CI_H=conf[2])





