####################2019-08-05##################
##DEG analysis for cohort studies###############
##remove batch effect by sva  ##################
################################################

BiocManager::install("sva")
BiocManager::install("zebrafishRNASeq")

library("edgeR")
library("sva")
#library("R.utils")
library(tidyverse)
library(purrr)
library(limma)




# args <- commandArgs(trailingOnly = TRUE);
# input.matrix = args[1];
input1<-"PRAD.tumors.counts.txt"
input2<-"PRAD.normals.counts.txt"
input3<-"Meta.counts.txt"
input4<-"gene_names_4sva.txt"
#input<-"normal.tumor.meta.txt"
tumor<-read.delim(input1,check.names=FALSE,stringsAsFactors = FALSE); #61484 502
normal<-read.delim(input2, stringsAsFactors = FALSE, check.names=FALSE); #61483 55
meta<-read.delim(input3, stringsAsFactors = FALSE, check.names=FALSE); #61484 73
genename<-read.delim(input4, stringsAsFactors = FALSE, check.names=FALSE)

genename.sub<-filter(genename,Category=="protein_coding")
tumor.pc<-tumor[tumor$Geneid%in%genename.sub$GeneID,] #19940   502
normal.pc<-normal[normal$Geneid%in%genename.sub$GeneID,] #19940    54
meta.pc<-meta[meta$Geneid%in%genename.sub$GeneID,] #19940    73

#all<-read.delim(input, stringsAsFactors = FALSE, check.names=FALSE)
#prad<-full_join(normal,meta, by = c("Geneid","Length"))
prad<-list(normal.pc,tumor.pc,meta.pc)%>%reduce(full_join,by = c("Geneid","Length")) #19940   625

#prad<-na.omit(prad) #59895   625
group<-as.factor(c(rep(1,500),rep(2,52),rep(3,71)))
row.names(prad)<-prad$Geneid

#src.data <- read.delim("./tumor_gtex.matrix",row.names="Geneid",check.names=FALSE);
#src.data <- read.delim(input.matrix, row.names="Geneid", check.names=FALSE);
data <- prad[,-c(1,2)];

#contract <- read.delim("../tumor_gtex.contract", check.names=FALSE, sep="\t", header=T);

#length(which(contract$Class ==1)); # 90
# 90 tumors samples and 118 gtex samples
#contract$condition = as.factor(ifelse(contract$Class == 1, 'tumor', 'normal')); # condition is factor in R

y <- DGEList(counts=data, group=group);

# filtering out low expressed genes
#filter = apply(zfGenes, 1, function(x) length(x[x>5])>=2)
#filtered = zfGenes[filter,]


keep <- rowSums(cpm(y)>1) >= 52
table(keep);
# FALSE  TRUE 
# 4587 15353 
y <- y[keep, keep.lib.sizes=FALSE];

# Note that the filtering does not use knowledge of what treatment corresponds to each sample, so
# the filtering does not bias the subsequent differential expression analysis.
# The TMM normalization is applied to account for the compositional biases:


# normalization by the library sizes
y <- calcNormFactors(y);
#y$samples;


# colors <- rep(c("blue","red"), 2)
# plotMDS(y, col=colors[contract$condition])


#plotMD(cpm(y, log=TRUE), column=1)
#abline(h=0, col="red", lty=2, lwd=2)

tmm <- cpm(y$counts);

write.table(tmm, file=paste0('normal.tumor.met.scott.txt'), sep='\t', row.names = TRUE, quote = FALSE);
# Inputing RNA-seq counts to clustering or heatmap routines designed for microarray data is not
# straight-forward, and the best way to do this is still a matter of research. To draw a heatmap
# of individual RNA-seq samples, we suggest using moderated log-counts-per-million. This can be
# calculated by cpm with positive values for prior.count, for example

#logcpm <- cpm(y, prior.count=1, log=TRUE);
#write.table(logcpm, file=paste0(args[1],'.tumor_gtex.logcpm.txt'), sep='\t',row.names = TRUE, quote = FALSE);
# where y is the normalized DGEList object. This produces a matrix of log2 counts-per-million
# (logCPM), with undefined values avoided and the poorly defined log-fold-changes for low counts
# shrunk towards zero. Larger values for prior.count produce stronger moderation of the values
# for low counts and more shrinkage of the corresponding log-fold-changes.


# SVA analysis
mod  <- model.matrix(~group);
mod0 <- cbind(mod[,1]);

##use "be" number of sv reduced to 3 from 619 by leek
#n.sv = num.sv(tmm,mod,method="be") #3
n.sv = num.sv(y$counts,mod,method="be") #3
#####################################################
svobj <- svaseq(y$counts, mod, mod0,n.sv = 3);
# svobj$n.sv == 26; $n.sv = 29 if no low expression filtering provided

svs<-svobj$sv;
colnames(svs) <- paste0("sv", 1:svobj$n.sv);

# design <- cbind(svs, mod); # ????

design.data <- as.data.frame(cbind(svs, group));
# the condition must be a factor
design.data$group <- as.factor(design.data[, svobj$n.sv + 1]);

rownames(design.data) <- colnames(y);

express <- paste(c(paste(colnames(svs), collapse="+"),"group"), collapse = "+");
expression <- as.formula(paste("~",express));

design <- model.matrix(expression, data=design.data);


##use limma package##
v <- voom(y,design,plot = TRUE)
# extract the normalised read counts
counts.voom <- v$E

# save normalised expression data into output dir
write.table(counts.voom,file="counts.voom.txt",row.names=T,quote=F,sep="\t");

# fit linear model for each gene given a series of libraries
fit <- lmFit(v, design)
# construct the contrast matrix corresponding to specified contrasts of a set of parameters
matrix.2vs1 <- makeContrasts(group2,levels=design)
# compute estimated coefficients and standard errors for a given set of contrasts
fit.2vs1 <- contrasts.fit(fit, matrix.2vs1)

# compute moderated t-statistics of differential expression by empirical Bayes moderation of the standard errors
fit.2vs1 <- eBayes(fit.2vs1)
options(digits=3)

#colnames(fit$coefficients)

summary(decideTests(fit.2vs1, p.value=0.05,lfc=1))
# group2 - group1
# Down               0
# NotSig           15353
# Up                 0

num = length(fit.2vs1$F.p.value)
degs.2vs1 <- topTable(fit.2vs1, coef="group2", confint=TRUE,number = num)
degs.2vs1$geneid<-row.names(degs.2vs1)
write.table(degs.2vs1, file=paste0('prad.tumor.vs.normal','.degs.pc.SD.txt'), sep='\t',row.names = TRUE, quote = FALSE);

##meta vs normal,3vs1##
# construct the contrast matrix corresponding to specified contrasts of a set of parameters
matrix.3vs1 <- makeContrasts(group3,levels=design)
# compute estimated coefficients and standard errors for a given set of contrasts
fit.3vs1 <- contrasts.fit(fit, matrix.3vs1)

# compute moderated t-statistics of differential expression by empirical Bayes moderation of the standard errors
fit.3vs1 <- eBayes(fit.3vs1)
options(digits=3)

#colnames(fit$coefficients)

summary(decideTests(fit.3vs1, p.value=0.05,lfc=1))
# group3 - group1
# Down     2025
# NotSig  12062
# Up      1266

num = length(fit.3vs1$F.p.value)
degs.3vs1 <- topTable(fit.3vs1, coef="group3", confint=TRUE,number = num)
degs.3vs1$geneid<-row.names(degs.3vs1)
write.table(degs.3vs1, file=paste0('prad.meta.vs.normal','.degs.pc.SD.txt'), sep='\t',row.names = TRUE, quote = FALSE);


##meta vs tumor,3vs2##
# construct the contrast matrix corresponding to specified contrasts of a set of parameters
matrix.3vs2 <- makeContrasts(group3-group2,levels=design)
# compute estimated coefficients and standard errors for a given set of contrasts
fit.3vs2 <- contrasts.fit(fit, matrix.3vs2)

# compute moderated t-statistics of differential expression by empirical Bayes moderation of the standard errors
fit.3vs2 <- eBayes(fit.3vs2)
options(digits=3)

#colnames(fit$coefficients)

summary(decideTests(fit.3vs2, p.value=0.05,lfc=1))
# group3 - group2
# Down     1920
# NotSig  12315
# Up       1118

num = length(fit.3vs2$F.p.value)
degs.3vs2 <- topTable(fit.3vs2, coef="group3 - group2", confint=TRUE,number = num)
degs.3vs2$geneid<-row.names(degs.3vs2)
write.table(degs.3vs2, file=paste0('prad.meta.vs.tumor','.degs.pc.SD.txt'), sep='\t',row.names = TRUE, quote = FALSE);


##gene name gene id##
genename<-read.delim("prad.meta.vs.normal.degs.txt",header=TRUE,stringsAsFactors = FALSE)
genename.filter<-filter(genename,gene_name!="")
genename.filter<-genename.filter[,9:10]

degs.2vs1.sub<-degs.2vs1[degs.2vs1$geneid%in%genename.filter$geneid,]
degs.2vs1.sub.join<-full_join(degs.2vs1.sub,genename.filter,by="geneid")
degs.3vs1.sub<-degs.3vs1[degs.3vs1$geneid%in%genename.filter$geneid,]
degs.3vs1.sub.join<-full_join(degs.3vs1.sub,genename.filter,by="geneid")
degs.3vs2.sub<-degs.3vs2[degs.3vs2$geneid%in%genename.filter$geneid,]
degs.3vs2.sub.join<-full_join(degs.3vs2.sub,genename.filter,by="geneid")





require(openxlsx)
list_of_datasets <- list("meta.vs.normal" = degs.3vs1.sub.join, "meta.vs.tumor" = degs.3vs2.sub.join, "tumor.vs.normal"=degs.2vs1.sub.join)
write.xlsx(list_of_datasets, file = "tcga.vs.su2c.proteincoding.SD.xlsx")




