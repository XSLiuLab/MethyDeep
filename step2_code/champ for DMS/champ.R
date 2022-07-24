################################################################ Loading Environment
rm(list = ls())
library(ChAMP)
library(stringr)
library(data.table, help, pos = 2, lib.loc = NULL)
library(impute)
library(tibble)
library(FactoMineR)
library(factoextra) 
library(dplyr, help, pos = 2, lib.loc = NULL)
library(EnhancedVolcano)
############################################################# Read in data
methydata<-fread("jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv.synapse_download_5096262.xena.gz",data.table=F)
print(dim(methydata))
############################################################ 数据预处理
#methydata<-d
methydata<-as.data.frame(methydata)
a = column_to_rownames(methydata,"sample")
a<-a[,!duplicated(colnames(a))]
b<-str_sub(colnames(a),14,15)
b1<-grep("01$",colnames(a))
b2<-grep("03$",colnames(a))
b3<-grep("06$",colnames(a))
b4<-grep("11$",colnames(a))
a<-a[,c(b1,b2,b3,b4)]                #########Only the original and transferred samples are retained
a<-a[,-c(8849,8850,8851)]            #########For samples with both primary and metastatic, retain the primary
########################################################################################### Read in sample information
cancer<-fread("~/TCGA/Methylation450K/TCGA_phenotype_denseDataOnlyDownload.tsv.gz",data.table=F)

normalsample<-cancer[which(cancer[,3]=="Solid Tissue Normal"),]
cancersample<-cancer[cancer[,3] %in% c("Additional Metastatic","Metastatic","Additional - New Primary","Primary Blood Derived Cancer - Peripheral Blood",
                                          "Recurrent Tumor","Primary Tumor"),]
normalsample<-normalsample[normalsample[,1] %in% colnames(a),]
cancersample<-cancersample[cancersample[,1] %in% colnames(a),]

m<-as.data.frame(which(table(cancersample[,4])<50))
cancersample<-cancersample[!(cancersample[,4] %in% rownames(m)),]

############################################################################# Extract cancer data
a<-a[,colnames(a) %in% cancersample[,1]]
a<-t(a)
a<-as.data.frame(a)
a[,ncol(a)+1]<-rownames(a)
colnames(a)[ncol(a)]<-colnames(cancersample)[1]
b<-inner_join(cancersample,a,by=colnames(cancersample)[1])
b<-b[,-c(1:3)]
a<-b[,-1]
a<-apply(a,2,function(x) {as.numeric(x)})
b[,-1]<-a
b<-na.omit(b)
colnames(b)[1]<-"CancerTypes"
########################################################################## Merge similar cancers, and remove cancer types with less than 50 samples
b[(b[,1] %in% c("colon adenocarcinoma","rectum adenocarcinoma")),1]<-"colorectal cancer"
b[(b[,1] %in% c("brain lower grade glioma","glioblastoma multiforme")),1]<-"Brain_Gliomas"

c<-as.data.frame(which(table(b[,1])<50))

b<-b[!(b[,1] %in% rownames(c)),]

a<-which(b[,1]=="esophageal carcinoma")     
b<-b[-a,]
a<-which(b[,1]=="cervical & endocervical cancer")
b<-b[-a,]

a<-which(b[,1]=="kidney chromophobe")
b<-b[-a,]
a<-which(b[,1]=="uterine carcinosarcoma")
b<-b[-a,]


#cancertype<-unique(cancersample[,4])
cancertype<-unique(b[,1])
differ_tag1<-c("wfe")
###################################################################### The difference between two cancer samples
for(i in 1:(length(cancertype)-1)){
   methysample1<-cancersample[which(cancersample[,4]==cancertype[i]),]
    for(j in (i+1):length(cancertype)){
    methysample2<-cancersample[which(cancersample[,4]==cancertype[j]),]
    patient = colnames(a)
    tp<-a[,patient %in% methysample1[,1]]
    np<-a[,patient %in% methysample2[,1]]
    a1<-cbind(tp,np)
#     methysample1[,4]<-"Tumor"
#     methysample2[,4]<-"Normal"
    methysample1[,4]<-cancertype[i]
    methysample2[,4]<-cancertype[j]
    pd<-rbind(methysample1,methysample2)
    beta=as.matrix(a1)
# Na value cannot be found in the beta signal value matrix  
beta=impute.knn(beta) 
sum(is.na(beta))

beta=beta$data
beta=beta+0.00001
myLoad=champ.filter(beta = beta,pd=pd)

myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=16,method="PBC")
num.na <- apply(myNorm,2,function(x)(sum(is.na(x))))
table(num.na)
names(num.na) = colnames(myNorm)
dt = names(num.na[num.na>0])
if(length(dt)!=0){
keep = setdiff(colnames(myNorm),dt)
myNorm = myNorm[,keep]
pd = myLoad$pd
pd = pd[pd$patient  %in% keep,]
}
########################### Draw a picture
dat <- t(myNorm)

group_list=pd$V4
table(group_list)
#> group_list
#> Normal  Tumor 
#>     29     29

dat.pca <- PCA(dat, graph = FALSE) 

p<-fviz_pca_ind(dat.pca,
             geom.ind = "point", 
             col.ind = group_list, 
             addEllipses = TRUE, 
             legend.title = "Groups")
ggsave(p,filename = paste(paste(paste("~/TCGA/Methylation450K/5_14_pca-picture/",cancertype[i],sep=""),cancertype[j],sep="_"),".png",sep=""),width = 12,height = 9)


# Correlation matrix heat map
cg=names(tail(sort(apply(myNorm,1,sd)),1000))
ac=data.frame(group=group_list)
rownames(ac)=colnames(myNorm) 
p<-pheatmap::pheatmap(cor(myNorm[cg,]),
                   annotation_col = ac,
                   show_rownames = F,
                   show_colnames = F)
ggsave(p,filename = paste(paste(paste("~/TCGA/Methylation450K/5_14_pheat-picture/",cancertype[i],sep=""),cancertype[j],sep="_"),".png",sep=""),width = 12,height = 9)
# ##############################


# group_list <- pd[,4]
# myDMP <- champ.DMP(beta = myNorm,pheno=group_list)
# head(myDMP$Tumor_to_Normal)


# df_DMP <- myDMP$Tumor_to_Normal
# df_DMP=df_DMP[df_DMP$gene!="",]

# write.table(df_DMP,file=paste(paste(paste("~/TCGA/Methylation450K/champ_differ/",cancertype[i],sep=""),cancertype[j],sep="_"),".txt",sep=""),quote=F)

# logFC_t <- 0.45
# P.Value_t <- 10^-15
# df_DMP$change <- ifelse(df_DMP$adj.P.Val < P.Value_t & abs(df_DMP$logFC) > logFC_t,
#                         ifelse(df_DMP$logFC > logFC_t ,'UP','DOWN'),'NOT') 
# differ_tag<-rownames(df_DMP)[which(df_DMP$change!="NOT")]
# differ_tag1<-c(differ_tag1,differ_tag)
# differ_tag1<-unique(differ_tag1)
# write.table(differ_tag1,file="differ_methy.txt")
# print(i)
# print(j)
 }
 }

# ################################################################### For the results of the two differences, take the top 50 as the next step
library(data.table)
a<-list.files("~/TCGA/Methylation450K/9_3_champ_两两差异/")
a<-a[-(which(a=="differ_methy.txt"))]
d<-data.frame()
for(i in 1:length(a)){
     b<-read.table(paste("~/TCGA/Methylation450K/9_3_champ_两两差异/",a[i],sep=""),fill=TRUE,header=T)
     b<-b[which(b[,2]>0.45 & b[,6]<10^-15),]
     b<-b[order(b[,6]),]
     b<-b[order(abs(b[,2]), na.last = TRUE, decreasing = T),]
     c<-b[1:500,1]
     d<-c(d,c)
     #d<-unique(d)
     print(i)
}
d<-unlist(d)
c<-table(c)
e<-c[which(c>8)]
c<-rownames(e)
write.table(d,file="~/TCGA/Methylation450K/cg500.txt",quote=F)
library(stringr)

# ########################################################################   Extract data line by line
1
library(data.table)
data<-fread("~/TCGA/Methylation450K/feature_importance_6_15/gridsearch_code/cg.txt",data.table=F)
a<-grep("^cg",data[1,])
data<-data[1,a]
for(i in 1:10){
x1 <- round(runif(50, 1, 394363) )
c<-data[x1]
c<-as.character(c)
con <- file("~/TCGA/Methylation450K/jhu-usc.edu_PANCAN_HumanMethylation450.betaValue_whitelisted.tsv.synapse_download_5096262.xena.gz", "r")
line=readLines(con,n=1)
j=0
b<-data.frame()
while( length(line) != 0 ) {
     if(j!=1){
     a<-str_sub(line,1,10)
     if(a %in% c){
     b<-rbind(b,line)
     }
     }
     line=readLines(con,n=1)
     j<-j+1
     print(j)
}
close(con)

d<-apply(b,1,function(x) {str_split_fixed(x,"\t",9665)})
d<-t(d)

methydata<-fread("~/TCGA/Methylation450K/test.txt",data.table=F)
colnames(d)<-colnames(methydata)

# b<-grep("A$",colnames(d))
# b1<-grep("A1$",colnames(d))
# b2<-grep("A2$",colnames(d))

# b<-c(b,b1,b2)
# d<-d[,c(1,b)]

methydata<-d
methydata<-as.data.frame(methydata)
a = column_to_rownames(methydata,"sample")
a<-a[,!duplicated(colnames(a))]

b<-str_sub(colnames(a),14,15)
b1<-grep("01$",colnames(a))
b2<-grep("03$",colnames(a))
b3<-grep("06$",colnames(a))
b4<-grep("11$",colnames(a))

a<-a[,c(b1,b2,b3,b4)]
a<-a[,-c(8849,8850,8851)]

cancer<-fread("~/TCGA/Methylation450K/TCGA_phenotype_denseDataOnlyDownload.tsv.gz",data.table=F)

normalsample<-cancer[which(cancer[,3]=="Solid Tissue Normal"),]
cancersample<-cancer[cancer[,3] %in% c("Additional Metastatic","Metastatic","Additional - New Primary","Primary Blood Derived Cancer - Peripheral Blood",
                                          "Recurrent Tumor","Primary Tumor"),]
normalsample<-normalsample[normalsample[,1] %in% colnames(a),]
cancersample<-cancersample[cancersample[,1] %in% colnames(a),]

m<-as.data.frame(which(table(cancersample[,4])<50))
cancersample<-cancersample[!(cancersample[,4] %in% rownames(m)),]

########################################################### Cancer data
a<-a[,colnames(a) %in% cancersample[,1]]
a<-t(a)
a<-as.data.frame(a)
a[,ncol(a)+1]<-rownames(a)
colnames(a)[ncol(a)]<-colnames(cancersample)[1]
b<-inner_join(cancersample,a,by=colnames(cancersample)[1])
sample_type<-b[,3]
b<-b[,-c(1:3)]
a<-b[,-1]
a<-apply(a,2,function(x) {as.numeric(x)})
b[,-1]<-a
b$sample_type<-sample_type
b<-na.omit(b)
colnames(b)[1]<-"CancerTypes"

b[(b[,1] %in% c("colon adenocarcinoma","rectum adenocarcinoma")),1]<-"colorectal cancer"
b[(b[,1] %in% c("brain lower grade glioma","glioblastoma multiforme")),1]<-"Brain_Gliomas"

c<-as.data.frame(which(table(b[,1])<50))

b<-b[!(b[,1] %in% rownames(c)),]

a<-which(b[,1]=="esophageal carcinoma")
if(length(a)>0){
            b<-b[-a,]
        }
a<-which(b[,1]=="cervical & endocervical cancer")
if(length(a)>0){
            b<-b[-a,]
        }

a<-which(b[,1]=="kidney chromophobe")
if(length(a)>0){
            b<-b[-a,]
        }
a<-which(b[,1]=="uterine carcinosarcoma")
if(length(a)>0){
            b<-b[-a,]
        }

#write.csv(b[,52],file="~/TCGA/Methylation450K/9_17_sample.csv",quote=F,row.names=F)

#write.csv(b,file="~/TCGA/Methylation450K/cg_7093_6_14cancertype.csv",quote=F,row.names=F)
write.csv(b,file=paste0("~/TCGA/Methylation450K/random_feature/random_",i,".csv"),quote=F,row.names=F)
}

c<-b[b[,1]%in%c("breast invasive carcinoma","colorectal cancer","esophageal carcinoma","lung adenocarcinoma","lung squamous cell carcinoma","stomach adenocarcinoma"),]
c<-c[,-7]
c<-c[-which(c[,1]=="esophageal carcinoma"),]
c[c[,1]%in%c("lung adenocarcinoma","lung squamous cell carcinoma"),1]<-"lung cancer"
write.csv(c,file="~/TCGA/Methylation450K/test/five_cancer.csv",row.names=F)

########################################
