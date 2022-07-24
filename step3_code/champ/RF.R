library(stringr)
library(data.table, help, pos = 2, lib.loc = NULL)
library(impute)
library(tibble)
#library(FactoMineR)
#library(factoextra) 
library(dplyr, help, pos = 2, lib.loc = NULL)

#c<-read.csv("~/TCGA/Methylation450K/Three_method/feature/RF.csv")
feature = read.csv("~/TCGA/Methylation450K/finall_result/RF_30.csv")
feature = colnames(feature)[-c(1,32)]
cg12<-c("cg14003974","cg09420439","cg04260891","cg06909254","cg14020652","cg15994467","cg04249706","cg26552733","cg05139187","cg19747465",
"cg03190578","cg11090352")
# feature[,3] = abs(feature[,3])
# feature = feature[order(feature[,3],decreasing=T),]
# c<-feature[1:30,2]
c<-c(feature,cg12)
#c<-c[1:50,2]    
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

########################################################### 癌症数据
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

#write.csv(b,paste0(file="~/TCGA/Methylation450K/Three_method/data/RF.csv"),quote=F,row.names=F)
write.csv(b,paste0(file="~/TCGA/Methylation450K/finall_result/cg12_RF.csv"),quote=F,row.names=F)