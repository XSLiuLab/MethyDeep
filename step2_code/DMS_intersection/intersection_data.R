#library(ChAMP)
library(stringr)
library(data.table, help, pos = 2, lib.loc = NULL)
library(impute)
library(tibble)
#library(FactoMineR)
#library(factoextra) 
library(dplyr, help, pos = 2, lib.loc = NULL)

feature_number<-read.table("~/TCGA/Methylation450K/cg500.txt")
feature_type<-table(feature_number)
for(i in 0:1){
        a<-which(feature_type>i)
        c<-rownames(feature_type)[a]
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

        ########################################################### cancer data
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

        write.csv(b,paste0(file="~/TCGA/Methylation450K/data_6_15/","feature_",i,".csv"),quote=F,row.names=F)
}