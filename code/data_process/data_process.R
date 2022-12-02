#################### Sample information processing
library(data.table)
library(stringr)
sample_info<-fread("/public/slst/home/ningwei/methylation/data/GDC-PANCAN.basic_phenotype.tsv.gz",data.table=F)
a<-str_split_fixed(sample_info[,1],"-",4)
sample_info<-sample_info[grep("A",a[,4]),]  ### Only sample A is retained

primary_sample<-sample_info[sample_info[,3]%in%c(1,3,5,9),c(1,5)] ### 14174 sample
metastasis_sample<-sample_info[sample_info[,3]%in%c(6,7),c(1,5)] ### 411 sample
Recurrent_sample<-sample_info[sample_info[,3]%in%c(2,4,40),c(1,5)] ### 310 sample


write.table(metastasis_sample,file="/public/slst/home/ningwei/methylation/data/process_data/metastasis_sample.txt")
write.table(Recurrent_sample,file="/public/slst/home/ningwei/methylation/data/process_data/Recurrent_sample.txt")

############ Referring to other literatures, we only reserved cancer types with sample number>50 for the training set
a<-as.data.frame(table(primary_sample[,2]))
b<-a[a[,2]<50,1] ###  TARGET-CCSK TARGET-RT   TCGA-DLBC
primary_sample<-primary_sample[!primary_sample[,2]%in%b,] ### 14080 sample

################################### Methylation data processing
con <- file("/public/slst/home/ningwei/methylation/data/GDC-PANCAN.methylation450.tsv.gz", "r")
line=readLines(con,n=1)
a<-str_split_fixed(line,"\t",9737)
a<-as.character(a)
write.table(a,file="/public/slst/home/ningwei/methylation/data/process_data/methydata_sample.txt")

b<-intersect(primary_sample[,1],a[-1])

primary_sample<-primary_sample[primary_sample[,1]%in%b,]
c<-as.data.frame(table(primary_sample[,2]))
c<-c[c[,2]<50,1] ###  TCGA-CHOL TCGA-OV  
primary_sample<-primary_sample[!primary_sample[,2]%in%c,] ### 8296 sample
write.csv(primary_sample,file="/public/slst/home/ningwei/methylation/data/process_data/primary_sample.csv")

####################################
args <- commandArgs(trailingOnly = TRUE)
file = args[1]
print(args)
print(file)
library(stringr)
primary_sample<-read.csv("/public/slst/home/ningwei/methylation/data/process_data/primary_sample.csv")
primary_sample<-primary_sample[,-1]
cancer<-unique(primary_sample[,2])
m<-cancer[30]
cancer_data<-data.frame()
cancer_sample<-primary_sample[primary_sample[,2]==m,1]
file <- list.files("/public/slst/home/ningwei/methylation/data/split_data")
a<-read.table("/public/slst/home/ningwei/methylation/data/process_data/methydata_sample.txt")
a<-a[,1]
c<-which(a%in%cancer_sample)
c<-c(1,c)
con <- file(paste0("/public/slst/home/ningwei/methylation/data/new_split/",file), "r")
line=readLines(con,n=1)
j=1
while( length(line) != 0 ) {
    d<-str_split_fixed(line,"\t",9737)
    d<-d[c]
    if(length(which(d==""))==0){
        cancer_data<-rbind(cancer_data,d)
        line=readLines(con,n=1)
    }else{
        line=readLines(con,n=1)
    }
    print(j)
    j = j+1
    }
close(con)
write.csv(cancer_data,file=paste0("/public/slst/home/ningwei/methylation/data/all_primary_data/",m,"/",m,"_",file,".csv"),row.names=F)

########################## After completing the difference, generate training data
library(stringr)
all_cg<-read.csv("/public/slst/home/ningwei/methylation/process_data/all_cg.csv")
con <- file("/public/slst/home/ningwei/methylation/data/GDC-PANCAN.methylation450.tsv.gz", "r")
line=readLines(con,n=1)
j=1
d<-line
line=readLines(con,n=1)
while( length(line) != 0 ) {
    a<-str_sub(line,1,10)
    if(a%in%all_cg[,1]){
        d<-rbind(d,line)
    }
    line=readLines(con,n=1)
    print(j)
    j = j+1
    }
close(con)

e<-apply(d,1,function(x){str_split_fixed(x,"\t",9737)})
e<-as.data.frame(e)
colnames(e)<-e[1,]
e<-e[-1,]

primary_sample<-read.csv("/public/slst/home/ningwei/methylation/data/process_data/primary_sample.csv")
primary_sample<-primary_sample[,-1]
colnames(primary_sample)[1]<-colnames(e)[1]<-"sample"

f<-inner_join(primary_sample,e,by="sample")
f1<-apply(f[,-c(1:2)],2,as.numeric)
f1[is.na(f1)]<-0
f[,-c(1:2)]<-f1
write.csv(f[,-1],file="/public/slst/home/ningwei/methylation/process_data/train_data/data_11_19.csv",row.names=F)
