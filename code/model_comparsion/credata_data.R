########################## Generate training data
library(stringr)
library(dplyr)
all_cg<-read.csv("/public/slst/home/ningwei/methylation/data/other_model/all_cg.csv") ### Feature set of six methods
all_cg<-all_cg[,1]

con <- file("/public/slst/home/ningwei/methylation/data/GDC-PANCAN.methylation450.tsv.gz", "r")
line=readLines(con,n=1)
j=1
d<-line
line=readLines(con,n=1)
while( length(line) != 0 ) {
    a<-str_sub(line,1,10)
    if(a%in%all_cg){
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

write.csv(f, file = "/public/slst/home/ningwei/methylation/data/cup_data/train_data.csv",row.names=F)