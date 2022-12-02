library(data.table)
library(stringr)
library(dplyr)


cancer_data<-list.dirs("/public/slst/home/ningwei/methylation/data/all_primary_data")
cancer_data<-cancer_data[-c(1,2)]
cancer<-str_split_fixed(cancer_data,"/",9)[,9]

args <- commandArgs(trailingOnly = TRUE)
number = as.numeric(args[1])

print(number)
for(i in number){
    for(j in 28){
        a1<-list.files(cancer_data[i])
        a2<-list.files(cancer_data[j])
        data1<-fread(paste0(cancer_data[i],"/",a1[1]),data.table=F)
        data2<-fread(paste0(cancer_data[j],"/",a2[1]),data.table=F)
        colnames(data1)<-data1[1,]
        colnames(data2)<-data2[1,]
        data1<-data1[-1,]
        data2<-data2[-1,]  
        data<-inner_join(data1,data2)
        a<-apply(data[,-1],1,function(x){wilcox.test(as.numeric(x[c(1:ncol(data1))]),as.numeric(x[(ncol(data1)+1):length(x)]))$p.value})
        differ_cg<-data[which(a<0.01),1]
        for(h in 2:length(a1)){
            data3<-fread(paste0(cancer_data[i],"/",a1[h]),data.table=F)
            data4<-fread(paste0(cancer_data[j],"/",a2[h]),data.table=F)
            data3<-na.omit(data3)
            data4<-na.omit(data4)
            if(nrow(data3)==0|nrow(data4)==0){
                next
            }
            colnames(data3)<-colnames(data1)
            colnames(data4)<-colnames(data2)
            data<-inner_join(data3,data4)
            a<-apply(data[,-1],1,function(x){wilcox.test(as.numeric(x[c(1:ncol(data1))]),as.numeric(x[(ncol(data1)+1):length(x)]))$p.value})
            differ_cg1<-data[which(a<0.05),1]
            differ_cg<-c(differ_cg,differ_cg1)
        }   
        print(j)
        write.table(differ_cg,file=paste("/public/slst/home/ningwei/methylation/process_data/wilcox_differen_result/",cancer[i],"_",cancer[j],".txt",sep=""),quote=F)
    }
}