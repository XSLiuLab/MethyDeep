library(data.table)
library(impute)
library(ChAMP)
library(stringr)
library(tibble)
library(dplyr)
library(FactoMineR)
library(factoextra) 

cancer_data<-list.dirs("/public/slst/home/ningwei/methylation/data/all_primary_data")
cancer_data<-cancer_data[-c(1,2)]
cancer<-str_split_fixed(cancer_data,"/",9)[,9]

for(i in cancer){
    dir.create(paste0("/public/slst/home/ningwei/methylation/data/all_primary_data/",i)) ### Methylation data for different cancers
}

args <- commandArgs(trailingOnly = TRUE)
number = as.numeric(args[1])

print(number)
for(i in number){
    for(j in 30){
        a1<-list.files(cancer_data[i])
        a2<-list.files(cancer_data[j])
        data1<-fread(paste0(cancer_data[i],"/",a1[1]),data.table=F)
        data2<-fread(paste0(cancer_data[j],"/",a2[1]),data.table=F)
        colnames(data1)<-data1[1,]
        colnames(data2)<-data2[1,]
        data1<-data1[-1,]
        data2<-data2[-1,]  
        data<-inner_join(data1,data2)
        data = column_to_rownames(data,"Composite Element REF")
        data3 <- apply(data,2,as.numeric)
        rownames(data3)<-rownames(data)  ########## Processed methylation matrix
        sample_info<-data.frame(c(rep(cancer[1],ncol(data1)-1),rep(cancer[2],ncol(data2)-1))) ###Sample information
        beta=as.matrix(data3)
        # The beta signal value matrix cannot contain NA values
        beta=impute.knn(beta) 
        beta=beta$data
        beta=beta+0.00001  #### Prevent 0 value from error during normalization
        myLoad=champ.filter(beta = beta ,pd = sample_info) #This step has automatically completed filtering
        myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=16,method="PBC") ### Normalization
        myNorm<- myLoad$beta
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
        group_list=sample_info[,1]
        group_list<-str_split_fixed(group_list,"-",2)[,2]
        myDMP <- champ.DMP(beta = myNorm,pheno=group_list)
        df_DMP <- myDMP[[1]]
        #df_DMP<-df_DMP[which(df_DMP[,1]>0.59|df_DMP[,1]<c(-0.59)),]
        df_DMP<-df_DMP[which(df_DMP[,5]<0.01),]
        for(h in 2:length(a1)){
            data3<-fread(paste0(cancer_data[i],"/",a1[h]),data.table=F)
            data4<-fread(paste0(cancer_data[j],"/",a2[h]),data.table=F)
            data3<-na.omit(data3)
            data4<-na.omit(data4)
            if(nrow(data3)==0|nrow(data4)==0){
                break
            }
            colnames(data3)<-colnames(data1)
            colnames(data4)<-colnames(data2)
            data1<-data3
            data2<-data4
            data<-inner_join(data1,data2)
            data = column_to_rownames(data,"Composite Element REF")
            data3 <- apply(data,2,as.numeric)
            rownames(data3)<-rownames(data) 
            sample_info<-data.frame(c(rep(cancer[1],ncol(data1)-1),rep(cancer[2],ncol(data2)-1))) 
            beta=as.matrix(data3)  
            beta=impute.knn(beta) 
            beta=beta$data
            beta=beta+0.00001  
            myLoad=champ.filter(beta = beta ,pd = sample_info) 
            myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=16,method="PBC")
            if(nrow(myLoad$beta)==0){
                next
            }else(
                myNorm<-myLoad$beta
            )
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
            group_list=sample_info[,1]
            group_list<-str_split_fixed(group_list,"-",2)[,2]
            myDMP <- champ.DMP(beta = myNorm,pheno=group_list)
            df_DMP1 <- myDMP[[1]]
            df_DMP1<-df_DMP1[which(df_DMP1[,1]>0.59|df_DMP1[,1]<c(-0.59)),]
            df_DMP1<-df_DMP1[which(df_DMP1[,5]<0.01),]
            df_DMP<-rbind(df_DMP,df_DMP1)
        }
    write.table(df_DMP,file=paste("/public/slst/home/ningwei/methylation/process_data/origin_differen_result/",cancer[i],"_",cancer[j],".txt",sep=""),quote=F)
 }
}