######################################
#################################额外验证数据
library(data.table)
#RF<-read.csv("F:/result/result_11_15/process_data/feature_importance_11_23/RF_feature_11_23.csv")
#EXtra<-read.csv("F:/result/result_11_15/process_data/feature_importance_11_23/Extratree_feature_11_23.csv")
#XGBOOST<-read.csv("F:/result/result_11_15/process_data/feature_importance_11_23/Xgboost_feature_11_23.csv")
all_cg<-read.csv(("F:/result/result_11_15/other_model/all_cg.csv"))

#methy_700<-t(XGBOOST[1:50,])

methy_700<-t(all_cg)
colnames(methy_700)<-methy_700[1,]
###################################第一个乳腺癌 34个样本
breast <-fread("F:/甲基化数据/GSE178235_series_matrix.txt",data.table=F)
a<-grep("*microdissected, cancer cell*",colnames(breast))
breast<-breast[,c(1,a)]
a<-sapply(colnames(methy_700),function(x,y){which(y==x)},breast[,1])
a<-unlist(a)
breast<-breast[a,]
colnames(breast)<-1:ncol(breast)
write.csv(breast,file = "F:/result/result_11_15/validata/primary/breast1.csv",quote = F,row.names = F)

######################################### 第二个乳腺癌 41个样本
breast3<-fread("F:\\甲基化数据/GSE141338_series_matrix.txt",data.table = F)
a<-sapply(colnames(methy_700),function(x,y){which(y==x)},breast3[,1])
a<-unlist(a)
breast3<-breast3[a,]

a<-grep("*normal*",colnames(breast3))
breast3<-breast3[,-a]
breast3<-breast3[,-c(43)]
write.csv(breast3,file = "F:/result/result_11_15/validata/primary/breast2.csv",row.names = F,quote = F)

########################################## 第三个乳腺癌 111个样本
breast4<-fread("F:/甲基化数据/GSE141441_series_matrix.txt",data.table = F)
a<-sapply(colnames(methy_700),function(x,y){which(y==x)},breast4[,1])
a<-unlist(a)
breast4<-breast4[a,]

a<-which(colnames(breast4)=="study_arm: Chemo")
breast4<-breast4[,-a]

write.csv(breast4,file = "F:/result/result_11_15/validata/primary/breast3.csv",row.names = F,quote = F)

########################################## 第四个乳腺癌 47个样本
breast5<-fread("F:/甲基化数据/GSE106360_series_matrix.txt",data.table = F)
a<-sapply(colnames(methy_700),function(x,y){which(y==x)},breast5[,1])
a<-unlist(a)
breast5<-breast5[a,]
write.csv(breast5,file = "F:/result/result_11_15/validata/primary/breast4.csv",row.names = F,quote = F)

########################  第5个乳腺癌
breast6<-fread("F:/甲基化数据/GSE78751_series_matrix.txt",data.table = F)
a<-sapply(colnames(methy_700),function(x,y){which(y==x)},breast6[,1])
a<-unlist(a)
breast6<-breast6[a,]

a<-grep("^genomic DNA from matched normal adjacent tissue Patient",colnames(breast6))
breast6<-breast6[,-a]

a<-grep("^genomic DNA from lymph node metastasis",colnames(breast6)) ### 转移样本
breast7<-breast6[,c(1,a)]
write.csv(breast7,file = "F:/result/result_11_15/validata/metastatic/breast1.csv",row.names = F,quote = F)

a<-grep("^genomic DNA from lymph node metastasis",colnames(breast6)) ### 转移样本
breast7<-breast6[,-a]
write.csv(breast7,file = "F:/result/result_11_15/validata/primary/breast5.csv",row.names = F,quote = F)


##################  胃癌 13个样本

stomach<-fread("F:/甲基化数据/GSE164988_series_matrix.txt",data.table = F)
a<-sapply(colnames(methy_700),function(x,y){which(y==x)},stomach[,1])
a<-unlist(a)
stomach<-stomach[a,]

a<-grep("normal stomach tissue sample",colnames(stomach))
stomach<-stomach[,-a]

write.csv(stomach,file = "F:/result/result_11_15/validata/primary/stomach1.csv",row.names = F,quote = F)

##################### 胃癌 （弥漫性_癌组织和肠癌组织） 25个样本

stomach2<-fread("F:/甲基化数据/GSE127857_series_matrix.txt",data.table = F)
a<-sapply(colnames(methy_700),function(x,y){which(y==x)},stomach2[,1])
a<-unlist(a)
stomach2<-stomach2[a,]

a<-grep("^Diffuse_cancer tissue_donor",colnames(stomach2))
a2<-grep("^Intestinal_cancer tissue_donor",colnames(stomach2))
stomach4<-stomach2[,c(1,a,a2)]
colnames(stomach4)<-1:ncol(stomach4)
write.csv(stomach4,file = "F:/result/result_11_15/validata/primary/stomach2.csv",row.names = F,quote = F)
################################################## 胰腺癌 82个样本

pancreatic<-fread("F:/甲基化数据/GSE155353_series_matrix.txt.gz",data.table=F,fill=T)
a<-sapply(colnames(methy_700),function(x,y){which(y==x)},pancreatic[,1])
a<-unlist(a)
pancreatic<-pancreatic[c(35,a),]
#pancreatic<-as.data.frame(pancreatic)

a<-grep("T$",pancreatic[1,])
pancreatic<-pancreatic[,c(1,a)]
pancreatic<-pancreatic[-1,]
write.csv(pancreatic,file = "F:/result/result_11_15/validata/primary/pancreatic.csv",row.names = F,quote = F)

############################################# 第二个胰腺癌  6个样本
pancreatic2<-fread("F:/甲基化数据/GSE149250_series_matrix.txt.gz",data.table = F,fill=T)
a<-sapply(colnames(methy_700),function(x,y){which(y==x)},pancreatic2[,1])
a<-unlist(a)
pancreatic2<-pancreatic2[c(59,a),]

a<-grep("^PDAC",pancreatic2[1,])
pancreatic2<-pancreatic2[,c(1,a)]
pancreatic2<-pancreatic2[-1,]
write.csv(pancreatic2,file = "F:/result/result_11_15/validata/primary/pancreatic2.csv",row.names = F,quote = F)

############################################# 第三个胰腺癌  6个样本
pancreatic3<-fread("F:/甲基化数据/GSE49149_series_matrix.txt",data.table = F)
a<-grep("^Tumor",colnames(pancreatic3))
pancreatic3<-pancreatic3[,c(1,a)]
a<-sapply(colnames(methy_700),function(x,y){which(y==x)},pancreatic3[,1])
a<-unlist(a)
pancreatic3<-pancreatic3[a,]

write.csv(pancreatic3,file = "F:/result/result_11_15/validata/primary/pancreatic3.csv",row.names = F,quote = F)

################################################ 肺
library(stringr)
a<-fread("F:/甲基化数据/GSE158075_series_matrix.txt.gz",data.table=F,fill=T)

a1<-which(a[37,]=="diagnosis (final clinical diagnosis for the patient): squamous cell carcinoma")
a1<-a[,c(1,a1)]
a1<-a1[35,-1]
a1<-apply(a1, 2, function(x){str_split_fixed(x,": ",2)[,2]})


a2<-which(a[37,]=="diagnosis (final clinical diagnosis for the patient): adenocarcinoma")
a2<-a[,c(1,a2)]
a2<-a2[35,-1]
a2<-apply(a2, 2, function(x){str_split_fixed(x,": ",2)[,2]})

####################################### 肺鳞状细胞癌 73个样本
lung<-fread("F:/甲基化数据/GSE158075_GA_illumina_methylation_MatrixProcessed_PB3.txt.gz",data.table = F)

a<-sapply(colnames(methy_700),function(x,y){which(y==x)},lung[,1])
a<-unlist(a)
lung<-lung[a,]

lung<-lung[,-grep("Pval$",colnames(lung))]
b<-c()
for(i in 1:length(a1)){
  a<-grep(a1[i],colnames(lung))
  b<-c(b,a)
}
b<-unique(b)
lung_scc<-lung[,c(1,b)]
write.csv(lung_scc,file = "F:/result/result_11_15/validata/primary/lung_scc.csv",row.names = F,quote = F)

##################################  肺腺癌 14个样本
b<-c()
for(i in 1:length(a2)){
  a<-grep(a2[i],colnames(lung))
  b<-c(b,a)
}
b<-unique(b)

lung_aden<-lung[,c(1,b)]

write.csv(lung_aden,file = "F:/result/result_11_15/validata/primary/lung_aden.csv",row.names = F,quote = F)

####################################  第二个肺鳞状细胞癌 20个样本
lung_scc2<-fread("F:/甲基化数据/GSE121849_series_matrix.txt",data.table=F)
a<-sapply(colnames(methy_700),function(x,y){which(y==x)},lung_scc2[,1])
a<-unlist(a)
lung_scc2 <- lung_scc2[a,]

a1<-which(colnames(lung_scc2)=="lung squamous cell carcinoma with idiopathic pulmonary fibrosis")
a2<-which(colnames(lung_scc2)=="lung squamous cell carcinoma without idiopathic pulmonary fibrosis")
lung_scc2<-lung_scc2[,c(1,a1,a2)]

write.csv(lung_scc2,file = "F:/result/result_11_15/validata/primary/lung_scc2.csv",row.names = F,quote = F)
########################################  肝细胞癌 47个样本
liver <- fread("F:/甲基化数据/GSE136319_series_matrix.txt.gz",data.table=F,fill=T)
a<-sapply(colnames(methy_700),function(x,y){which(y==x)},liver[,1])
a<-unlist(a)
liver<-liver[c(37,a),]
liver<-liver[,c(1,grep("^HCC",liver[1,]))]
liver<-liver[-1,]
write.csv(liver,file = "F:/result/result_11_15/validata/primary/liver.csv",row.names = F,quote = F)

####################################  第二个肝细胞癌  15个样本
liver2 <- fread("F:/甲基化数据/GSE99036_series_matrix.txt",data.table=F)
a<-sapply(colnames(methy_700),function(x,y){which(y==x)},liver2[,1])
a<-unlist(a)
liver2<-liver2[a,]

a1<-which(colnames(liver2)=="Hepatocellular carcinoma")
liver2<-liver2[,c(1,a1)]
write.csv(liver2,file = "F:/result/result_11_15/validata/primary/liver2.csv",row.names = F,quote = F)
####################################  第三个肝细胞癌  37个样本
liver3 <- fread("F:/甲基化数据/GSE113019_series_matrix.txt",data.table=F)
a<-sapply(colnames(methy_700),function(x,y){which(y==x)},liver3[,1])
a<-unlist(a)
liver3<-liver3[a,]

a1<-which(colnames(liver3)=="tissue: Liver cancer, non-tumorous tissue")
liver3<-liver3[,-a1]

a<-grep("tissue: Liver cancer, recurrent tumorous tissue",colnames(liver3))
liver4<-liver3[,c(1,a)]
colnames(liver4)<-1:ncol(liver4)
write.csv(liver4,file = "F:/result/result_11_15/validata/recurren/liver3_recurrent_18.csv",row.names = F,quote = F)

a<-grep("tissue: Liver cancer, recurrent tumorous tissue",colnames(liver3))
liver4<-liver3[,-a]
colnames(liver4)<-1:ncol(liver4)
write.csv(liver4,file = "F:/result/result_11_15/validata/primary/liver3.csv",row.names = F,quote = F)
####################################  第四个肝细胞癌  37个样本
liver4 <- fread("F:/甲基化数据/GSE113017_series_matrix.txt",data.table=F)
a<-sapply(colnames(methy_700),function(x,y){which(y==x)},liver4[,1])
a<-unlist(a)
liver4<-liver4[a,]

a1<-grep("N$",colnames(liver4))
liver4<-liver4[,-a1]
a1<-grep("^L022",colnames(liver4))
liver4<-liver4[,-a1]
write.csv(liver4,file = "F:/result/result_11_15/validata/primary/liver4.csv",row.names = F,quote = F)
####################################  第五个肝细胞癌  37个样本
liver5 <- fread("F:/甲基化数据/GSE89852_series_matrix.txt",data.table=F)
a<-sapply(colnames(methy_700),function(x,y){which(y==x)},liver5[,1])
a<-unlist(a)
liver5<-liver5[a,]

a1<-grep("^Cancerous tissue",colnames(liver5))
liver5<-liver5[,c(1,a1)]

write.csv(liver5,file = "F:/result/result_11_15/validata/primary/liver5.csv",row.names = F,quote = F)


####################################  皮肤黑色素瘤 38个样本
melanoma<-fread("F:/甲基化数据/GSE140169_series_matrix.txt.gz",data.table = F,fill=T)
a<-sapply(colnames(methy_700),function(x,y){which(y==x)},melanoma[,1])
a<-unlist(a)
melanoma<-melanoma[c(45,a),]

a<-which(melanoma[1,]=="tissue: skin")
melanoma<-melanoma[,c(1,a)]
melanoma<-melanoma[-1,]
write.csv(melanoma,file = "F:/result/result_11_15/validata/primary/melanoma.csv",row.names = F,quote = F)
################################## 第二个皮肤黑色素 89个   
melanoma2<-fread("F:/甲基化数据/GSE120878_series_matrix.txt",data.table = F)
a<-sapply(colnames(methy_700),function(x,y){which(y==x)},melanoma2[,1])
a<-unlist(a)
melanoma2<-melanoma2[a,]

a<-which(colnames(melanoma2)=="nevus sample")
melanoma2<-melanoma2[,-a]
write.csv(melanoma2,file = "F:/result/result_11_15/validata/primary/melanoma2.csv",row.names = F,quote = F)
# ############################# 肉瘤 44个样本 不是软组织肉瘤
# sarcoma<-fread("f:/甲基化数据/GSE125645_series_matrix.txt",data.table = F)
# a<-sapply(colnames(methy_700),function(x,y){which(y==x)},sarcoma[,1])
# a<-unlist(a)
# sarcoma<-sarcoma[a,]
# 
# a1<-grep("^Osteosarcoma metastatic tumor tissue",colnames(sarcoma))
# a2<-grep("^Osteosarcoma primary tumor tissue",colnames(sarcoma))
# a3<-grep("^Ewing sarcoma primary tumor tissue",colnames(sarcoma))
# 
# sarcoma1<-sarcoma[,c(1,a2,a3)]
# write.csv(sarcoma1,file = "F:/result/result_11_15/validata/primary/sarcoma.csv",row.names = F,quote = F)
# 
# sarcoma1<-sarcoma[,c(1,a1)]
# write.csv(sarcoma1,file = "F:/result/result_11_15/validata/metastatic/sarcoma.csv",row.names = F,quote = F)

##################################### 结直肠癌 40个样本

Colorectal<-fread("F:/甲基化数据/GSE139404_series_matrix.txt",data.table=F)

a<-sapply(colnames(methy_700),function(x,y){which(y==x)},Colorectal[,1])
a<-unlist(a)
Colorectal<-Colorectal[a,]

a<-which(colnames(Colorectal)=="disease: normal")
Colorectal<-Colorectal[,-a]
write.csv(Colorectal,file = "F:/result/result_11_15/validata/primary/colorectal.csv",row.names = F,quote = F)

######################################## 第二个结直肠癌  69个样本
Colorectal2<-fread("F:/甲基化数据/GSE129364_series_matrix.txt",data.table=F)

a<-sapply(colnames(methy_700),function(x,y){which(y==x)},Colorectal2[,1])
a<-unlist(a)
Colorectal2<-Colorectal2[a,]

a<-which(colnames(Colorectal2)=="normal colon mucosa")
Colorectal2<-Colorectal2[,-a]

a<-which(colnames(Colorectal2)=="recurrent colorectal adenoma")
Colorectal3<-Colorectal2[,c(1,a)]
write.csv(Colorectal3,file = "F:/result/result_11_15/validata/recurren/colon_recurrent_1.csv",row.names = F,quote = F)

Colorectal3<-Colorectal2[,-a]
write.csv(Colorectal3,file = "F:/result/result_11_15/validata/primary/colorectal2.csv",row.names = F,quote = F)
######################################## 第三个结直肠癌  10个样本
Colorectal3<-fread("F:/甲基化数据/GSE106556_series_matrix.txt",data.table=F)

a<-sapply(colnames(methy_700),function(x,y){which(y==x)},Colorectal3[,1])
a<-unlist(a)
Colorectal3<-Colorectal3[a,]

a<-grep("LSTP$",colnames(Colorectal3))
Colorectal3<-Colorectal3[,-a]
write.csv(Colorectal3,file = "F:/result/result_11_15/validata/primary/colorectal3.csv",row.names = F,quote = F)

#######################################  胶质瘤 122个样本

lGG<-fread("F:/甲基化数据/GSE129477_series_matrix.txt",data.table=F)
a<-sapply(colnames(methy_700),function(x,y){which(y==x)},lGG[,1])
a<-unlist(a)
lGG<-lGG[a,]

write.csv(lGG,file = "F:/result/result_11_15/validata/primary/LGG.csv",row.names = F,quote = F)

#################################  多形性胶质母细胞瘤 35个

gbm<-fread("F:/甲基化数据/GSE128654_series_matrix.txt",data.table=F)

a<-sapply(colnames(methy_700),function(x,y){which(y==x)},gbm[,1])
a<-unlist(a)
gbm<-gbm[a,]

a<-grep("tumor$",colnames(gbm))
gbm<-gbm[,c(1,a)]

write.csv(gbm,file = "F:/result/result_11_15/validata/primary/gbm.csv",row.names = F,quote = F)

################################ 肾透明细胞癌 144个
ccRcc<-fread("F:/甲基化数据/GSE113501_series_matrix.txt",data.table=F)

a<-sapply(colnames(methy_700),function(x,y){which(y==x)},ccRcc[,1])
a<-unlist(a)
ccRcc<-ccRcc[a,]
colnames(ccRcc)<-1:ncol(ccRcc)
write.csv(ccRcc,file = "F:/result/result_11_15/validata/primary/ccRcc.csv",row.names = F,quote = F)

############################### 肾上腺皮质癌
Adrenal_gland<-fread("F:/甲基化数据/GSE77871_series_matrix.txt",data.table=F)

a<-sapply(colnames(methy_700),function(x,y){which(y==x)},Adrenal_gland[,1])
a<-unlist(a)
Adrenal_gland<-Adrenal_gland[a,]

a<-grep("^normal",colnames(Adrenal_gland))
Adrenal_gland<-Adrenal_gland[,-a]

write.csv(Adrenal_gland,file = "F:/result/result_11_15/validata/primary/Adrenal_gland.csv",row.names = F,quote = F)


############################## 头颈癌
HNC<-fread("F:/甲基化数据/GSE38268_series_matrix.txt",data.table=F)
a<-sapply(colnames(methy_700),function(x,y){which(y==x)},HNC[,1])
a<-unlist(a)
HNC<-HNC[a,]

write.csv(HNC,file = "F:/result/result_11_15/validata/primary/HNC.csv",row.names = F,quote = F)

############################## 前列腺癌
prostate<-fread("F:/甲基化数据/GSE112047_series_matrix.txt",data.table=F)

a<-sapply(colnames(methy_700),function(x,y){which(y==x)},prostate[,1])
a<-unlist(a)
prostate<-prostate[a,]

a<-grep("normal$",colnames(prostate))
prostate<-prostate[,-a]

write.csv(prostate,file = "F:/result/result_11_15/validata/primary/prostate.csv",row.names = F,quote = F)
############################## 第二个前列腺癌
prostate2<-fread("F:/甲基化数据/GSE38240_series_matrix.txt",data.table=F)

a<-sapply(colnames(methy_700),function(x,y){which(y==x)},prostate2[,1])
a<-unlist(a)
prostate2<-prostate2[a,]

a<-grep("Normal",colnames(prostate2))
prostate2<-prostate2[,-a]
#write.csv(prostate2,file = "F:/result/新转移数据/data/prostate2_8.csv",row.names = F,quote = F)
write.csv(prostate2,file = "F:/result/result_11_15/validata/metastatic/prostate2.csv",row.names = F,quote = F)
############################## 第三个前列腺癌
prostate3<-fread("F:/甲基化数据/GSE73549_series_matrix.txt",data.table=F)

a<-sapply(colnames(methy_700),function(x,y){which(y==x)},prostate3[,1])
a<-unlist(a)
prostate3<-prostate3[a,]

a<-grep("^disease state: Tumor",colnames(prostate3))
a1<-grep("^disease state: Metastasis",colnames(prostate3))


prostate4<-prostate3[,c(1,a1)]
write.csv(prostate4,file = "F:/result/result_11_15/validata/metastatic/prostate3.csv",row.names = F,quote = F)

prostate4<-prostate3[,c(1,a)]
write.csv(prostate4,file = "F:/result/result_11_15/validata/primary/prostate3.csv",row.names = F,quote = F)

################################## 甲状腺癌

thyroid<-fread("F:/甲基化数据/GSE86961_series_matrix.txt",data.table=F)

a<-sapply(colnames(methy_700),function(x,y){which(y==x)},thyroid[,1])
a<-unlist(a)
thyroid<-thyroid[a,]

a<-grep("normal$",colnames(thyroid))
thyroid<-thyroid[,-a]
write.csv(thyroid,file = "F:/result/result_11_15/validata/primary/thyroid.csv",row.names = F,quote = F)
################################ 综合数据
three<-fread("F:/甲基化数据/GSE52955_series_matrix.txt",data.table=F)

a<-sapply(colnames(methy_700),function(x,y){which(y==x)},three[,1])
a<-unlist(a)
three<-three[a,]

a<-which(colnames(three)=="Renal Tumor")
renal<-three[,c(1,a)]
write.csv(renal,file = "F:/result/result_11_15/validata/primary/renal.csv",row.names = F,quote = F)

a<-which(colnames(three)=="Bladder Tumor")
Bladder<-three[,c(1,a)]
write.csv(Bladder,file = "F:/result/result_11_15/validata/primary/Bladder.csv",row.names = F,quote = F)

a<-which(colnames(three)=="Prostate Tumor")
Prostate<-three[,c(1,a)]
write.csv(Prostate,file = "F:/result/result_11_15/validata/primary/Prostate.csv",row.names = F,quote = F)

######################### 处理甲基化数据
library(stringr)
library(data.table)
library(dplyr)
############################

#RF<-read.csv("F:/result/result_11_15/process_data/feature_importance_11_23/RF_feature_11_23.csv")
#EXtra<-read.csv("F:/result/result_11_15/process_data/feature_importance_11_23/Extratree_feature_11_23.csv")
all_cg<-read.csv(("F:/result/result_11_15/other_model/all_cg.csv"))

#methy_700<-t(XGBOOST[1:50,])

# methy_700<-t(all_cg)
# colnames(methy_700)<-methy_700[1,]
cg50<-all_cg
##################################
data<-fread("D:/job-two/new_data/txt/breast_34_cancer.txt")
data<-as.data.frame(data)
a<-data[,1] %in% cg50[,1]
data<-data[a,]
write.csv(data,file = "F:/result/result_11_15/validata/primary/breast_34_cancer.csv",row.names = F)


methylation<-fread("D:/job-two/new_data/GSE148748_series_matrix.txt.gz",data.table = F,fill=T)
methylation<-methylation[c(58:nrow(methylation)),]
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]
write.csv(methylation,file = "F:/result/result_11_15/validata/primary/breast_82_cancer.csv",row.names = F)

methylation<-fread("D:/job-two/new_data/GSE184159_beta_detp.csv.gz",data.table = F,fill=T)
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]
a<-grep("det",colnames(methylation))
methylation<-methylation[,-a]
write.csv(methylation,file = "F:/result/result_11_15/validata/primary/breast_63_cancer.csv",row.names = F)
##########################  GBM
data<-fread("D:/job-two/new_data/txt/GBM_380_cancer.txt")
data<-as.data.frame(data)
a<-data[,1] %in% cg50[,1]
data<-data[a,]
write.csv(data,file = "F:/result/result_11_15/validata/primary/GBM_380_cancer.csv",row.names = F)
data<-fread("D:/job-two/new_data/txt/GBM_60_cancer.txt")

data<-fread("D:/job-two/new_data/txt/GBM_60_cancer.txt")
data<-as.data.frame(data)
a<-data[,1] %in% cg50[,1]
data<-data[a,]
write.csv(data,file = "F:/result/result_11_15/validata/primary/GBM_60_cancer.csv",row.names = F)

methylation<-fread("D:/job-two/new_data/GSE188547_series_matrix.txt.gz",data.table = F,blank.lines.skip=T,fill=T)
methylation<-methylation[c(59:865918),]
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]
write.csv(methylation,file = "F:/result/result_11_15/validata/primary/GBM_112_cancer.csv",row.names = F)

methylation<-fread("D:/job-two/new_data/GSE122920-GPL13534_series_matrix.txt.gz",data.table = F,blank.lines.skip=T,fill=T)
test<-methylation[1:100,]
methylation<-methylation[c(76:nrow(methylation)),]
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]
write.csv(methylation,file = "F:/result/result_11_15/validata/primary/GBM_16_cancer.csv",row.names = F)

methylation<-fread("D:/job-two/new_data/GSE122920-GPL21145_series_matrix.txt.gz",data.table = F,blank.lines.skip=T,fill=T)
methylation<-methylation[c(52,76:nrow(methylation)),]
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]
write.csv(methylation,file = "F:/result/result_11_15/validata/primary/GBM_88_cancer.csv",row.names = F)

methylation<-fread("D:/job-two/new_data/GSE119774_GBM_GSC_i850K_beta_formatted.txt.gz",data.table = F,blank.lines.skip=T,fill=T)
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]
a<-grep("NSC",colnames(methylation))
methylation<-methylation[,-a]
a<-grep("GSC",colnames(methylation))
methylation<-methylation[,-a]
write.csv(methylation,file = "F:/result/result_11_15/validata/primary/GBM_41_cancer.csv",row.names = F)

methylation<-fread("D:/job-two/new_data/GSE147391_series_matrix.txt.gz",data.table = F,blank.lines.skip=T,fill=T)
methylation<-methylation[c(60:nrow(methylation)),]
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]
write.csv(methylation,file = "F:/result/result_11_15/validata/primary/GBM_16_cancer.csv",row.names = F)


methylation<-fread("D:/job-two/new_data/GSE143842-GPL13534_series_matrix.txt.gz",data.table = F,blank.lines.skip=T,fill=T)
test<-methylation[1:100,]
methylation<-methylation[c(59:nrow(methylation)),]
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]
write.csv(methylation,file = "F:/result/result_11_15/validata/primary/GBM_14_cancer.csv",row.names = F)


methylation<-fread("D:/job-two/new_data/GSE143842-GPL21145_series_matrix.txt.gz",data.table = F,blank.lines.skip=T,fill=T)
methylation<-methylation[c(59:nrow(methylation)),]
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]
write.csv(methylation,file = "F:/result/result_11_15/validata/primary/GBM_4_cancer.csv",row.names = F)
################### Acute myeloid leukemia

methylation<-fread("D:/job-two/new_data/GSE153347_series_matrix.txt.gz",data.table = F,blank.lines.skip=T,fill=T)
test<-methylation[1:100,]
methylation<-methylation[c(64:nrow(methylation)),]
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]
write.csv(methylation,file = "F:/result/result_11_15/validata/primary/AML_105_cancer.csv",row.names = F)
#write.csv(methylation,file = "F:/result/12cg/validata/AML_105_cancer.csv",row.names = F)


methylation<-fread("D:/job-two/new_data/GSE133986_series_matrix.txt.gz",data.table = F,blank.lines.skip=T,fill=T)
#test<-methylation[1:100,]
methylation<-methylation[c(81:nrow(methylation)),]
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]
write.csv(methylation,file = "F:/result/result_11_15/validata/primary/AML_64_cancer.csv",row.names = F)
#write.csv(methylation,file = "F:/result/12cg/validata/AML_64_cancer.csv",row.names = F)


methylation<-fread("D:/job-two/new_data/GSE113545_series_matrix.txt.gz",data.table = F,blank.lines.skip=T,fill=T)
test<-methylation[1:100,]
methylation<-methylation[c(61:nrow(methylation)),]
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]
write.csv(methylation,file = "F:/result/result_11_15/validata/primary/AML_31_cancer.csv",row.names = F)
#write.csv(methylation,file = "F:/result/12cg/validata/AML_31_cancer.csv",row.names = F)

########################## Lung adenocarcinoma

methylation<-fread("D:/job-two/new_data/GSE180060_processed_data.txt.gz",data.table = F,fill=T)
# test<-methylation[1:100,]
# test<-str_split_fixed(test[,1],"\t",200)
# test<-test[,1:89]
methylation<-str_split_fixed(methylation[,1],"\t",90)
methylation<-methylation[,1:89]
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]
write.csv(methylation,file = "F:/result/result_11_15/validata/primary/LUNG_ADEN__88_cancer.csv",row.names = F)
#write.csv(methylation,file = "F:/result/12cg/validata/LUNG_ADEN__88_cancer.csv",row.names = F)


methylation<-fread("D:/job-two/new_data/GSE114989_Matrixprocessed.txt.gz",data.table = F,fill=T)
#test<-methylation[1:100,]
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]
a<-grep("Pval",colnames(methylation))
methylation<-methylation[,-a]
normal<-c("SAMPLE 1","SAMPLE 12","SAMPLE 17","SAMPLE 28","SAMPLE 32","SAMPLE 35","SAMPLE 39")
methylation<-methylation[,!colnames(methylation)%in%normal]

metastasis<-c("SAMPLE 2" , "SAMPLE 3","SAMPLE 7","SAMPLE 21","SAMPLE 22","SAMPLE 29")
methylation1<-methylation[,c(1,which(colnames(methylation)%in%metastasis))]
write.csv(methylation1,file = "F:/result/result_11_15/validata/metastatic/LUNG_ADEN__6_cancer.csv",row.names = F)

methylation2<-methylation[,-which(colnames(methylation)%in%metastasis)]
write.csv(methylation,file = "F:/result/result_11_15/validata/primary/LUNG_ADEN__27_cancer.csv",row.names = F)
#write.csv(methylation,file = "F:/result/12cg/validata/LUNG_ADEN__40_cancer.csv",row.names = F)

######################### 不确定是不是腺癌
# methylation<-fread("D:/job-two/new_data/GSE115246_series_matrix.txt.gz",data.table = F,fill=T)
# test<-methylation[1:100,]
# methylation<-methylation[c(65:nrow(methylation)),]
# a<-methylation[,1] %in% cg50[,1]
# methylation<-methylation[a,]
# write.csv(methylation,file = "F:/result/result_11_15/validata/primary/LUNG_81_cancer.csv",row.names = F)
#write.csv(methylation,file = "F:/result/12cg/validata/LUNG_ADEN__81_cancer.csv",row.names = F)
###################################### rectal cancer
methylation<-fread("D:/job-two/new_data/GSE132668_series_matrix.txt.gz",data.table = F,fill=T)
#test<-methylation[1:100,]
methylation<-methylation[c(66:nrow(methylation)),1:33]
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]
write.csv(methylation,file = "F:/result/result_11_15/validata/primary/rectal_32_cancer.csv",row.names = F)
#write.csv(methylation,file = "F:/result/12cg/validata/CRC__14_cancer.csv",row.names = F)


methylation<-fread("D:/job-two/new_data/GSE148766_series_matrix.txt.gz",data.table = F,fill=T)
#test<-methylation[1:100,]
methylation<-methylation[c(64:nrow(methylation)),]
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]
#write.csv(methylation,file = "F:/result/result_11_15/validata/primary/CRC_36_cancer.csv",row.names = F)
write.csv(methylation,file = "F:/result/result_11_15/validata/metastatic/CRC_36_cancer.csv",row.names = F)


methylation<-fread("D:/job-two/new_data/GSE149282_series_matrix.txt.gz",data.table = F,fill=T)
test<-methylation[1:100,]
methylation<-methylation[c(70:nrow(methylation)),-grep("normal",test[45,])]
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]
write.csv(methylation,file = "F:/result/result_11_15/validata/primary/CRC_12_cancer.csv",row.names = F)
#write.csv(methylation,file = "F:/result/12cg/validata/CRC__12_cancer.csv",row.names = F)

################## meloma
methylation<-fread("D:/job-two/new_data/GSE133395_series_matrix.txt.gz",data.table = F,fill=T)
#test<-methylation[1:100,]
methylation<-methylation[c(66:nrow(methylation)),]
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]
write.csv(methylation,file = "F:/result/result_11_15/validata/primary/meloma_53_cancer.csv",row.names = F)
#write.csv(methylation,file = "F:/result/12cg/validata/meloma_53_cancer.csv",row.names = F)


methylation<-fread("D:/job-two/new_data/GSE194042_series_matrix.txt.gz",data.table = F,fill=T)
#test<-methylation[1:100,]
methylation<-methylation[c(55:nrow(methylation)),]
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]
write.csv(methylation,file = "F:/result/result_11_15/validata/primary/meloma_30_cancer.csv",row.names = F)
#write.csv(methylation,file = "F:/result/12cg/validata/meloma_30_cancer.csv",row.names = F)


methylation<-fread("D:/job-two/new_data/GSE140171-GPL13534_series_matrix.txt.gz",data.table = F,fill=T)
test<-methylation[1:100,]
methylation<-methylation[c(66:nrow(methylation)),]
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]
write.csv(methylation,file = "F:/result/result_11_15/validata/primary/meloma_46_cancer.csv",row.names = F)
#write.csv(methylation,file = "F:/result/12cg/validata/meloma_46_cancer.csv",row.names = F)


methylation<-fread("D:/job-two/new_data/GSE140171-GPL21145_series_matrix.txt.gz",data.table = F,fill=T)
#test<-methylation[1:100,]
methylation<-methylation[c(64:nrow(methylation)),]
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]
write.csv(methylation,file = "F:/result/result_11_15/validata/primary/meloma_9_cancer.csv",row.names = F)
#write.csv(methylation,file = "F:/result/12cg/validata/meloma_9_cancer.csv",row.names = F)


methylation<-fread("D:/job-two/new_data/GSE144487_Tumor_beta_normalized.txt.gz",data.table = F,fill=T)
#test<-methylation[1:100,]
#methylation<-methylation[c(64:nrow(methylation)),]
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]
write.csv(methylation,file = "F:/result/result_11_15/validata/metastatic/meloma_196_cancer.csv",row.names = F)
#write.csv(methylation,file = "F:/result/12cg/validata/meloma_196_cancer.csv",row.names = F)

#################### Urothelial carcinoma
methylation<-fread("D:/job-two/new_data/GSE161651_utuc_epic_beta_matrix.csv.gz",data.table = F,fill=T)
#test<-methylation[1:100,]
#methylation<-methylation[c(64:nrow(methylation)),]
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]
samples<-fread("D:/job-two/new_data/GSE161651_series_matrix.txt.gz",data.table = F,fill=T)
methylation<-methylation[,-grep("normal",samples[39,])]
write.csv(methylation,file = "F:/result/result_11_15/validata/primary//Urothelial_35_cancer.csv",row.names = F)
#write.csv(methylation,file = "F:/result/12cg/validata/Urothelial_35_cancer.csv",row.names = F)

####################  thyroid
methylation<-fread("D:/job-two/new_data/GSE121377_series_matrix.txt.gz",data.table = F,fill=T)
test<-methylation[1:100,]
methylation<-methylation[c(65:nrow(methylation)),-grep("normal",test[40,])]
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]
write.csv(methylation,file = "F:/result/result_11_15/validata/primary/thyroid_27_cancer.csv",row.names = F)
#write.csv(methylation,file = "F:/result/12cg/validata/thyroid_27_cancer.csv",row.names = F)


methylation<-fread("D:/job-two/new_data/GSE146003_series_matrix.txt.gz",data.table = F,fill=T)
test<-methylation[1:100,]
methylation<-methylation[c(65:nrow(methylation)),-grep("normal",test[41,])]
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]
write.csv(methylation,file = "F:/result/result_11_15/validata/primary/thyroid_10_cancer.csv",row.names = F)
#write.csv(methylation,file = "F:/result/12cg/validata/thyroid_10_cancer.csv",row.names = F)

###################### liver
methylation<-fread("D:/job-two/new_data/GSE132399_series_matrix.txt.gz",data.table = F,fill=T)
test<-methylation[1:150,]
methylation<-methylation[c(113:nrow(methylation)),-c(2:27)]
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]

write.csv(methylation,file = "F:/result/result_11_15/validata/primary/liver_28_cancer.csv",row.names = F)
#write.csv(methylation,file = "F:/result/12cg/validata/liver_28_cancer.csv",row.names = F)

################# mouth cancer

methylation<-fread("D:/job-two/new_data/GSE124633-GPL13534_series_matrix.txt.gz",data.table = F,fill=T)
test<-methylation[1:100,]
methylation<-methylation[c(62:nrow(methylation)),-c(2:6)]
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]
write.csv(methylation,file = "F:/result/result_11_15/validata/primary/mouth_64_cancer.csv",row.names = F)
#write.csv(methylation,file = "F:/result/12cg/validata/mouth_64_cancer.csv",row.names = F)


methylation<-fread("D:/job-two/new_data/GSE124633-GPL21145_series_matrix.txt.gz",data.table = F,fill=T)
test<-methylation[1:100,]
methylation<-methylation[c(62:nrow(methylation)),]
a<-methylation[,1] %in% cg50[,1]
methylation<-methylation[a,]
write.csv(methylation,file = "F:/result/result_11_15/validata/primary/mouth_25_cancer.csv",row.names = F)
#write.csv(methylation,file = "F:/result/12cg/validata/mouth_25_cancer.csv",row.names = F)

###################
library(dplyr)
files<-list.files("F:/result/result_11_15/validata/primary/")
files<-files[grep("csv$",files)]
cg_type<-cg50
a<-c("Adrenal_gland","AML","Bladder","breast","ccRcc","colorectal","CRC","gbm","GBM","HNC","LGG","liver","LUNG","lung_aden","lung_scc","melanoma",
     "mouth","pancreatic","prostate","sarcoma","stomach","thyroid","Urothelial","rectal","meloma")
b<-c("TCGA-ACC","TCGA-LAML","TCGA-BLCA","TCGA-BRCA","TCGA-KIRC","TCGA-COAD","TCGA-COAD","TCGA-GBM","TCGA-GBM","TCGA-HNSC","TCGA-GBM","TCGA-LIHC",
     "TCGA-LUAD","TCGA-LUAD","TCGA-LUSC","TCGA-SKCM","TCGA-HNSC","TCGA-PAAD",
     "TCGA-PRAD","TCGA-SARC","TCGA-STAD","TCGA-THCA","TCGA-BLCA","TCGA-COAD","TCGA-SKCM")
c<-data.frame(a,b)

fin_data<-data.frame()
for(i in 1:length(a)){
  m<-paste0("^",a[i])
  c<-grep(m,files)
  if(i==21){
    next
  } ##### 去除胃癌
  if(i==20){
    next
  } ### 去除肉瘤
  if(i==15){
    next
  } ### 去除肺鳞癌
  if(i==1){
    next
  } ##### 去除胃癌
  # if(i==4){
  #   c<-c[-c(2,3,6,8)]
  # } #### 去除三阴性乳腺癌
  # 
  # if(i==13){
  #   c<-c[-c(1)]
  # } ### 去除不确定的肺腺癌
  # if(i==5){
  #   next
  # } ### 去除ccRcc
  # 
  for(j in 1:length(c)){
    data<-read.csv(paste0("F:/result/result_11_15/validata/primary/",files[c[j]]))
    if(nrow(data)<=47){
      next
    }
    colnames(data)[1]<-colnames(cg_type)[1]
    data<-left_join(cg_type,data,by=colnames(cg_type)[1])
    data<-data[,-2]
    test<-data[,-1]
    test<-apply(test,2,function(x){as.numeric(x)})
    data[,-1]<-test
    #data[is.na(data)]<-0  ###### 缺失值用0代替
    ######### 缺失值用均值代替
    for(h in 2:ncol(data)){
      if(length(which(is.na(data[,h])=="TRUE"))>3){
        next;
      }
      data[is.na(data[,h]),h]<- mean(data[!is.na(data[,h]),h])
    }
    data<-t(data)
    data<-as.data.frame(data)
    colnames(data)<-data[1,]
    data<-data[-1,]
    
    data$CancerTypes<-b[i]
    
    data<-na.omit(data)
    print(files[c[j]])
    print(nrow(data))
    fin_data<-rbind(fin_data,data)
  }
}

#write.csv(data[,-51],file="~/methylation_data/validation_data/test/breast_validata.csv",row.names=F)
#write.csv(data[,51],file="~/methylation_data/validation_data/test/breast_validata_cancertype.csv",row.names=F)

write.csv(fin_data[,-51],file="F:/result/result_11_15/validata/finall/other_method_GEO_validata.csv",row.names=F)
write.csv(fin_data[,51],file="F:/result/result_11_15/validata/finall/other_method_validata_cancertype.csv",row.names=F)

##################
library(dplyr)
files<-list.files("F:/result/result_11_15/validata/metastatic//")
files<-files[grep("csv$",files)]
cg_type<-cg50
a<-c("Adrenal_gland","AML","Bladder","breast","ccRcc","colorectal","CRC","gbm","GBM","HNC","LGG","liver","LUNG","lung_aden","lung_scc","melanoma",
     "mouth","pancreatic","prostate","sarcoma","stomach","thyroid","Urothelial","rectal","meloma")
b<-c("TCGA-ACC","TCGA-LAML","TCGA-BLCA","TCGA-BRCA","TCGA-KIRC","TCGA-COAD","TCGA-COAD","TCGA-GBM","TCGA-GBM","TCGA-HNSC","TCGA-GBM","TCGA-LIHC",
     "TCGA-LUAD","TCGA-LUAD","TCGA-LUSC","TCGA-SKCM","TCGA-HNSC","TCGA-PAAD",
     "TCGA-PRAD","TCGA-SARC","TCGA-STAD","TCGA-THCA","TCGA-BLCA","TCGA-COAD","TCGA-SKCM")
c<-data.frame(a,b)

fin_data<-data.frame()
for(i in 1:length(a)){
  m<-paste0("^",a[i])
  c<-grep(m,files)
  if(length(c)==0){
    next
  }
  if(i==21){
    next
  } ##### 去除胃癌
  if(i==20){
    next
  } ### 去除肉瘤
  if(i==15){
    next
  } ### 去除肺鳞癌
  if(i==1){
    next
  } ##### 去除胃癌
  # if(i==4){
  #   c<-c[-c(2,3,6,8)]
  # } #### 去除三阴性乳腺癌
  # 
  # if(i==13){
  #   c<-c[-c(1)]
  # } ### 去除不确定的肺腺癌
  # if(i==5){
  #   next
  # } ### 去除ccRcc
  # 
  for(j in 1:length(c)){
    data<-read.csv(paste0("F:/result/result_11_15/validata/metastatic/",files[c[j]]))
    if(nrow(data)<=47){
      next
    }
    colnames(data)[1]<-colnames(cg_type)[1]
    data<-left_join(cg_type,data,by=colnames(cg_type)[1])
    data<-data[,-2]
    test<-data[,-1]
    test<-apply(test,2,function(x){as.numeric(x)})
    data[,-1]<-test
    #data[is.na(data)]<-0  ###### 缺失值用0代替
    ######### 缺失值用均值代替
    for(h in 2:ncol(data)){
      if(length(which(is.na(data[,h])=="TRUE"))>3){
        next;
      }
      data[is.na(data[,h]),h]<- mean(data[!is.na(data[,h]),h])
    }
    data<-t(data)
    data<-as.data.frame(data)
    colnames(data)<-data[1,]
    data<-data[-1,]
    
    data$CancerTypes<-b[i]
    
    data<-na.omit(data)
    print(files[c[j]])
    print(nrow(data))
    fin_data<-rbind(fin_data,data)
  }
}

#write.csv(data[,-51],file="~/methylation_data/validation_data/test/breast_validata.csv",row.names=F)
#write.csv(data[,51],file="~/methylation_data/validation_data/test/breast_validata_cancertype.csv",row.names=F)

write.csv(fin_data[,-51],file="F:/result/result_11_15/validata/finall/other_method_GEO_metasta_validata.csv",row.names=F)
write.csv(fin_data[,51],file="F:/result/result_11_15/validata/finall/other_method_metasta_validata_cancertype.csv",row.names=F)

